// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include "kern_interp/skel_factorization/skel_factorization.h"
#include "kern_interp/kernel/kernel.h"

namespace kern_interp {


std::vector<int> big_to_small(const std::vector<int>& big,
                              const std::unordered_map<int,
                              int>& map) {
  std::vector<int> small;
  for (int idx : big) {
    small.push_back(map.at(idx));
  }
  return small;
}


SkelFactorization::SkelFactorization(double id_tol, int fact_threads) {
  assert(id_tol > 0 && "id_tol must be greater than one to init tools.");
  this->id_tol = id_tol;
  this->fact_threads = fact_threads;
}


int SkelFactorization::id_compress(const Kernel& kernel,
                                   const QuadTree* tree, QuadTreeNode* node) {

  assert(node != nullptr && "InterpolativeDecomposition fails on null node.");
  assert(node->dof_lists.active_box.size() > 0 &&
         "Num of DOFs must be positive in InterpolativeDecomposition.");

  ki_Mat pxy = kernel.get_id_mat(tree, node);
  if (pxy.height() == 0) {
    return 0;
  }
  std::vector<int> p;
  int numskel = pxy.id(&p, &node->T, id_tol);
  if (numskel == 0) {
    return 0;
  }
  node->dof_lists.set_rs_ranges(p, node->T.height(), node->T.width());
  node->dof_lists.set_skelnear_range();

  return node->T.width();
}


void SkelFactorization::decouple(const Kernel& kernel, QuadTreeNode* node) {

  // height of Z is number of skeleton columns
  int num_redundant = node->T.width();
  int num_skel      = node->T.height();
  // GENERATE K_BN,BN
  std::vector<int> BN;
  for (int idx : node->dof_lists.active_box) {
    BN.push_back(idx);
  }

  // Note that BN has all currently deactivated DoFs removed.
  ki_Mat update(BN.size(), BN.size());
  get_all_schur_updates(&update, BN, node);
  ki_Mat K_BN = kernel(BN, BN) - update;

  // Generate various index ranges within BN
  std::vector<int> s, r, n, sn;
  for (int i = 0; i < num_skel; i++) {
    s.push_back(node->dof_lists.permutation[i]);
    sn.push_back(node->dof_lists.permutation[i]);
  }

  for (int i = 0; i < num_redundant; i++) {
    r.push_back(node->dof_lists.permutation[i + num_skel]);
  }
  ki_Mat K_BN_r_sn = K_BN(r, s) - node->T.transpose() * K_BN(s, s);
  node->X_rr = K_BN(r, r) - node->T.transpose() * K_BN(s, r)
               - K_BN_r_sn * node->T;
  ki_Mat K_BN_sn_r = K_BN(s, r) - K_BN(s, s) * node->T;

  node->X_rr.LU_factorize(&node->X_rr_lu, &node->X_rr_piv);

  node->X_rr_is_LU_factored = true;
  node->L = ki_Mat(sn.size(),  num_redundant);
  node->U = ki_Mat(num_redundant, sn.size());
  node->X_rr_lu.right_multiply_inverse(K_BN_sn_r, node->X_rr_piv, &node->L);
  node->X_rr_lu.left_multiply_inverse(K_BN_r_sn, node->X_rr_piv,  &node->U);
  node->schur_update = node->L * K_BN_r_sn;
  node->compressed = true;

}


void SkelFactorization::skeletonize(const Kernel& kernel, QuadTree* tree) {
  int node_counter = 0;
  int lvls = tree->levels.size();

  double start, end;
  start = omp_get_wtime();
  int nodes_left = kernel.boundary_points_.size();
  // std::cout << "Starting with " << nodes_left << " nodes " << std::endl;
  int prev_nodes_left = nodes_left;
  for (int level = lvls - 1; level >  1; level--) {
    end = omp_get_wtime();
     std::cout << "level " << level << " " << (end - start) << std::endl;
    prev_nodes_left = nodes_left;
    start = end;
    tree->remove_inactive_dofs_at_level(level);
    QuadTreeLevel* current_level = tree->levels[level];
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];
      if (current_node->compressed || current_node->dof_lists.active_box.size()
          < MIN_DOFS_TO_COMPRESS) {
        continue;
      }
      double ids = omp_get_wtime();
      if (id_compress(kernel, tree, current_node) == 0) {
        continue;
      }
      double ide = omp_get_wtime();
      // std::cout << "took " << ide - ids << " to get rid of " <<
      //           current_node->T.width() << " out of "<<(current_node->T.width()+current_node->T.height())<<" center "<<current_node->center[0]<<" "<<current_node->center[1]<<" "<<current_node->center[2]<<std::endl;
      nodes_left -= current_node->T.width();
      decouple(kernel, current_node);
      node_counter++;
    }
  }

  end = omp_get_wtime();
   std::cout << "Last level " << (end - start) << std::endl;
  // std::cout << "Nodes left " << nodes_left << std::endl;
  // If the above breaks due to a cap, we need to manually propagate active
  // boxes up the tree.
  tree->remove_inactive_dofs_at_all_boxes();
  std::vector<int> allskel = tree->root->dof_lists.active_box;

  if (allskel.size() > 0) {
    ki_Mat allskel_updates = ki_Mat(allskel.size(), allskel.size());
    get_all_schur_updates(&allskel_updates, allskel, tree->root);
    start = omp_get_wtime();
    tree->allskel_mat = kernel(allskel, allskel, false, true) - allskel_updates;
    end = omp_get_wtime();
     std::cout << "allskel get time " << end - start << std::endl;
  }

  if (tree->U.width() == 0) {
    double lufs = omp_get_wtime();
    openblas_set_num_threads(fact_threads);
    tree->allskel_mat.LU_factorize(&tree->allskel_mat_lu,
                                   &tree->allskel_mat_piv);
    openblas_set_num_threads(1);
    double lufe = omp_get_wtime();
    // std::cout << "allskel lu " << (lufe - lufs) << std::endl;
    return;
  }
// std::cout<<"all skel cond "<<tree->allskel_mat.condition_number()<<std::endl;
  std::vector<QuadTreeNode*> all_nodes;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = tree->levels[level];
    for (int n = 0; n < current_level->nodes.size(); n++) {
      all_nodes.push_back(current_level->nodes[n]);
    }
  }

  std::vector<int> sorted_allskel = allskel;
  std::sort(sorted_allskel.begin(), sorted_allskel.end());
  int skel_idx = 0;

  std::vector<int> allredundant;
  for (int i = 0;
       i < (kernel.boundary_points_.size() / kernel.domain_dimension)
       *kernel.solution_dimension;
       i++) {
    if (skel_idx < sorted_allskel.size() && i == sorted_allskel[skel_idx]) {
      skel_idx++;
    } else {
      allredundant.push_back(i);
    }
  }
  if (allredundant.size() == 0) {
    std::cout << "No compression possible" << std::endl;
    exit(0);
  }

  // In our bordered linear system, the skel and redundant indices are
  // partitioned so we create a map from their original index into their
  // partition
  std::unordered_map<int, int> skel_big2small, red_big2small;
  for (int i = 0; i < allskel.size(); i++) {
    skel_big2small[allskel[i]] = i;
  }
  for (int i = 0; i < allredundant.size(); i++) {
    red_big2small[allredundant[i]] = i;
  }
  ki_Mat modified_Psi = tree->Psi.transpose();
  ki_Mat modified_U = tree->U;

  // First apply the sweep matrices to x and U to modify them.
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = tree->levels[level];

    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->compressed) {
        continue;
      }
      apply_sweep_matrix(-current_node->T, &modified_U,
                         current_node->dof_lists.skel,
                         current_node->dof_lists.redundant, true);
      apply_sweep_matrix(-current_node->L, &modified_U,
                         current_node->dof_lists.redundant,
                         current_node->dof_lists.skelnear, false);
    }
  }

  // Now apply the other sweep matrices to Psi to modify it.
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = tree->levels[level];
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->compressed) {
        continue;
      }
      apply_sweep_matrix(-current_node->T, &modified_Psi,
                         current_node->dof_lists.skel,
                         current_node->dof_lists.redundant,
                         true);
      apply_sweep_matrix(-current_node->U, &modified_Psi,
                         current_node->dof_lists.redundant,
                         current_node->dof_lists.skelnear,
                         true);
    }
  }

  modified_Psi = modified_Psi.transpose();
  // Again, C is mostly 0s, so we just apply Dinv to the nonzero block
  ki_Mat Dinv_C_nonzero = modified_U(allredundant, 0, modified_U.width());
  #pragma omp parallel for num_threads(fact_threads)
  for (int n = 0; n < all_nodes.size(); n++) {
    QuadTreeNode* current_node = all_nodes[n];

    if (current_node->dof_lists.redundant.size() == 0) continue;
    if (!current_node->compressed) {
      continue;
    }
    std::vector<int> small_redundants = big_to_small(
                                          current_node->dof_lists.redundant,
                                          red_big2small);
    assert(current_node->X_rr_is_LU_factored);
    apply_diag_inv_matrix(current_node->X_rr_lu, current_node->X_rr_piv,
                          &Dinv_C_nonzero,
                          small_redundants);
  }
  ki_Mat ident(tree->Psi.height(), tree->Psi.height());
  if (kernel.domain_dimension == 2) {
    ident.eye(tree->Psi.height());
  }
  ki_Mat S(allskel.size() + tree->Psi.height(),
           allskel.size() + tree->Psi.height());

  S.set_submatrix(0, allskel.size(),
                  0, allskel.size(), tree->allskel_mat);
  S.set_submatrix(allskel.size(), S.height(),
                  0, allskel.size(),
                  modified_Psi(0, tree->Psi.height(), allskel));
  S.set_submatrix(0, allskel.size(),
                  allskel.size(), S.width(),
                  modified_U(allskel, 0, tree->U.width()));

  S.set_submatrix(allskel.size(), S.height(), allskel.size(), S.width(),
                  - ident - (modified_Psi(0, modified_Psi.height(),
                                          allredundant) * Dinv_C_nonzero));
  double slustart = omp_get_wtime();
  // std::cout << "S height " << S.height() << std::endl;
  // std::cout << "slustart " << std::endl;
  openblas_set_num_threads(fact_threads);
  S.LU_factorize(&tree->S_LU, &tree->S_piv);
  openblas_set_num_threads(1);
  double sluend = omp_get_wtime();
  // std::cout << "slu " << (sluend - slustart) << std::endl;

}


void SkelFactorization::get_all_schur_updates(ki_Mat * updates,
    const std::vector<int>& BN, const QuadTreeNode * node) const {
  assert(node != nullptr && "get_all_schur_updates fails on null node.");
  assert(BN.size() > 0 && "get_all_schur_updates needs positive num of DOFs");
  if (!node->is_leaf) get_descendents_updates(updates, BN, node);
}


void SkelFactorization::get_descendents_updates(ki_Mat * updates,
    const std::vector<int>& BN, const QuadTreeNode * node)  const {
  assert(node != nullptr && "get_descendents_updates fails on null node.");
  assert(!node->is_leaf &&
         "get_descendents_updates must be called on non-leaf.");
  // by assumption, node is not a leaf
  for (QuadTreeNode* child : node->children) {
    if (child->compressed) get_update(updates, BN, child);
    if (!child->is_leaf) get_descendents_updates(updates, BN, child);
  }
}


void SkelFactorization::get_update(ki_Mat * update,
                                   const std::vector<int>& BN,
                                   const QuadTreeNode * node)  const {
  // node needs to check all its dofs against BN, enter interactions into
  // corresponding locations
  // node only updated its own BN dofs, and the redundant ones are no longer
  // relevant, so we only care about child's SN dofs
  // First create a list of Dofs that are also in node's skelnear,
  // and with each one give the index in skelnear and the index in BN
  std::vector<int> BN_;
  std::vector<int> sn_;
  for (int sn_idx = 0; sn_idx < node->dof_lists.skelnear.size();
       sn_idx++) {
    for (int bn_idx = 0; bn_idx < BN.size(); bn_idx++) {
      if (BN[bn_idx] == node->dof_lists.skelnear[sn_idx]) {
        sn_.push_back(sn_idx);
        BN_.push_back(bn_idx);
      }
    }
  }
  // For every pair of dofs shared by both, update their interaction
  int num_shared_by_both = BN_.size();
  for (int i = 0; i < num_shared_by_both; i++) {
    for (int j = 0; j < num_shared_by_both; j++) {
      update->addset(BN_[i], BN_[j], node->schur_update.get(sn_[i], sn_[j]));
    }
  }
}


///////////////////////////////////////////////////////////////////////////////


// Sets vec(b) = vec(b) + mat*vec(a)
void SkelFactorization::apply_sweep_matrix(const ki_Mat & mat, ki_Mat * vec,
    const std::vector<int>& a,
    const std::vector<int>& b,
    bool transpose = false) const {
  if (a.size()*b.size() == 0) return;
  if (transpose) {
    assert(mat.height() == a.size());
  } else {
    assert(mat.width() == a.size());
  }
  ki_Mat product;
  if (transpose) {
    product = mat.transpose() * (*vec)(a, 0, vec->width());
  } else {
    product = mat * (*vec)(a, 0, vec->width());
  }
  vec->set_submatrix(b, 0, vec->width(), product + (*vec)(b, 0, vec->width()));
}


// Sets vec(range) = mat * vec(range)
void SkelFactorization::apply_diag_matrix(const ki_Mat & mat, ki_Mat * vec,
    const std::vector<int>& range)
const {
  if (range.size() == 0) return;
  vec->set_submatrix(range,  0, vec->width(),  mat * (*vec)(range, 0,
                     vec->width()));
}


void SkelFactorization::apply_diag_inv_matrix(const ki_Mat & mat,
    const std::vector<lapack_int>& piv, ki_Mat * vec,
    const std::vector<int>& range) const {
  if (range.size() == 0) return;
  ki_Mat product(range.size(),  vec->width());
  mat.left_multiply_inverse((*vec)(range,  0, vec->width()), piv, &product);
  vec->set_submatrix(range,  0, vec->width(), product);
}


void SkelFactorization::solve(const QuadTree & quadtree, ki_Mat * x,
                              const ki_Mat & b) const {
  assert(x->height() == b.height());
  int lvls = quadtree.levels.size();
  *x = b;
  std::vector<QuadTreeNode*> all_nodes;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (int n = 0; n < current_level->nodes.size(); n++) {
      all_nodes.push_back(current_level->nodes[n]);
    }
  }
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->compressed) {
        continue;
      }
      apply_sweep_matrix(-current_node->T, x,
                         current_node->dof_lists.skel,
                         current_node->dof_lists.redundant, true);
      apply_sweep_matrix(-current_node->L, x,
                         current_node->dof_lists.redundant,
                         current_node->dof_lists.skelnear, false);
    }
  }
  #pragma omp parallel for num_threads(fact_threads)
  for (int n = 0; n < all_nodes.size(); n++) {
    QuadTreeNode* current_node = all_nodes[n];
    if (current_node->dof_lists.redundant.size() == 0) continue;
    if (!current_node->compressed) {
      continue;
    }
    assert(current_node->X_rr_is_LU_factored);

    apply_diag_inv_matrix(current_node->X_rr_lu, current_node->X_rr_piv, x,
                          current_node->dof_lists.redundant);
  }
  std::vector<int> allskel = quadtree.root->dof_lists.active_box;
  if (allskel.size() > 0) {
    apply_diag_inv_matrix(quadtree.allskel_mat_lu, quadtree.allskel_mat_piv, x,
                          allskel);
  }
  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->compressed) {
        continue;
      }
      apply_sweep_matrix(-current_node->U, x,
                         current_node->dof_lists.skelnear,
                         current_node->dof_lists.redundant, false);
      apply_sweep_matrix(-current_node->T, x,
                         current_node->dof_lists.redundant,
                         current_node->dof_lists.skel,
                         false);
    }
  }
}


void SkelFactorization::multiply_connected_solve(const QuadTree & quadtree,
    ki_Mat * mu, ki_Mat * alpha, const ki_Mat & b) const {

  assert(mu->height() == b.height());
  int lvls = quadtree.levels.size();
  std::vector<QuadTreeNode*> all_nodes;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    for (int n = 0; n < current_level->nodes.size(); n++) {
      all_nodes.push_back(current_level->nodes[n]);
    }
  }
  std::vector<int> allskel = quadtree.root->dof_lists.active_box;
  std::vector<int> sorted_allskel = allskel;
  std::sort(sorted_allskel.begin(), sorted_allskel.end());
  int skel_idx = 0;
  std::vector<int> allredundant;
  for (int i = 0; i < b.height(); i++) {
    if (skel_idx < sorted_allskel.size() && i == sorted_allskel[skel_idx]) {
      skel_idx++;
    } else {
      allredundant.push_back(i);
    }
  }
  // In our bordered linear system, the skel and redundant indices are
  // partitioned so we create a map from their original index into their
  // partition
  std::unordered_map<int, int> skel_big2small, red_big2small;
  for (int i = 0; i < allskel.size(); i++) {
    skel_big2small[allskel[i]] = i;
  }
  for (int i = 0; i < allredundant.size(); i++) {
    red_big2small[allredundant[i]] = i;
  }

  *mu = b;
  ki_Mat modified_Psi = quadtree.Psi.transpose();
  ki_Mat modified_U = quadtree.U;
  // First apply the sweep matrices to x and U to modify them.
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];

    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->compressed) {
        continue;
      }
      apply_sweep_matrix(-current_node->T, mu,
                         current_node->dof_lists.skel,
                         current_node->dof_lists.redundant, true);
      apply_sweep_matrix(-current_node->L, mu,
                         current_node->dof_lists.redundant,
                         current_node->dof_lists.skelnear, false);
      apply_sweep_matrix(-current_node->T, &modified_U,
                         current_node->dof_lists.skel,
                         current_node->dof_lists.redundant, true);
      apply_sweep_matrix(-current_node->L, &modified_U,
                         current_node->dof_lists.redundant,
                         current_node->dof_lists.skelnear, false);
    }
  }
  // After the result of the first sweep matrices, grab w and z.
  ki_Mat w = (*mu)(allredundant, 0, 1);
  ki_Mat z = (*mu)(allskel, 0, 1);
  // Now apply the other sweep matrices to Psi to modify it.
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->compressed) {
        continue;
      }
      apply_sweep_matrix(-current_node->T, &modified_Psi,
                         current_node->dof_lists.skel,
                         current_node->dof_lists.redundant,
                         true);
      apply_sweep_matrix(-current_node->U, &modified_Psi,
                         current_node->dof_lists.redundant,
                         current_node->dof_lists.skelnear,
                         true);
    }
  }
  modified_Psi = modified_Psi.transpose();
  ki_Mat Dinv_w = w;
  #pragma omp parallel for num_threads(fact_threads)
  for (int n = 0; n < all_nodes.size(); n++) {
    QuadTreeNode* current_node = all_nodes[n];
    if (current_node->dof_lists.redundant.size() == 0) continue;
    if (!current_node->compressed) {
      continue;
    }
    std::vector<int> small_redundants = big_to_small(
                                          current_node->dof_lists.redundant,
                                          red_big2small);
    assert(current_node->X_rr_is_LU_factored);
    apply_diag_inv_matrix(current_node->X_rr_lu, current_node->X_rr_piv,
                          &Dinv_w, small_redundants);
  }
  ki_Mat M(allskel.size() + quadtree.Psi.height(), 1);
  M.set_submatrix(0, allskel.size(), 0, 1, z);
  M.set_submatrix(allskel.size(), M.height(), 0, 1, -(modified_Psi(0,
                  modified_Psi.height(),
                  allredundant) * Dinv_w));
  ki_Mat y(quadtree.S_LU.height(), 1);
  quadtree.S_LU.left_multiply_inverse(M, quadtree.S_piv, &y);
  *alpha =  y(allskel.size(), y.height(), 0, 1);
  ki_Mat N = w - modified_U(allredundant, 0, modified_U.width()) * (*alpha);
  ki_Mat Dinv_N = N;
  #pragma omp parallel for num_threads(fact_threads)
  for (int n = 0; n < all_nodes.size(); n++) {
    QuadTreeNode* current_node = all_nodes[n];
    if (current_node->dof_lists.redundant.size() == 0) continue;
    if (!current_node->compressed) {
      continue;
    }
    std::vector<int> small_redundants = big_to_small(
                                          current_node->dof_lists.redundant,
                                          red_big2small);
    assert(current_node->X_rr_is_LU_factored);

    apply_diag_inv_matrix(current_node->X_rr_lu, current_node->X_rr_piv,
                          &Dinv_N, small_redundants);
  }
  mu->set_submatrix(allredundant, 0, 1, Dinv_N);
  mu->set_submatrix(allskel, 0, 1, y(0, allskel.size(), 0, 1));

  for (int level = 0; level < lvls; level++) {
    QuadTreeLevel* current_level = quadtree.levels[level];
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = current_level->nodes.size() - 1; n >= 0; n--) {
      QuadTreeNode* current_node = current_level->nodes[n];
      if (!current_node->compressed) {
        continue;
      }
      apply_sweep_matrix(-current_node->U, mu,
                         current_node->dof_lists.skelnear,
                         current_node->dof_lists.redundant, false);
      apply_sweep_matrix(-current_node->T, mu,
                         current_node->dof_lists.redundant,
                         current_node->dof_lists.skel,
                         false);
    }
  }
}


}  // namespace kern_interp

