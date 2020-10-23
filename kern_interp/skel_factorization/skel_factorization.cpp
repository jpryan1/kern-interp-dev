// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include "kern_interp/skel_factorization/skel_factorization.h"
#include "kern_interp/kernel/kernel.h"

#define USING_HIF 0

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

// TODO(John) again this is just redef of above for new type
int SkelFactorization::id_compress(const Kernel& kernel,
                                   const QuadTree* tree,
                                   MidLevelNode* node) {

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
  get_descendents_updates(&update, BN, BN, node, nullptr,
                          nullptr);
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
  // std::cout << "posting schur updates for 1node " << node->center[0] << " " <<
  //           node->center[1] << " " << node->center[2] << " to " << std::endl;
  // for (int i = 0; i < num_skel; i++) {
  //   std::cout << node->dof_lists.skel[i] << ", ";
  // } std::cout << std::endl;
}





void SkelFactorization::decouple(const Kernel& kernel, MidLevelNode* node) {
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
  get_mid_level_schur_updates(&update, BN,  BN, node, nullptr, nullptr,
                              nullptr);

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
  // std::cout << "posting schur updates for 3node " << node->center[0] << " " <<
  //           node->center[1] << " " << node->center[2] << " to " << std::endl;
  // for (int i = 0; i < num_skel; i++) {
  //   std::cout << node->dof_lists.skel[i] << ", ";
  // } std::cout << std::endl;
  node->schur_update = node->L * K_BN_r_sn;
  node->compressed = true;
}

void SkelFactorization::skeletonize(const Kernel& kernel, QuadTree* tree) {
  int node_counter = 0;
  int lvls = tree->levels.size();
  double start, end;
  start = omp_get_wtime();

  int nodes_left = tree->solution_dimension * (kernel.boundary_points_.size() /
                   tree->domain_dimension);
  int prev_nodes_left = nodes_left;
  for (int level = lvls - 1; level >= 1 ; level--) {
    end = omp_get_wtime();
    prev_nodes_left = nodes_left;
    std::cout << "Nodes left " << nodes_left << std::endl;
    start = end;
    tree->remove_inactive_dofs_at_level(level);
    QuadTreeLevel* current_level = tree->levels[level];
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_level->nodes.size(); n++) {
      QuadTreeNode* current_node = current_level->nodes[n];

      if (current_node->compressed || current_node->dof_lists.active_box.size()
          < MIN_DOFS_TO_COMPRESS) {
        if (current_node->compressed) nodes_left -=
            current_node->dof_lists.redundant.size();
        continue;
      }
      double ids = omp_get_wtime();
      if (id_compress(kernel, tree, current_node) == 0) {
        continue;
      }
      double ide = omp_get_wtime();
      nodes_left -= current_node->T.width();
      // std::cout << "Took " << current_node->T.width() << " from node " << std::endl;
      decouple(kernel, current_node);
      node_counter++;
    }
    std::cout << "nodes left after full level " << nodes_left << std::endl;
    if (!USING_HIF) continue;
    tree->populate_half_level_dofs(level);
    MidLevel* current_half_level = tree->levels[level]->half_level;

    int halfdofs = 0;
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_half_level->nodes.size(); n++) {
      MidLevelNode* current_half_level_node = current_half_level->nodes[n];
      halfdofs += current_half_level_node->dof_lists.active_box.size();
      if (current_half_level_node->compressed
          || current_half_level_node->dof_lists.active_box.size() <
          MIN_DOFS_TO_COMPRESS / 4) {
        if (current_half_level_node->compressed) nodes_left -=
            current_half_level_node->dof_lists.redundant.size();

        continue;
      }
      // std::cout<<"Tried recompression on halfbox with "<<current_half_level_node->dof_lists.active_box.size() <<" dofs"<<std::endl;
      if (id_compress(kernel, tree, current_half_level_node) == 0) {
        // std::cout<<"failed"<<std::endl;
        continue;
      }        
      // std::cout<<"Successfully removed "<<current_half_level_node->T.width()<<std::endl;
      nodes_left -= current_half_level_node->T.width();

      decouple(kernel, current_half_level_node);

      node_counter++;
    }
    std::cout << "nodes left after halflevel " << nodes_left << " found " <<
              halfdofs << " to compress" << std::endl;

    if (kernel.domain_dimension == 3) {
      tree->populate_third_level_dofs(level);
      MidLevel* current_third_level = tree->levels[level]->third_level;
      int thirddofs = 0;
      #pragma omp parallel for num_threads(fact_threads)
      for (int n = 0; n < current_third_level->nodes.size(); n++) {
        MidLevelNode* current_third_level_node = current_third_level->nodes[n];
        thirddofs += current_third_level_node->dof_lists.active_box.size();
        if (current_third_level_node->compressed
            || current_third_level_node->dof_lists.active_box.size() <
            MIN_DOFS_TO_COMPRESS / 4) {
          continue;
        }
        if (id_compress(kernel, tree, current_third_level_node) == 0) {
          continue;
        }
        nodes_left -= current_third_level_node->T.width();
        decouple(kernel, current_third_level_node);
        node_counter++;
      }
      std::cout << "nodes left after thirdlevel " << nodes_left << " found " <<
                thirddofs << " to compress" << std::endl;
    }
  }
  end = omp_get_wtime();
  // std::cout << "Last level " << (end - start) << std::endl;
  // std::cout << "Nodes left " << nodes_left << std::endl;
  // If the above breaks due to a cap, we need to manually propagate active
  // boxes up the tree.
  tree->remove_inactive_dofs_at_all_boxes();
  std::vector<int> allskel = tree->root->dof_lists.active_box;

  if (allskel.size() > 0) {
    ki_Mat allskel_updates = ki_Mat(allskel.size(), allskel.size());
    get_descendents_updates(&allskel_updates, allskel, allskel, tree->root, nullptr,
                            nullptr);
    start = omp_get_wtime();
    tree->allskel_mat = kernel(allskel, allskel, false, true) - allskel_updates;
    end = omp_get_wtime();
    // std::cout << "allskel get time " << end - start << std::endl;
  }
  std::cout << "Allskelsize " << allskel.size() << std::endl;
  if (tree->U.width() == 0) {
    double lufs = omp_get_wtime();
    openblas_set_num_threads(fact_threads);
    tree->allskel_mat.LU_factorize(&tree->allskel_mat_lu,
                                   &tree->allskel_mat_piv);
    openblas_set_num_threads(1);
    double lufe = omp_get_wtime();
    std::cout << "allskel lu done"<< std::endl;
    return;
  }
// std::cout<<"all skel cond "<<tree->allskel_mat.condition_number()<<std::endl;
  std::vector<QuadTreeNode*> all_nodes;
  std::vector<MidLevelNode*> all_half_nodes;
  std::vector<MidLevelNode*> all_third_nodes;
  for (int level = lvls - 1; level >= 0; level--) {
    QuadTreeLevel* current_level = tree->levels[level];
    for (int n = 0; n < current_level->nodes.size(); n++) {
      all_nodes.push_back(current_level->nodes[n]);
    }
    if (!USING_HIF) continue;
    for (int n = 0; n < current_level->half_level->nodes.size(); n++) {
      all_half_nodes.push_back(current_level->half_level->nodes[n]);
    }
    if (kernel.domain_dimension == 3) {
      for (int n = 0; n < current_level->third_level->nodes.size(); n++) {
        all_third_nodes.push_back(current_level->third_level->nodes[n]);
      }
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

    if (!USING_HIF) continue;
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_level->half_level->nodes.size(); n++) {
      MidLevelNode* current_node = current_level->half_level->nodes[n];
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

    if (kernel.domain_dimension == 3) {

      #pragma omp parallel for num_threads(fact_threads)
      for (int n = 0; n < current_level->third_level->nodes.size(); n++) {
        MidLevelNode* current_node = current_level->third_level->nodes[n];
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
    if (!USING_HIF) continue;
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_level->half_level->nodes.size(); n++) {
      MidLevelNode* current_node = current_level->half_level->nodes[n];
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
    if (kernel.domain_dimension == 3) {

      #pragma omp parallel for num_threads(fact_threads)
      for (int n = 0; n < current_level->third_level->nodes.size(); n++) {
        MidLevelNode* current_node = current_level->third_level->nodes[n];
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
  if (USING_HIF) {
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < all_half_nodes.size(); n++) {
      MidLevelNode* current_node = all_half_nodes[n];

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
    if (kernel.domain_dimension == 3) {

      #pragma omp parallel for num_threads(fact_threads)
      for (int n = 0; n < all_third_nodes.size(); n++) {
        MidLevelNode* current_node = all_third_nodes[n];

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
    }
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
  openblas_set_num_threads(fact_threads);
  S.LU_factorize(&tree->S_LU, &tree->S_piv);
  openblas_set_num_threads(1);
  double sluend = omp_get_wtime();
  // std::cout << "S LU factorization: " << (sluend - slustart) << std::endl;

}


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
  std::vector<MidLevelNode*> all_half_nodes;
  std::vector<MidLevelNode*> all_third_nodes;

  if (USING_HIF) {
    for (int level = lvls - 1; level >= 0; level--) {
      MidLevel* current_level = quadtree.levels[level]->half_level;
      for (int n = 0; n < current_level->nodes.size(); n++) {
        all_half_nodes.push_back(current_level->nodes[n]);
      }
    }
    if (quadtree.domain_dimension == 3) {

      for (int level = lvls - 1; level >= 0; level--) {
        MidLevel* current_level = quadtree.levels[level]->third_level;
        for (int n = 0; n < current_level->nodes.size(); n++) {
          all_third_nodes.push_back(current_level->nodes[n]);
        }
      }
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
    if (!USING_HIF) continue;
    MidLevel* current_half_level = quadtree.levels[level]->half_level;
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_half_level->nodes.size(); n++) {
      MidLevelNode* current_node = current_half_level->nodes[n];
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
    if (quadtree.domain_dimension == 3) {

      MidLevel* current_third_level = quadtree.levels[level]->third_level;
      #pragma omp parallel for num_threads(fact_threads)
      for (int n = 0; n < current_third_level->nodes.size(); n++) {
        MidLevelNode* current_node = current_third_level->nodes[n];
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
  if (USING_HIF) {
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < all_half_nodes.size(); n++) {
      MidLevelNode* current_node = all_half_nodes[n];
      if (current_node->dof_lists.redundant.size() == 0) continue;
      if (!current_node->compressed) {
        continue;
      }
      assert(current_node->X_rr_is_LU_factored);

      apply_diag_inv_matrix(current_node->X_rr_lu, current_node->X_rr_piv, x,
                            current_node->dof_lists.redundant);
    }
    if (quadtree.domain_dimension == 3) {

      #pragma omp parallel for num_threads(fact_threads)
      for (int n = 0; n < all_third_nodes.size(); n++) {
        MidLevelNode* current_node = all_third_nodes[n];
        if (current_node->dof_lists.redundant.size() == 0) continue;
        if (!current_node->compressed) {
          continue;
        }
        assert(current_node->X_rr_is_LU_factored);

        apply_diag_inv_matrix(current_node->X_rr_lu, current_node->X_rr_piv, x,
                              current_node->dof_lists.redundant);
      }
    }

  }
  std::vector<int> allskel = quadtree.root->dof_lists.active_box;
  if (allskel.size() > 0) {
    apply_diag_inv_matrix(quadtree.allskel_mat_lu, quadtree.allskel_mat_piv, x,
                          allskel);
  }


  for (int level = 0; level < lvls; level++) {
    if (quadtree.domain_dimension == 3) {
      if (USING_HIF) {
        MidLevel* current_third_level = quadtree.levels[level]->third_level;
        #pragma omp parallel for num_threads(fact_threads)
        for (int n = current_third_level->nodes.size() - 1; n >= 0; n--) {
          MidLevelNode* current_node = current_third_level->nodes[n];
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
      MidLevel* current_half_level = quadtree.levels[level]->half_level;
      #pragma omp parallel for num_threads(fact_threads)
      for (int n = current_half_level->nodes.size() - 1; n >= 0; n--) {
        MidLevelNode* current_node = current_half_level->nodes[n];
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


  std::vector<MidLevelNode*> all_half_nodes;
  std::vector<MidLevelNode*> all_third_nodes;

  if (USING_HIF) {
    for (int level = lvls - 1; level >= 0; level--) {
      MidLevel* current_half_level = quadtree.levels[level]->half_level;
      for (int n = 0; n < current_half_level->nodes.size(); n++) {
        all_half_nodes.push_back(current_half_level->nodes[n]);
      }
    }


    if (quadtree.domain_dimension == 3) {

      for (int level = lvls - 1; level >= 0; level--) {
        MidLevel* current_third_level = quadtree.levels[level]->third_level;
        for (int n = 0; n < current_third_level->nodes.size(); n++) {
          all_third_nodes.push_back(current_third_level->nodes[n]);
        }
      }
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
    if (!USING_HIF) continue;

    MidLevel* current_half_level = quadtree.levels[level]->half_level;

    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_half_level->nodes.size(); n++) {
      MidLevelNode* current_node = current_half_level->nodes[n];
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
    if (quadtree.domain_dimension == 3) {

      MidLevel* current_third_level = quadtree.levels[level]->third_level;

      #pragma omp parallel for num_threads(fact_threads)
      for (int n = 0; n < current_third_level->nodes.size(); n++) {
        MidLevelNode* current_node = current_third_level->nodes[n];
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
    if (!USING_HIF) continue;
    MidLevel* current_half_level = quadtree.levels[level]->half_level;
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < current_half_level->nodes.size(); n++) {
      MidLevelNode* current_node = current_half_level->nodes[n];
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
    if (quadtree.domain_dimension == 3) {

      MidLevel* current_third_level = quadtree.levels[level]->third_level;
      #pragma omp parallel for num_threads(fact_threads)
      for (int n = 0; n < current_third_level->nodes.size(); n++) {
        MidLevelNode* current_node = current_third_level->nodes[n];
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
  if (USING_HIF) {
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < all_half_nodes.size(); n++) {
      MidLevelNode* current_node = all_half_nodes[n];
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
    if (quadtree.domain_dimension == 3) {

      #pragma omp parallel for num_threads(fact_threads)
      for (int n = 0; n < all_third_nodes.size(); n++) {
        MidLevelNode* current_node = all_third_nodes[n];
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
    }
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
  if (USING_HIF) {
    #pragma omp parallel for num_threads(fact_threads)
    for (int n = 0; n < all_half_nodes.size(); n++) {
      MidLevelNode* current_node = all_half_nodes[n];
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

    if (quadtree.domain_dimension == 3) {

      #pragma omp parallel for num_threads(fact_threads)
      for (int n = 0; n < all_third_nodes.size(); n++) {
        MidLevelNode* current_node = all_third_nodes[n];
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
    }
  }
  mu->set_submatrix(allredundant, 0, 1, Dinv_N);
  mu->set_submatrix(allskel, 0, 1, y(0, allskel.size(), 0, 1));

  for (int level = 0; level < lvls; level++) {
    if (USING_HIF) {
      if (quadtree.domain_dimension == 3) {

        MidLevel* current_third_level = quadtree.levels[level]->third_level;
        #pragma omp parallel for num_threads(fact_threads)
        for (int n = current_third_level->nodes.size() - 1; n >= 0; n--) {
          MidLevelNode* current_node = current_third_level->nodes[n];
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
      MidLevel* current_half_level = quadtree.levels[level]->half_level;
      #pragma omp parallel for num_threads(fact_threads)
      for (int n = current_half_level->nodes.size() - 1; n >= 0; n--) {
        MidLevelNode* current_node = current_half_level->nodes[n];
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

