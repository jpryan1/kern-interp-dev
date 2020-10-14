// Copyright 2019 John Paul Ryan
#ifndef KERN_INTERP_QUADTREE_QUADTREE_H_
#define KERN_INTERP_QUADTREE_QUADTREE_H_

#include <cassert>
#include <set>
#include <vector>
#include "kern_interp/ki_mat.h"
#include "kern_interp/boundaries/boundary.h"

#define MAX_LEAF_DOFS 128

namespace kern_interp {

struct InteractionLists {
  std::vector<int> original_box,
      active_box,
      redundant,
      skel,
      near,
      skelnear,
      permutation;

  void set_rs_ranges(const std::vector<int>& prm, int sk,
                     int rd) {
    assert(prm.size() == sk + rd);

    for (int i = 0; i < sk; i++) {
      skel.push_back(active_box[prm[i]]);
      permutation.push_back(prm[i]);
    }
    for (int i = sk; i < sk + rd; i++) {
      redundant.push_back(active_box[prm[i]]);
      permutation.push_back(prm[i]);
    }
  }

  void set_skelnear_range() {
    for (int i = 0; i < skel.size(); i++) {
      skelnear.push_back(skel[i]);
    }
  }
};

struct QuadTreeLevel;
struct HalfLevelNode;

struct QuadTreeNode {
  int level, dofs_below;
  bool is_leaf,
       X_rr_is_LU_factored = false,
       compressed = false,
       p_marked = false;
  double side_length;
  QuadTreeNode* parent;
  std::vector<QuadTreeNode*> children;
  std::vector<QuadTreeNode*> neighbors;

  InteractionLists dof_lists;
  // For inverse operator
  ki_Mat T, L, U, X_rr, schur_update, X_rr_lu;

  std::vector<lapack_int> X_rr_piv;

  std::vector<HalfLevelNode*> recompressor_nodes;

  std::vector<double> center;
  QuadTreeNode() {
    is_leaf = true;
  }
};


struct HalfLevelNode {
  int partner_level;
  bool X_rr_is_LU_factored = false, compressed = false;

  double side_length;
  std::vector<QuadTreeNode*> containing_nodes;

  InteractionLists dof_lists;
  // For inverse operator
  ki_Mat T, L, U, X_rr, schur_update, X_rr_lu;
  std::vector<lapack_int> X_rr_piv;
  std::vector<double> center;
  HalfLevelNode() {}
};


struct HalfLevel {
  std::vector<HalfLevelNode*> nodes;
  HalfLevel* half_level;
  QuadTreeLevel* partner_level;
  ~HalfLevel() {
    for (HalfLevelNode* node : nodes) {
      delete node;
    }
  }
};


struct QuadTreeLevel {
  std::vector<QuadTreeNode*> nodes;
  HalfLevel* half_level;
  ~QuadTreeLevel() {
    for (QuadTreeNode* node : nodes) {
      delete node;
    }
  }
};


class QuadTree {
 public:
  int solution_dimension, domain_dimension;
  int no_proxy_level = 0;
  double min, max;
  std::vector<double> boundary_points;
  QuadTreeNode* root;
  ki_Mat allskel_mat, allskel_mat_lu, U, Psi, S_LU;
  std::vector<lapack_int> allskel_mat_piv, S_piv;
  std::vector<QuadTreeLevel*> levels;
  ~QuadTree();
  void initialize_tree(Boundary* boundary, int solution_dimension_,
                       int domain_dimension_);
  void compute_neighbor_lists();
  void compute_half_levels();
  void recursive_add(QuadTreeNode* node, std::vector<double> pt,
                     int mat_ind);
  void get_descendent_neighbors(QuadTreeNode* big, QuadTreeNode* small);
  void node_subdivide(QuadTreeNode* node);
  void consolidate_node(QuadTreeNode* node);
  void reset();
  void reset(Boundary* boundary_);
  void copy_into(QuadTree* new_tree) const;
  void mark_neighbors_and_parents(QuadTreeNode* node);
  void perturb(const Boundary& new_boundary);
  void sort_leaves();

  void remove_inactive_dofs_at_level(int level);
  void remove_inactive_dofs_at_all_boxes();
  void remove_inactive_dofs_at_box(QuadTreeNode* node);
  void remove_hif_deactivated_dofs(QuadTreeNode* node);
  void populate_half_level_dofs(int level);
};



void get_descendents_updates(ki_Mat* updates,
                             const std::vector<int>& BN,
                             const QuadTreeNode* node,
                             std::set<const QuadTreeNode*>* visited_nodes = nullptr,
                             std::set<const HalfLevelNode*>* visited_halfnodes = nullptr);

void get_update(ki_Mat* updates, const std::vector<int>& BN,
                const QuadTreeNode* node,
                std::set<const QuadTreeNode*>* visited_nodes);

void get_update(ki_Mat* updates, const std::vector<int>& BN,
                const HalfLevelNode* node,
                std::set<const HalfLevelNode*>* visited_halfnodes);
void get_half_level_schur_updates(ki_Mat * updates,
                                  const std::vector<int>& BN,
                                  const HalfLevelNode * node,
                                  std::set<const QuadTreeNode*>* visited_nodes,
                                  std::set<const HalfLevelNode*>* visited_halfnodes) ;
}  // namespace kern_interp

#endif  // KERN_INTERP_QUADTREE_QUADTREE_H_
