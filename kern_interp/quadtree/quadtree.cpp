// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <iostream>
#include <queue>
#include <utility>
#include <boost/functional/hash.hpp>
#include "kern_interp/quadtree/quadtree.h"

namespace kern_interp {

// typedef std::pair<double, double> pair;

QuadTree::~QuadTree() {
  for (QuadTreeLevel* level : levels) {
    if (level) {
      delete level;
    }
  }
  levels.clear();
}


bool almost_align(const std::vector<double>& a,
                  const std::vector<double>& b) {
  for (int i = 0; i < a.size(); i++) {
    if (std::abs(a[i] - b[i]) > 1e-9) {
      return false;
    }
  }
  return true;
}


void QuadTree::compute_half_levels() {
  // Every node (except root) gets recompressor node for each side.
  // Before recompressor created, check if it's already been made.
  // Once recompressor created, check neighbors and add accordingly

  for (int level = 0; level < levels.size(); level++) {
    QuadTreeLevel* current_level = levels[level];
    for (int k = 0; k < current_level->nodes.size(); k++) {
      QuadTreeNode* node_a = current_level->nodes[k];

      std::vector<std::vector<double>> existent_recompressors;
      for (HalfLevelNode* recompressor : node_a->recompressor_nodes) {
        existent_recompressors.push_back(recompressor->center);
      }
      // iterate over sides
      for (int d = 0; d < domain_dimension; d++) {
        // each dimension has two sides
        for (int side_parity = -1; side_parity <= 1; side_parity += 2) {
          std::vector<double> newcenter = node_a->center;
          newcenter[d] += side_parity * node_a->side_length / 2.0;
          bool already_exists = false;
          for (std::vector<double> existent_recompressor : existent_recompressors) {
            if (almost_align(existent_recompressor, newcenter)) {
              already_exists = true;
              break;
            }
          }
          if (already_exists) continue;
          QuadTreeNode* sharing_container;
          HalfLevelNode* recompressor;
          bool node_created = false;
          // Check same-level neighbors to ensure one exists to share the node
          // with
          for (QuadTreeNode* neighbor : node_a->neighbors) {
            if (neighbor->level != node_a->level) continue;
            std::vector<double> othernewcenter = neighbor->center;
            othernewcenter[d] -= side_parity * neighbor->side_length / 2.0;
            if (almost_align(newcenter, othernewcenter)) {
              // Create the new halfnode
              recompressor = new HalfLevelNode();
              recompressor->center = newcenter;
              recompressor->side_length = node_a->side_length / sqrt(2);
              recompressor->partner_level = level;
              current_level->half_level->nodes.push_back(recompressor);
              node_a->recompressor_nodes.push_back(recompressor);
              neighbor->recompressor_nodes.push_back(recompressor);
              recompressor->containing_nodes.push_back(neighbor);
              recompressor->containing_nodes.push_back(node_a);
              node_created = true;
              sharing_container = neighbor;
              break;
            }
          }
        }
      }
    }
  }
}


void QuadTree::compute_neighbor_lists() {
  for (int level = 0; level < levels.size(); level++) {
    QuadTreeLevel* current_level = levels[level];
    for (int k = 0; k < current_level->nodes.size(); k++) {
      QuadTreeNode* node_a = current_level->nodes[k];
      // Each node is neighbors with all its siblings
      if (node_a->parent != nullptr) {
        for (QuadTreeNode* sibling : node_a->parent->children) {
          if (sibling != node_a) node_a->neighbors.push_back(sibling);
        }
        // Now check all parents' neighbors' children
        for (QuadTreeNode* parents_neighbor : node_a->parent->neighbors) {
          for (QuadTreeNode* cousin : parents_neighbor->children) {
            if (cousin == nullptr) continue;
            if (cousin->level != node_a->level) continue;
            double dist = 0.;
            for (int d = 0; d < domain_dimension; d++) {
              dist += pow(node_a->center[d] - cousin->center[d], 2);
            }
            dist = sqrt(dist);
            // just need to check if the distance of the bl corners
            // is <=s*sqrt(2)
            if (dist < node_a->side_length * sqrt(domain_dimension) + 1e-5) {
              node_a->neighbors.push_back(cousin);
            }
          }
        }
      }
      // now if it is a leaf, check against nodes in all subsequent levels
      if (node_a->is_leaf) {
        for (int n = 0; n < node_a->neighbors.size(); n++) {
          QuadTreeNode* neighbor =  node_a->neighbors[n];
          // make sure this isn't a neighbor from a higher level
          if (neighbor->level != node_a->level) {
            continue;
          }
          for (QuadTreeNode* child : neighbor->children) {
            if (child != nullptr) {
              get_descendent_neighbors(node_a, child);
            }
          }
        }
      }
    }
  }
}


void QuadTree::initialize_tree(Boundary* boundary,
                               int solution_dimension_,
                               int domain_dimension_) {
  assert(boundary->points.size() > 0
         && "number of boundary_points to init tree cannot be 0.");
  this->boundary_points = boundary->points;
  this->solution_dimension = solution_dimension_;
  this->domain_dimension = domain_dimension_;
  min = boundary_points[0];
  max = boundary_points[0];

  for (double point : boundary_points) {
    if (point < min) min = point;
    if (point > max) max = point;
  }

  double tree_min = min - 0.01;
  double tree_max = max + 0.01;

  root = new QuadTreeNode();

  root->level = 0;
  root->parent = nullptr;

  for (int i = 0; i < domain_dimension; i++) {
    root->center.push_back((tree_min + tree_max) / 2.0);
  }
  root->side_length = tree_max - tree_min;
  QuadTreeLevel* level_one = new QuadTreeLevel();
  level_one->half_level = new HalfLevel();
  level_one->nodes.push_back(root);
  levels.push_back(level_one);

  for (int i = 0; i < boundary_points.size(); i += domain_dimension) {
    std::vector<double> pt;
    for (int j = i; j < i + domain_dimension; j++) {
      pt.push_back(boundary_points[j]);
    }
    recursive_add(this->root, pt, i / domain_dimension);
  }
  compute_neighbor_lists();
  compute_half_levels();
  sort_leaves();
}


// adds neighbors to leaf which are on lower levels, by recursing and checking
// if the corners are along the wall.
void QuadTree::get_descendent_neighbors(QuadTreeNode* big,
                                        QuadTreeNode* small) {
  assert(big->level < small->level);
  bool are_neighbors = true;
  for (int d = 0; d < domain_dimension; d++) {
    double dist =  abs(big->center[d] - small->center[d]);
    double pred_dist = (big->side_length + small->side_length) / 2.;
    if (dist - pred_dist > 1e-14) {
      are_neighbors = false;
      break;
    }
  }
  if (are_neighbors) {
    big->neighbors.push_back(small);
    small->neighbors.push_back(big);
  }
  for (QuadTreeNode* child : small->children) {
    if (child != nullptr) {
      get_descendent_neighbors(big, child);
    }
  }
}


void QuadTree::recursive_add(QuadTreeNode* node, std::vector<double> pt,
                             int point_ind) {
  assert(node != nullptr && "recursive_add fails on null node.");
  for (int i = 0; i < solution_dimension; i++) {
    node->dof_lists.original_box.push_back(solution_dimension * point_ind + i);
  }
  // // figure out which child
  if (!node->is_leaf) {
    int child_idx = 0;
    for (int i = 0; i < domain_dimension; i++) {
      if (pt[i] >= node->center[i]) {
        child_idx += pow(2, domain_dimension - i - 1);
      }
    }
    recursive_add(node->children[child_idx], pt, point_ind);
  } else {
    // do we need one?
    // if this node is exploding and needs children

    if (node->is_leaf
        && node->dof_lists.original_box.size() > MAX_LEAF_DOFS) {
      node_subdivide(node);
    }
  }
}


// 1) gives node its four children
// 2) puts these children in their proper level
// 3) gives these children their corners
void QuadTree::node_subdivide(QuadTreeNode* node) {
  assert(node != nullptr && "node_subdivide fails on null node.");
  node->is_leaf = false;
  for (int child_idx = 0; child_idx < pow(2, domain_dimension); child_idx++) {
    std::vector<double> child_center;
    int tmp = child_idx;
    for (int d = 0; d < domain_dimension; d++) {
      if (tmp >= pow(2, domain_dimension - d - 1)) {
        child_center.push_back(node->center[d] + (node->side_length / 4.));
        tmp -= pow(2, domain_dimension - d - 1);
      } else {
        child_center.push_back(node->center[d] - (node->side_length / 4.));
      }
    }
    assert(tmp == 0);
    QuadTreeNode* child = new QuadTreeNode();
    child->center = child_center;
    node->children.push_back(child);
  }


  for (QuadTreeNode* child : node->children) {
    child->level = node->level + 1;
    child->side_length = node->side_length / 2.0;
    child->parent = node;
  }
  if (levels.size() < node->level + 2) {
    QuadTreeLevel* new_level = new QuadTreeLevel();
    new_level->half_level = new HalfLevel();
    levels.push_back(new_level);
    for (QuadTreeNode* child : node->children) {
      new_level->nodes.push_back(child);
    }
  } else {
    for (QuadTreeNode* child : node->children) {
      levels[node->level + 1]->nodes.push_back(child);
    }
  }
  // Now we bring the indices from the parent's box down into its childrens
  // boxes
  for (int index = 0; index < node->dof_lists.original_box.size();
       index += solution_dimension) {
    int matrix_index = node->dof_lists.original_box[index];
    int points_vec_index = (matrix_index / solution_dimension) *
                           domain_dimension;

    // double x = boundary_points[points_vec_index];
    // double y = boundary_points[points_vec_index + 1];

    int child_idx = 0;
    for (int i = 0; i < domain_dimension; i++) {
      if (boundary_points[points_vec_index + i] >= node->center[i]) {
        child_idx += pow(2, domain_dimension - i - 1);
      }
    }
    for (int i = 0; i < solution_dimension; i++) {
      node->children[child_idx]->dof_lists.original_box.push_back(matrix_index + i);
    }
  }
  for (QuadTreeNode* child : node->children) {
    if (child->dof_lists.original_box.size()  > MAX_LEAF_DOFS) {
      node_subdivide(child);
    }
  }
}


void QuadTree::mark_neighbors_and_parents(QuadTreeNode * node) {
  if (node == nullptr) return;
  if (node->p_marked) return;
  node->p_marked = true;
  node->compressed = false;
  node->X_rr_is_LU_factored = false;
  node->dof_lists.active_box.clear();
  node->dof_lists.skel.clear();
  node->dof_lists.skelnear.clear();
  node->dof_lists.redundant.clear();
  node->dof_lists.permutation.clear();
  node->T = ki_Mat(0, 0);
  node->U = ki_Mat(0, 0);
  node->L = ki_Mat(0, 0);

  for (QuadTreeNode* neighbor : node->neighbors) {
    // These aren't p_marked as they don't contain perturbations within them,
    // but need updating due to their neighbors.
    neighbor->compressed = false;
    neighbor->X_rr_is_LU_factored = false;
    neighbor->dof_lists.active_box.clear();
    neighbor->dof_lists.skel.clear();
    neighbor->dof_lists.skelnear.clear();
    neighbor->dof_lists.redundant.clear();
    neighbor->dof_lists.permutation.clear();
    neighbor->T = ki_Mat(0, 0);
    neighbor->U = ki_Mat(0, 0);
    neighbor->L = ki_Mat(0, 0);
    mark_neighbors_and_parents(neighbor->parent);
  }
  mark_neighbors_and_parents(node->parent);
}


void QuadTree::consolidate_node(QuadTreeNode* node) {
  // Need to
  //  Move leaf child dofs into my original box
  //  erase all descendents from levels
  //  delete immediate descentdents
  node->dof_lists.original_box.clear();
  std::vector<QuadTreeNode*> remove_from_lvl;
  std::vector<QuadTreeNode*> queue;
  queue.push_back(node);
  for (int i = 0; i < queue.size(); i++) {
    QuadTreeNode* current = queue[i];
    if (current->is_leaf) {
      node->dof_lists.original_box.insert(
        node->dof_lists.original_box.end(),
        current->dof_lists.original_box.begin(),
        current->dof_lists.original_box.end());
    } else {
      for (QuadTreeNode* child : current->children) {
        queue.push_back(child);
      }
    }
    if (current != node) {
      remove_from_lvl.push_back(current);
    }
  }
  for (QuadTreeNode* erase : remove_from_lvl) {
    QuadTreeLevel* erase_level = levels[erase->level];
    for (int i = 0; i < erase_level->nodes.size(); i++) {
      if (erase_level->nodes[i] == erase) {
        QuadTreeNode* del = erase_level->nodes[i];
        erase_level->nodes.erase(erase_level->nodes.begin() + i);
        delete del;
        break;
      }
    }
  }
  node->children.clear();
  node->is_leaf = true;
}


void QuadTree::sort_leaves() {
  for (int level = 0; level < levels.size(); level++) {
    QuadTreeLevel* current_level = levels[level];
    for (int k = 0; k < current_level->nodes.size(); k++) {
      QuadTreeNode* node = current_level->nodes[k];
      if ((!node->is_leaf) || node->compressed) continue;
      std::vector<std::pair<double, int>> origboxtmp;
      for (int i = 0; i < node->dof_lists.original_box.size(); i++) {
        int points_vec_index = domain_dimension * (node->dof_lists.original_box[i] /
                               solution_dimension);
        double x = boundary_points[points_vec_index];
        double y = boundary_points[points_vec_index + 1];
        double dist = 0.;
        for (int d = 0; d < domain_dimension; d++) {
          dist += pow(node->center[d] - boundary_points[points_vec_index + d], 2);
        }
        dist = sqrt(dist);
        origboxtmp.push_back(std::pair<double, int>(dist,
                             node->dof_lists.original_box[i]));
      }
      std::sort(origboxtmp.begin(),
                origboxtmp.end());
      node->dof_lists.original_box.clear();
      for (int i = 0; i < origboxtmp.size(); i++) {
        node->dof_lists.original_box.push_back(origboxtmp[i].second);
      }
    }
  }
}



void QuadTree::perturb(const Boundary & perturbed_boundary) {
  // 1) create mapping, storing vectors of additions/deletions
  // 2) go to every node, marking those with additions and deletions
  // these are vectors of point indices (p_0, p_1, etc)
  std::vector<double> additions;
  std::vector<double> deletions;

  std::vector<double> old_points = boundary_points;
  std::vector<double> new_points = perturbed_boundary.points;
  // now create mapping of new_points to their point index in the new vec
  std::unordered_map<std::vector<double>, int, boost::hash<std::vector<double>>>
  point_to_new_index;
  for (int i = 0; i < new_points.size(); i += domain_dimension) {
    std::vector<double> new_point;
    for (int j = 0; j < domain_dimension; j++) {
      new_point.push_back(new_points[i + j]);
    }
    point_to_new_index[new_point] = i / domain_dimension;
  }

  std::vector<bool> found_in_old(new_points.size() / domain_dimension);
  for (int i = 0; i < found_in_old.size(); i++) {
    found_in_old[i] = false;
  }
  // Mapping from point index in old points vec to point index in new points vec
  std::unordered_map<int, int> old_index_to_new_index;
  for (int i = 0; i < old_points.size(); i += domain_dimension) {
    std::vector<double> old_point;
    for (int j = 0; j < domain_dimension; j++) {
      old_point.push_back(old_points[i + j]);
    }
    // Is this point also in the new points vec?
    std::unordered_map<std::vector<double>, int, boost::hash<std::vector<double>>>::const_iterator
    element =
      point_to_new_index.find(old_point);
    if (element != point_to_new_index.end()) {
      old_index_to_new_index[i / domain_dimension] = element->second;
      found_in_old[element->second] = true;
    } else {
      // If it's in the old vec but not the new vec, it was deleted
      deletions.push_back(i / domain_dimension);
    }
  }
  for (int i = 0; i < found_in_old.size(); i++) {
    if (!found_in_old[i]) {
      additions.push_back(i);
    }
  }

  // go through all leaf original box vectors and apply mapping.
  // (if there is a deletion it will be processed later)
  // each node will be one of three things
  //   1) unmarked, in which case the below is a perfectly good mapping
  //   2) marked non-leaf, the below is irrelevant, everything will be dumped
  //   3) marked leaf, only the leaf portion of the below is relevant.
  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      std::vector<int> ob, ab, s, r, sn, n;
      if (node->is_leaf) {
        for (int idx : node->dof_lists.original_box) {
          int point_index = idx / solution_dimension;
          std::unordered_map<int, int>::const_iterator element =
            old_index_to_new_index.find(point_index);
          if (element != old_index_to_new_index.end()) {
            ob.push_back(solution_dimension * element->second
                         + idx % solution_dimension);
          }
        }
        node->dof_lists.original_box = ob;
      }
      for (int idx : node->dof_lists.active_box) {
        int point_index = idx / solution_dimension;
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(point_index);
        if (element != old_index_to_new_index.end()) {
          ab.push_back(solution_dimension * element->second
                       + idx % solution_dimension);
        }
      }
      for (int idx : node->dof_lists.skel) {
        int point_index = idx / solution_dimension;
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(point_index);
        if (element != old_index_to_new_index.end()) {
          s.push_back(solution_dimension * element->second
                      + idx % solution_dimension);
        }
      }
      for (int idx : node->dof_lists.redundant) {
        int point_index = idx / solution_dimension;
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(point_index);
        if (element != old_index_to_new_index.end()) {
          r.push_back(solution_dimension * element->second
                      + idx % solution_dimension);
        }
      }
      for (int idx : node->dof_lists.skelnear) {
        int point_index = idx / solution_dimension;
        std::unordered_map<int, int>::const_iterator element =
          old_index_to_new_index.find(point_index);
        if (element != old_index_to_new_index.end()) {
          sn.push_back(solution_dimension * element->second
                       + idx % solution_dimension);
        }
      }
      node->dof_lists.active_box = ab;
      node->dof_lists.skel = s;
      node->dof_lists.skelnear = sn;
      node->dof_lists.redundant = r;
    }
  }

  // go through all additions, find their leaves, make addition and call mark
  // function
  std::vector<QuadTreeNode*> maybe_bursting;
  for (int i = 0; i < additions.size(); i++) {
    // double newx = new_points[2 * additions[i]];
    // double newy = new_points[2 * additions[i] + 1];
    std::vector<double> newpt;
    for (int j = 0; j < domain_dimension; j++) {
      newpt.push_back(new_points[domain_dimension * additions[i] + j]);
    }
    QuadTreeNode* current = root;
    while (!current->is_leaf) {
      int child_idx = 0;
      for (int j = 0; j < domain_dimension; j++) {
        if (newpt[j] >= current->center[j]) {
          child_idx += pow(2, domain_dimension - j - 1);
        }
      }
      current = current->children[child_idx];
    }
    for (int j = 0; j < solution_dimension; j++) {
      current->dof_lists.original_box.push_back(solution_dimension
          * additions[i] + j);
    }
    maybe_bursting.push_back(current);
    mark_neighbors_and_parents(current);
  }

  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      if (node->is_leaf) {
        node->dofs_below = node->dof_lists.original_box.size();
      } else {
        node->dofs_below = 0;
      }
    }
  }
  for (int l = levels.size() - 1; l >= 1; l--) {
    QuadTreeLevel* level = levels[l];
    for (QuadTreeNode* node : level->nodes) {
      node->parent->dofs_below += node->dofs_below;
    }
  }

  // go through all deletions, find their leaves, make deletion and call mark
  // function
  std::unordered_map<QuadTreeNode*, bool> sparse;
  for (int i = 0; i < deletions.size(); i++) {
    QuadTreeNode* current = root;
    // oldpt is the current point that's being deleted
    std::vector<double> oldpt;
    for (int j = 0; j < domain_dimension; j++) {
      oldpt.push_back(old_points[domain_dimension * deletions[i] + j]);
    }

    bool path_marked = false;
    // Starting from the root, mark all boxes containing deleted point, identify
    // first that's sparse and can be consolidated.
    while (!current->is_leaf) {
      if (current->dofs_below < MAX_LEAF_DOFS && !path_marked) {
        path_marked = true;
        sparse[current] = true;
      }

      int child_idx = 0;
      for (int j = 0; j < domain_dimension; j++) {
        if (oldpt[j] >= current->center[j]) {
          child_idx += pow(2, domain_dimension - j - 1);
        }
      }
      current = current->children[child_idx];
    }
    // Now current is the leaf containing the deleted point
    mark_neighbors_and_parents(current);
  }

  this->boundary_points = perturbed_boundary.points;
  // If any nodes are bursting now, subdivide them.
  for (QuadTreeNode* node : maybe_bursting) {
    if (node->is_leaf
        && node->dof_lists.original_box.size() > MAX_LEAF_DOFS) {
      node_subdivide(node);
      // Is this necessary?
      for (QuadTreeNode* child : node->children) {
        mark_neighbors_and_parents(child);
      }
    }
  }
// If we can consolidate nodes into their parent, do that.
  for (auto it = sparse.begin(); it != sparse.end(); ++it) {
    consolidate_node(it->first);
  }
  for (QuadTreeLevel* level : levels) {
    for (QuadTreeNode* node : level->nodes) {
      node->p_marked = false;
      node->neighbors.clear();
    }
  }
  compute_neighbor_lists();

  sort_leaves();
}


void copy_info(QuadTreeNode* old_node, QuadTreeNode* new_node) {
  new_node->level = old_node->level;
  new_node->dofs_below = old_node->dofs_below;
  new_node->is_leaf = old_node->is_leaf;
  new_node->X_rr_is_LU_factored = old_node->X_rr_is_LU_factored;
  new_node->compressed = old_node->compressed;
  new_node->side_length = old_node->side_length;
  new_node->dof_lists = old_node->dof_lists;
  new_node->T = old_node->T;
  new_node->L = old_node->L;
  new_node->U = old_node->U;
  new_node->X_rr = old_node->X_rr;
  new_node->schur_update = old_node->schur_update;
  new_node->X_rr_lu = old_node->X_rr_lu;
  new_node->X_rr_piv = old_node->X_rr_piv;
  new_node->center = old_node->center;
  // for (int i = 0; i < 8; i++) {
  //   new_node->corners[i] = old_node->corners[i];
  // }
}

// TODO(HIF)
void QuadTree::copy_into(QuadTree* new_tree) const {
  // The strategy here is going to be to create a new node for every old node,
  // then keep a mapping from new to old. With that, we'll copy all the data
  // over, including connections, levels, and matrices.
  *new_tree = QuadTree();
  std::vector < QuadTreeNode*> new_nodes;
  std::unordered_map<QuadTreeNode*, QuadTreeNode*> old_to_new;
  std::unordered_map<QuadTreeNode*, QuadTreeNode*> new_to_old;
  for (int lvl = 0; lvl < levels.size(); lvl++) {
    QuadTreeLevel* level = levels[lvl];
    for (int n = 0; n < level->nodes.size(); n++) {
      QuadTreeNode* old_node = level->nodes[n];
      QuadTreeNode* new_node = new QuadTreeNode();
      new_nodes.push_back(new_node);
      copy_info(old_node, new_node);
      old_to_new[old_node] = new_node;
      new_to_old[new_node] = old_node;
    }
  }

  for (int n = 0; n < new_nodes.size(); n++) {
    QuadTreeNode* new_node = new_nodes[n];
    new_node->parent = old_to_new[new_to_old[new_node]->parent];
    if (!new_to_old[new_node]->is_leaf) {
      for (int c = 0; c < pow(2, domain_dimension); c++) {
        new_node->children.push_back(old_to_new[new_to_old[new_node]->children[c]]);
      }
    }
    for (int nbr = 0; nbr < new_to_old[new_node]->neighbors.size(); nbr++) {
      QuadTreeNode* neighbor = new_to_old[new_node]->neighbors[nbr];
      new_node->neighbors.push_back(old_to_new[neighbor]);
    }
  }

  new_tree->root = old_to_new[root];
  new_tree->solution_dimension = solution_dimension;
  new_tree->domain_dimension = domain_dimension;
  new_tree->no_proxy_level = no_proxy_level;
  new_tree->min = min;
  new_tree->max = max;
  new_tree->allskel_mat = allskel_mat;
  new_tree->allskel_mat_lu = allskel_mat_lu;
  new_tree->U = U;
  new_tree->Psi = Psi;
  new_tree->S_LU = S_LU;
  new_tree->allskel_mat_piv = allskel_mat_piv;
  new_tree->S_piv = S_piv;

  for (int lvl = 0; lvl < levels.size(); lvl++) {
    QuadTreeLevel* old_level = levels[lvl];
    QuadTreeLevel* new_level = new QuadTreeLevel();
    for (int n = 0; n < old_level->nodes.size(); n++) {
      new_level->nodes.push_back(old_to_new[old_level->nodes[n]]);
    }
    new_tree->levels.push_back(new_level);
  }
}


void QuadTree::reset(Boundary * boundary_) {
  for (QuadTreeLevel* level : levels) {
    if (level) {
      delete level;
    }
  }
  levels.clear();
  initialize_tree(boundary_, solution_dimension,
                  domain_dimension);
}


void QuadTree::remove_inactive_dofs_at_all_boxes() {
  int lvls = levels.size();
  for (int level = lvls - 1; level >= 0; level--) {
    remove_inactive_dofs_at_level(level);
  }
}


void QuadTree::remove_inactive_dofs_at_level(int level) {
  QuadTreeLevel* current_level = levels[level];
  // First, get all active dofs from children
  for (QuadTreeNode* node : current_level->nodes) {
    if (node->compressed) continue;
    remove_inactive_dofs_at_box(node);
  }
  // Next, get all active near dofs from neighbors
  for (QuadTreeNode* node_a : current_level->nodes) {
    node_a->dof_lists.near.clear();
    for (QuadTreeNode* neighbor : node_a->neighbors) {
      // Some neighbors are smaller boxes from higher levels, we don't
      // care about those, their parents have the updated information.
      if (neighbor->level > node_a->level) {
        continue;
      }
      if (neighbor->is_leaf) {
        for (int idx : neighbor->dof_lists.original_box) {
          node_a->dof_lists.near.push_back(idx);
        }
      } else {
        for (int idx : neighbor->dof_lists.active_box) {
          node_a->dof_lists.near.push_back(idx);
        }
      }
    }
  }
}


void QuadTree::remove_inactive_dofs_at_box(QuadTreeNode* node) {
  // this function removes from the box any DoFs which have already been made
  // redundant. It involves a bunch of annoying C++ functions and probably
  // would look nicer in matlab.

  // populate active_box
  node->dof_lists.skel.clear();
  node->dof_lists.skelnear.clear();
  node->dof_lists.redundant.clear();
  node->dof_lists.active_box.clear();

  if (!node->is_leaf) {
    for (QuadTreeNode* child : node->children) {
      if (child->compressed) {
        for (int i : child->dof_lists.skel) {
          node->dof_lists.active_box.push_back(i);
        }
      } else {
        for (int i : child->dof_lists.active_box) {
          node->dof_lists.active_box.push_back(i);
        }
      }
    }
    remove_hif_deactivated_dofs(node);
  } else {
    node->dof_lists.active_box = node->dof_lists.original_box;
  }
}

void QuadTree::remove_hif_deactivated_dofs(QuadTreeNode* node) {
  std::vector<int> hif_deactivated_dofs;
  std::set<HalfLevelNode*> visited_halfnodes;
  for (QuadTreeNode* child : node->children) {
    for (HalfLevelNode* recompressor_node : child->recompressor_nodes) {
      if (visited_halfnodes.find(recompressor_node) != visited_halfnodes.end()) {
        continue;
      }
      visited_halfnodes.insert(recompressor_node);
      hif_deactivated_dofs.insert(hif_deactivated_dofs.end(),
                                  recompressor_node->dof_lists.redundant.begin(),
                                  recompressor_node->dof_lists.redundant.end());
    }
  }

  std::sort(hif_deactivated_dofs.begin(), hif_deactivated_dofs.end());
  std::sort(node->dof_lists.active_box.begin(),
            node->dof_lists.active_box.end());
  std::vector<int> difference;
  std::set_difference(
    node->dof_lists.active_box.begin(),  node->dof_lists.active_box.end(),
    hif_deactivated_dofs.begin(), hif_deactivated_dofs.end(),
    std::back_inserter(difference)
  );
  node->dof_lists.active_box = difference;
}


void QuadTree::populate_half_level_dofs(int level) {
  double diam = levels[level]->nodes[0]->side_length;
  int correct = 0;
  // Go through boxes on this level, assign dofs to relevant half nodes

  for (QuadTreeNode* node : levels[level]->nodes) {
    if (node->recompressor_nodes.empty()) continue;
    std::vector<int> dofs;
    if (node->compressed) {
      dofs = node->dof_lists.skel;
    } else if (node->is_leaf) {
      dofs = node->dof_lists.original_box;
    } else {
      dofs = node->dof_lists.active_box;
    }
    for (int i = 0; i < dofs.size(); i++) {
      int idx = dofs[i];
      int points_vec_index = (idx / solution_dimension) *
                             domain_dimension;
      double mindist = node->side_length * 10;
      std::vector<double> closest_face_center;
      // Check face centers
      for (int d = 0; d < domain_dimension; d++) {
        // each dimension has two sides
        for (int side_parity = -1; side_parity <= 1; side_parity += 2) {
          std::vector<double> facecenter = node->center;
          facecenter[d] += side_parity * node->side_length / 2.0;
          double facedist = 0;
          for (int dd = 0; dd < domain_dimension; dd++) {
            facedist += pow(facecenter[dd] - boundary_points[points_vec_index + dd], 2);
          }
          if (sqrt(facedist) < mindist) {
            closest_face_center = facecenter;
            mindist = sqrt(facedist);
          }
        }
      }
      // Now check if closest facecenter corresponds to a recomp node
      for (HalfLevelNode* recompressor_node : node->recompressor_nodes) {
        if (almost_align(closest_face_center, recompressor_node->center)) {
          recompressor_node->dof_lists.active_box.push_back(idx);
          break;
        }
      }
    }
  }

  // Near should be any dof in containers or their same or larger sized neighbors
  for (HalfLevelNode* halfnode : levels[level]->half_level->nodes) {
    std::set<QuadTreeNode*> visited_nodes;
    // TODO(HIF) set complexity bad? Use hash?
    std::set<int> halfnode_contained_dofs;
    for (int idx : halfnode->dof_lists.active_box) {
      halfnode_contained_dofs.insert(idx);
    }
    for (QuadTreeNode* containing_node : halfnode->containing_nodes) {
      for (QuadTreeNode* neighbor : containing_node->neighbors) {
        // skip neighbors from higher levels, parents have info
        if (neighbor->level > halfnode->partner_level) continue;
        if (visited_nodes.find(neighbor) != visited_nodes.end()) {
          continue;
        }
        visited_nodes.insert(neighbor);

        // If neighbor is same level, either active or skel has relevant dofs
        // If it is from lower level, it is leaf with original dofs
        std::vector<int> dof_list;
        if (neighbor->compressed) {
          dof_list = neighbor->dof_lists.skel;
        } else  if (neighbor->is_leaf) {
          dof_list = neighbor->dof_lists.original_box;
        }  else {
          dof_list = neighbor->dof_lists.active_box;
        }
        for (int idx : dof_list) {
          if (halfnode_contained_dofs.find(idx) != halfnode_contained_dofs.end()) {
            continue;
          }
          halfnode->dof_lists.near.push_back(idx);
          halfnode_contained_dofs.insert(idx);
        }
      }
    }
  }
}


void get_half_level_schur_updates(ki_Mat * updates,
                                  const std::vector<int>& BN, const HalfLevelNode * node,
                                  std::set<const QuadTreeNode*>* visited_nodes,
                                  std::set<const HalfLevelNode*>* visited_halfnodes) {
  assert(node != nullptr &&
         "get_half_level_schur_updates fails on null node.");
  assert(BN.size() > 0 &&
         "get_half_level_schur_updates needs positive num of DOFs");

  std::set<const QuadTreeNode*> visited_nodes_ = std::set<const QuadTreeNode*>();
  std::set<const HalfLevelNode*> visited_halfnodes_ =
    std::set<const HalfLevelNode*>();
  if (visited_nodes == nullptr) {
    visited_nodes = &visited_nodes_;
    visited_halfnodes = &visited_halfnodes_;
  }

  for (QuadTreeNode* containing_node : node->containing_nodes) {
    if (containing_node->compressed) get_update(updates, BN, containing_node,
          visited_nodes);
    get_descendents_updates(updates, BN,
                            containing_node, visited_nodes, visited_halfnodes);

  }
}


void get_descendents_updates(ki_Mat * updates,
                             const std::vector<int>& BN, const QuadTreeNode * node,
                             std::set<const QuadTreeNode*>* visited_nodes,
                             std::set<const HalfLevelNode*>* visited_halfnodes) {
  assert(node != nullptr && "get_descendents_updates fails on null node.");

  std::set<const QuadTreeNode*> visited_nodes_ = std::set<const QuadTreeNode*>();
  std::set<const HalfLevelNode*> visited_halfnodes_ =
    std::set<const HalfLevelNode*>();
  if (visited_nodes == nullptr) {
    visited_nodes = &visited_nodes_;
    visited_halfnodes = &visited_halfnodes_;
  }

  for (QuadTreeNode* child : node->children) {
    if (child->compressed) get_update(updates, BN, child, visited_nodes);
    for (HalfLevelNode* recompressor_node : child->recompressor_nodes) {
      if (recompressor_node->compressed) {
        get_update(updates, BN, recompressor_node, visited_halfnodes);
      }
    }
    get_descendents_updates(updates, BN, child,
                            visited_nodes, visited_halfnodes);
  }
}


void get_update(ki_Mat * update,
                const std::vector<int>& BN,
                const QuadTreeNode * node,
                std::set<const QuadTreeNode*>* visited_nodes)  {
  // Node needs to check all its dofs against BN, enter interactions into
  // corresponding locations
  // Node only updated its own BN dofs, and the redundant ones are no longer
  // relevant, so we only care about child's SN dofs
  // First create a list of Dofs that are also in node's skelnear,
  // and with each one give the index in skelnear and the index in BN

  if (visited_nodes->find(node) != visited_nodes->end()) {
    return;
  }
  visited_nodes->insert(node);
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

void get_update(ki_Mat * update,
                const std::vector<int>& BN,
                const HalfLevelNode * node,
                std::set<const HalfLevelNode*>* visited_halfnodes)  {
  // Node needs to check all its dofs against BN, enter interactions into
  // corresponding locations
  // Node only updated its own BN dofs, and the redundant ones are no longer
  // relevant, so we only care about child's SN dofs
  // First create a list of Dofs that are also in node's skelnear,
  // and with each one give the index in skelnear and the index in BN
  if (visited_halfnodes->find(node) != visited_halfnodes->end()) {
    return;
  }
  visited_halfnodes->insert(node);
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


}  // namespace kern_interp
