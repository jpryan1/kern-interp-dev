// Copyright 2019 John Paul Ryan
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <algorithm>
#include "kern_interp/kernel/kernel.h"

namespace kern_interp {


inline double Kernel::boundary_normals(int i) const {
  return (boundary_normals_.empty() ? 0. : boundary_normals_[i]);
}

inline double Kernel::boundary_weights(int i) const {
  return (boundary_weights_.empty() ? 0. : boundary_weights_[i]);
}

inline double Kernel::boundary_curvatures(int i) const {
  return (boundary_curvatures_.empty() ? 0. : boundary_curvatures_[i]);
}


void Kernel::update_data(Boundary* boundary) {
  boundary_points_ = boundary->points;
  boundary_normals_ = boundary->normals;
  boundary_weights_ = boundary->weights;
  boundary_curvatures_ = boundary->curvatures;
}


void Kernel::one_d_kern(int mat_idx, ki_Mat* ret, double r1, double r2,
                        double tn1, double tn2, double sn1, double sn2,
                        double sw, double sc, bool forward) const {
  if (pde == Pde::LAPLACE) {
    if (r1 == 0. && r2 == 0.) {
      ret->mat[mat_idx] = 0.5 + 0.5 * sc * sw * (1.0 / (2 * M_PI));
    } else {
      ret->mat[mat_idx] =  -sw * (1.0 / (2 * M_PI)) *
                           (r1 * sn1 + r2 * sn2) /
                           (r1 * r1 + r2 * r2);
    }
  } else if (pde == Pde::LAPLACE_NEUMANN) {
    if (forward) {
      ret->mat[mat_idx] =  sw * (1.0 / (2 * M_PI)) * log(sqrt(
                             r1 * r1 + r2 * r2));
    } else if (r1 == 0. && r2 == 0.) {
      ret->mat[mat_idx] = sw - 0.5 + 0.5 * sc * sw * (1.0 / (2 * M_PI));
    } else {
      ret->mat[mat_idx] =  sw + sw * (1.0 / (2 * M_PI)) *
                           (r1 * tn1 + r2 * tn2) /
                           (r1 * r1 + r2 * r2);
    }
  } else if (pde == Pde::GAUSS) {
    if (r1 == 0. && r2 == 0.) {
      ret->mat[mat_idx] = 2.;
    } else {
      ret->mat[mat_idx] = exp(-(pow(r1, 2) + pow(r2, 2)));
    }
  }
}


void Kernel::two_d_kern(int mat_idx, int tgt_parity, int src_parity,
                        ki_Mat* ret, double r1, double r2,
                        double tn1, double tn2, double sn1, double sn2,
                        double sw, double sc, bool forward) const {
  double sn[2] = {sn1, sn2};
  double tn[2] = {tn1, tn2};
  double r[2] = {r1, r2};
  double scale = 1 / M_PI;
  int is_diag = (tgt_parity + src_parity + 1) % 2;
  if (pde == Pde::STOKES) {
    if (r1 == 0. && r2 == 0.) {
      ret->mat[mat_idx] = - 0.5 * is_diag  + (1 - 2 * is_diag) * 0.5 * sc
                          * sw * scale * sn[(tgt_parity + 1) % 2] *
                          sn[(src_parity + 1) % 2 ] + sw * tn[tgt_parity]
                          * sn[src_parity];
    } else {
      ret->mat[mat_idx] = sw * scale * (r1 * sn[0] + r2 * sn[1]) /
                          (pow(r1 * r1 + r2 * r2, 2))
                          * r[tgt_parity] * r[src_parity]
                          + sw * tn[tgt_parity] * sn[src_parity ];
    }
  }
}

void Kernel::three_d_stokes(int mat_idx, int tgt_parity, int src_parity,
                            ki_Mat* ret, double r1, double r2, double r3, const ki_Mat& diag,
                            double tn1, double tn2, double tn3, double sn1, double sn2, double sn3,
                            double sw) const {
  double sn[3] = {sn1, sn2, sn3};
  double tn[3] = {tn1, tn2, tn3};
  double r[3] = {r1, r2, r3};
  double scale = -3.0 / (4.0 * M_PI);

  if (r1 == 0. && r2 == 0. && r3 == 0.) {
    ret->mat[mat_idx] = diag.get(tgt_parity, src_parity);//TODO
  } else {
    ret->mat[mat_idx] = sw * scale * (r1 * sn1 + r2 * sn2 + r3 * sn3) /
                        (pow(r1 * r1 + r2 * r2 + r3 * r3, 2.5))
                        * r[tgt_parity] * r[src_parity]
                        + sw * tn[tgt_parity] * sn[src_parity ];
  }

}

void Kernel::three_d_laplace(int mat_idx, ki_Mat* ret, double r1, double r2,
                             double r3,
                             double diag,
                             double sn1, double sn2, double sn3,
                             double sw) const {
  if (r1 == 0. && r2 == 0. && r3 == 0.) {
    ret->mat[mat_idx] = diag;
  } else {
    ret->mat[mat_idx] =  -sw * (1.0 / (4 * M_PI)) *
                         (r1 * sn1 + r2 * sn2 + r3 * sn3) /
                         pow(r1 * r1 + r2 * r2 + r3 * r3, 1.5);
  }
}


void Kernel::compute_diag_entries_3dlaplace(Boundary* boundary) {
  boundary_diags = std::vector<double>(boundary->weights.size());
  #pragma omp parallel for num_threads(8)
  for (int pt_idx = 0; pt_idx < boundary->num_outer_nodes; pt_idx++) {
    double tp1 = boundary->points[3 * pt_idx];
    double tp2 = boundary->points[3 * pt_idx + 1];
    double tp3 = boundary->points[3 * pt_idx + 2];

    double rowsum = 0.0;
    for (int other_pt_idx = 0; other_pt_idx < boundary->num_outer_nodes;
         other_pt_idx++) {
      if (pt_idx == other_pt_idx) continue;

      double sp1 = boundary->points[3 * other_pt_idx];
      double sp2 = boundary->points[3 * other_pt_idx + 1];
      double sp3 = boundary->points[3 * other_pt_idx + 2];
      double sn1 = boundary->normals[3 * other_pt_idx];
      double sn2 = boundary->normals[3 * other_pt_idx + 1];
      double sn3 = boundary->normals[3 * other_pt_idx + 2];
      double sw = boundary->weights[other_pt_idx];

      double r1 = tp1 - sp1;
      double r2 = tp2 - sp2;
      double r3 = tp3 - sp3;

      rowsum += (-sw * (1.0 / (4 * M_PI)) *
                 (r1 * sn1 + r2 * sn2 + r3 * sn3) /
                 pow(r1 * r1 + r2 * r2 + r3 * r3, 1.5));
    }
    boundary_diags[pt_idx] = 1 - rowsum;
  }

  // #pragma omp parallel for num_threads(8)

  int curr_idx = boundary->num_outer_nodes;
  for (Hole hole : boundary->holes) {
    #pragma omp parallel for num_threads(8)
    for (int pt_idx = curr_idx; pt_idx < curr_idx + hole.num_nodes; pt_idx++) {
      double tp1 = boundary->points[3 * pt_idx];
      double tp2 = boundary->points[3 * pt_idx + 1];
      double tp3 = boundary->points[3 * pt_idx + 2];

      double rowsum = 0.0;
      for (int other_pt_idx = curr_idx; other_pt_idx < curr_idx + hole.num_nodes;
           other_pt_idx++) {
        if (pt_idx == other_pt_idx) continue;

        double sp1 = boundary->points[3 * other_pt_idx];
        double sp2 = boundary->points[3 * other_pt_idx + 1];
        double sp3 = boundary->points[3 * other_pt_idx + 2];
        double sn1 = boundary->normals[3 * other_pt_idx];
        double sn2 = boundary->normals[3 * other_pt_idx + 1];
        double sn3 = boundary->normals[3 * other_pt_idx + 2];
        double sw = boundary->weights[other_pt_idx];

        double r1 = tp1 - sp1;
        double r2 = tp2 - sp2;
        double r3 = tp3 - sp3;

        rowsum += (-sw * (1.0 / (4 * M_PI)) *
                   (r1 * sn1 + r2 * sn2 + r3 * sn3) /
                   pow(r1 * r1 + r2 * r2 + r3 * r3, 1.5));
      }
      boundary_diags[pt_idx] = -rowsum;
    }
    curr_idx += hole.num_nodes;
  }
}


void Kernel::compute_diag_entries_3dstokes(Boundary* boundary) {



  boundary_diag_tensors = std::vector<ki_Mat>(boundary->weights.size());

  for (int i = 0; i < boundary->weights.size(); i++) {
    boundary_diag_tensors[i] = ki_Mat(3, 3);
    //   boundary_diag_tensors[i].set(0,0,1);
    //   boundary_diag_tensors[i].set(1,1,1);
    //   boundary_diag_tensors[i].set(2,2,1);
  }
  // return;
  #pragma omp parallel for num_threads(8)
  for (int dof = 0; dof < boundary->num_outer_nodes * domain_dimension; dof++) {
    int pt_idx = dof / 3;
    double tp1 = boundary->points[3 * pt_idx];
    double tp2 = boundary->points[3 * pt_idx + 1];
    double tp3 = boundary->points[3 * pt_idx + 2];
    double tn1 = boundary->normals[3 * pt_idx];
    double tn2 = boundary->normals[3 * pt_idx + 1];
    double tn3 = boundary->normals[3 * pt_idx + 2];
    double tn[3] = {tn1, tn2, tn3};

    for (int other_dof = 0;
         other_dof < boundary->num_outer_nodes * domain_dimension; other_dof++) {
      int other_pt_idx = other_dof / 3;
      if (pt_idx == other_pt_idx) continue;

      double sp1 = boundary->points[3 * other_pt_idx];
      double sp2 = boundary->points[3 * other_pt_idx + 1];
      double sp3 = boundary->points[3 * other_pt_idx + 2];
      double sn1 = boundary->normals[3 * other_pt_idx];
      double sn2 = boundary->normals[3 * other_pt_idx + 1];
      double sn3 = boundary->normals[3 * other_pt_idx + 2];
      double sw = boundary->weights[other_pt_idx];
      double sn[3] = {sn1, sn2, sn3};
      double scale = -3.0 / (4.0 * M_PI);
      double r1 = tp1 - sp1;
      double r2 = tp2 - sp2;
      double r3 = tp3 - sp3;
      double r[3] = {r1, r2, r3};
      double pot = sw * scale * (r1 * sn1 + r2 * sn2 + r3 * sn3) /
                   (pow(r1 * r1 + r2 * r2 + r3 * r3, 2.5))
                   * r[dof % 3] * r[other_dof % 3]
                   + sw * tn[dof % 3] * sn[other_dof % 3 ];


      double tmp = boundary_diag_tensors[pt_idx].get(dof % 3, other_dof % 3);
      boundary_diag_tensors[pt_idx].set(dof % 3, other_dof % 3, tmp + pot);
    }
  }
  int curr_idx = boundary->num_outer_nodes * domain_dimension;

  for (Hole hole : boundary->holes) {
    #pragma omp parallel for num_threads(8)
    for (int dof = curr_idx; dof < curr_idx + hole.num_nodes * domain_dimension;
         dof++) {
      int pt_idx = dof / 3;
      double tp1 = boundary->points[3 * pt_idx];
      double tp2 = boundary->points[3 * pt_idx + 1];
      double tp3 = boundary->points[3 * pt_idx + 2];
      double tn1 = boundary->normals[3 * pt_idx];
      double tn2 = boundary->normals[3 * pt_idx + 1];
      double tn3 = boundary->normals[3 * pt_idx + 2];
      double tn[3] = {tn1, tn2, tn3};

      for (int other_dof = curr_idx;
           other_dof < curr_idx + hole.num_nodes * domain_dimension; other_dof++) {
        int other_pt_idx = other_dof / 3;
        if (pt_idx == other_pt_idx) continue;

        double sp1 = boundary->points[3 * other_pt_idx];
        double sp2 = boundary->points[3 * other_pt_idx + 1];
        double sp3 = boundary->points[3 * other_pt_idx + 2];
        double sn1 = boundary->normals[3 * other_pt_idx];
        double sn2 = boundary->normals[3 * other_pt_idx + 1];
        double sn3 = boundary->normals[3 * other_pt_idx + 2];
        double sw = boundary->weights[other_pt_idx];
        double sn[3] = {sn1, sn2, sn3};
        double scale = -3.0 / (4.0 * M_PI);
        double r1 = tp1 - sp1;
        double r2 = tp2 - sp2;
        double r3 = tp3 - sp3;
        double r[3] = {r1, r2, r3};
        double pot = sw * scale * (r1 * sn1 + r2 * sn2 + r3 * sn3) /
                     (pow(r1 * r1 + r2 * r2 + r3 * r3, 2.5))
                     * r[dof % 3] * r[other_dof % 3]
                     + sw * tn[dof % 3] * sn[other_dof % 3 ];
        ;

        double tmp = boundary_diag_tensors[pt_idx].get(dof % 3, other_dof % 3);
        boundary_diag_tensors[pt_idx].set(dof % 3, other_dof % 3, tmp + pot);

      }
    }
    curr_idx += hole.num_nodes * domain_dimension;
  }

  ki_Mat ident(3, 3);
  ident.eye(3);
  for (int pt_idx = 0; pt_idx < boundary->num_outer_nodes; pt_idx++) {
    boundary_diag_tensors[pt_idx] = (ident) - boundary_diag_tensors[pt_idx];
  }
  for (int pt_idx = boundary->num_outer_nodes; pt_idx < boundary->weights.size();
       pt_idx++) {
    boundary_diag_tensors[pt_idx] = - boundary_diag_tensors[pt_idx];
  }
}


ki_Mat Kernel::operator()(const std::vector<int>& tgt_inds,
                          const std::vector<int>& src_inds, bool forward, bool parallel) const {
  if (domain_dimension == 3) {
    return get_3d(tgt_inds, src_inds, forward, parallel);
  }
  ki_Mat ret(tgt_inds.size(), src_inds.size());
  int olda_ = tgt_inds.size();
  int threads = 1;
  if (parallel) {
    threads = 1;
  }
  #pragma omp parallel for num_threads(threads)
  for (int j = 0; j < src_inds.size(); j++) {
    int src_ind = src_inds[j];
    int j_point_index = src_ind / solution_dimension;
    int j_points_vec_index = j_point_index * domain_dimension;

    double sp[2] = {boundary_points_[j_points_vec_index],
                    boundary_points_[j_points_vec_index + 1]
                   };
    double sn[2] =  {boundary_normals(j_points_vec_index),
                     boundary_normals(j_points_vec_index + 1)
                    };
    double sw =  boundary_weights(j_point_index);
    double sc = boundary_curvatures(j_point_index);
    for (int i = 0; i < tgt_inds.size(); i++) {
      int tgt_ind = tgt_inds[i];
      int i_point_index = tgt_ind / solution_dimension;
      int i_points_vec_index = i_point_index * domain_dimension;

      double tp[2], tn[2];
      if (forward) {
        tp[0] = domain_points[i_points_vec_index];
        tp[1] = domain_points[i_points_vec_index + 1];
        tn[0] = 0;
        tn[1] = 0;
      } else {
        tp[0] = boundary_points_[i_points_vec_index];
        tp[1] = boundary_points_[i_points_vec_index + 1];
        tn[0] = boundary_normals(i_points_vec_index);
        tn[1] = boundary_normals(i_points_vec_index + 1);
      }
      double r[2] = {tp[0] - sp[0], tp[1] - sp[1]};
      switch (solution_dimension) {
        case 1:
          one_d_kern(i + olda_ * j, &ret, r[0], r[1], tn[0], tn[1],
                     sn[0], sn[1], sw, sc, forward);
          break;
        case 2:
          two_d_kern(i + olda_ * j, tgt_ind % 2, src_ind % 2, &ret, r[0], r[1],
                     tn[0], tn[1], sn[0], sn[1],  sw, sc);
          break;
      }
    }
  }
  return ret;
}


ki_Mat Kernel::get_3d(const std::vector<int>& tgt_inds,
                      const std::vector<int>& src_inds, bool forward, bool parallel) const {

  ki_Mat ret(tgt_inds.size(), src_inds.size());
  int olda_ = tgt_inds.size();
  int threads = 1;
  if (parallel) {
    threads = 1;
  }
  #pragma omp parallel for num_threads(threads)
  for (int j = 0; j < src_inds.size(); j++) {
    int src_ind = src_inds[j];
    int j_point_index = src_ind / solution_dimension;
    int j_points_vec_index = j_point_index * domain_dimension;

    double sp[3] = {boundary_points_[j_points_vec_index],
                    boundary_points_[j_points_vec_index + 1],
                    boundary_points_[j_points_vec_index + 2]
                   };
    double sn[3] =  {boundary_normals(j_points_vec_index),
                     boundary_normals(j_points_vec_index + 1),
                     boundary_normals(j_points_vec_index + 2)
                    };
    double sw =  boundary_weights(j_point_index);

    for (int i = 0; i < tgt_inds.size(); i++) {


      int tgt_ind = tgt_inds[i];

      int i_point_index = tgt_ind / solution_dimension;
      int i_points_vec_index = i_point_index * domain_dimension;
      double tp[3], tn[3];
      if (forward) {
        tp[0] = domain_points[i_points_vec_index];
        tp[1] = domain_points[i_points_vec_index + 1];
        tp[2] = domain_points[i_points_vec_index + 2];

        tn[0] = 0;
        tn[1] = 0;
        tn[2] = 0;
      } else {
        tp[0] = boundary_points_[i_points_vec_index];
        tp[1] = boundary_points_[i_points_vec_index + 1];
        tp[2] = boundary_points_[i_points_vec_index + 2];

        tn[0] = boundary_normals(i_points_vec_index);
        tn[1] = boundary_normals(i_points_vec_index + 1);
        tn[2] = boundary_normals(i_points_vec_index + 2);
      }
      double r[3] = {tp[0] - sp[0], tp[1] - sp[1], tp[2] - sp[2]};

      if (solution_dimension == 1) {
        three_d_laplace(i + olda_ * j, &ret, r[0], r[1], r[2],
                        boundary_diags[j_point_index],
                        sn[0], sn[1], sn[2], sw);
      } else {
        three_d_stokes(i + olda_ * j, tgt_ind % 3, src_ind % 3, &ret, r[0], r[1], r[2],
                       boundary_diag_tensors[j_point_index],
                       tn[0], tn[1], tn[2], sn[0], sn[1], sn[2], sw);

      }
    }
  }
  return ret;
}


ki_Mat Kernel::get_id_mat(const QuadTree* tree,
                          const QuadTreeNode* node) const {

  std::vector<int> active_box = node->dof_lists.active_box;
  // Grab all points inside the proxy circle which are outside the box
  std::vector<int> inner_circle, outside_box;

  // So if we're at level 2 or 1, we don't use the proxy trick
  // If at level 1, just grab active from neighbors
  if (node->level == 1) {
    for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
      if (level_node != node) {
        for (int matrix_index :
             level_node->dof_lists.active_box) {
          outside_box.push_back(matrix_index);
        }
      }
    }
    ki_Mat mat(2 * outside_box.size(), active_box.size());
    mat.set_submatrix(0, outside_box.size(), 0, active_box.size(),
                      (*this)(outside_box, active_box), false, true);
    mat.set_submatrix(outside_box.size(), 2 * outside_box.size(),
                      0, active_box.size(),
                      (*this)(active_box, outside_box), true, true);
    return mat;
  }

  // If at level 2, grab active from all on level, plus from leaves of level 1
  // Commented to allow exact updating, note this may affect acc if 
  // discretization is particularly heavy in low intrinsic dimension
  // if (node->level ==  && domain_dimension < 3) {
  //   for (QuadTreeNode* level_node : tree->levels[node->level]->nodes) {
  //     if (level_node != node) {
  //       for (int matrix_index :
  //            level_node->dof_lists.active_box) {
  //         outside_box.push_back(matrix_index);
  //       }
  //     }
  //   }
  //   for (int lvl = node->level - 1; lvl >= 0; lvl--) {
  //     for (QuadTreeNode* level_node : tree->levels[lvl]->nodes) {
  //       //TODO shouldnt this check for compression?
  //       if (level_node->is_leaf) {
  //         for (int matrix_index :
  //              level_node->dof_lists.original_box) {

  //           outside_box.push_back(matrix_index);
  //         }
  //       }
  //     }
  //   }

  //   std::vector<std::pair<double, int>> outboxtmp;
  //   for (int i = 0; i < outside_box.size(); i++) {
  //     int points_vec_index = domain_dimension * (outside_box[i] / solution_dimension);
  //     double x = boundary_points_[points_vec_index];
  //     double y = boundary_points_[points_vec_index + 1];
  //     double dist = 0.;
  //     for (int d = 0; d < domain_dimension; d++) {
  //       dist += pow(node->center[d] - boundary_points_[points_vec_index + d], 2);
  //     }
  //     dist = sqrt(dist);
  //     outboxtmp.push_back(std::pair<double, int>(dist,
  //                         outside_box[i]));
  //   }
  //   std::sort(outboxtmp.begin(),
  //             outboxtmp.end());
  //   outside_box.clear();
  //   for (int i = 0; i < outboxtmp.size(); i++) {
  //     outside_box.push_back(outboxtmp[i].second);


  //     // if ((std::abs(node->center[0]+0.6453917096) < 1e-5) &&
  //     //     (std::abs(node->center[1] -0.11820276) < 1e-5)) {
  //     //    int points_vec_index = domain_dimension * (outside_box[i] / solution_dimension);
  //     //   double x = boundary_points_[points_vec_index];
  //     // double y = boundary_points_[points_vec_index + 1];
  //     //   std::cout << "outside " << outboxtmp[i].first << " "<<outboxtmp[i].second<<" "<<x<<" "<<y<<
  //     // std::endl;
  //     // }
  //   }
  //   ki_Mat mat(2 * outside_box.size(), active_box.size());
  //   mat.set_submatrix(0, outside_box.size(), 0, active_box.size(),
  //                     (*this)(outside_box, active_box), false, true);
  //   mat.set_submatrix(outside_box.size(), 2 * outside_box.size(),
  //                     0, active_box.size(),
  //                     (*this)(active_box, outside_box), true, true);

  //   return mat;
  // }

  for (int matrix_index : node->dof_lists.near) {
    int point_index = matrix_index / solution_dimension;
    int points_vec_index = point_index * domain_dimension;
    double x = boundary_points_[points_vec_index];
    double y = boundary_points_[points_vec_index + 1];
    double dist = 0.;
    for (int d = 0; d < domain_dimension; d++) {
      dist += pow(node->center[d] - boundary_points_[points_vec_index + d], 2);
    }
    dist = sqrt(dist);
    if (dist < RADIUS_RATIO * node->side_length) {
      inner_circle.push_back(matrix_index);
    }
  }

  std::vector<std::pair<double, int>> innercirctmp;
  for (int i = 0; i < inner_circle.size(); i++) {
    int points_vec_index = domain_dimension * (inner_circle[i] /
                           solution_dimension);
    double x = boundary_points_[points_vec_index];
    double y = boundary_points_[points_vec_index + 1];
    double dist = 0.;
    for (int d = 0; d < domain_dimension; d++) {
      dist += pow(node->center[d] - boundary_points_[points_vec_index + d], 2);
    }
    dist = sqrt(dist);
    innercirctmp.push_back(std::pair<double, int>(dist,
                           inner_circle[i]));
  }
  std::sort(innercirctmp.begin(),
            innercirctmp.end());
  inner_circle.clear();
  for (int i = 0; i < innercirctmp.size(); i++) {
    inner_circle.push_back(innercirctmp[i].second);

  }

  int num_p_points = NUM_PROXY_POINTS;
  ki_Mat pxy = get_proxy_mat(node->center, NUM_PROXY_POINTS, node->side_length
                             * RADIUS_RATIO, tree, active_box);
  // Now all the matrices are gathered, put them into mat.
  ki_Mat mat(2 * inner_circle.size() + pxy.height(), active_box.size());
  mat.set_submatrix(0, inner_circle.size(),
                    0, active_box.size(), (*this)(inner_circle, active_box),
                    false, true);

  mat.set_submatrix(inner_circle.size(), 2 * inner_circle.size(),
                    0, active_box.size(), (*this)(active_box, inner_circle),
                    true, true);
  mat.set_submatrix(2 * inner_circle.size(),  pxy.height()
                    + 2 * inner_circle.size(), 0,
                    active_box.size(),  pxy, false, true);
  return mat;
}

// TODO(John) the general proxy stuff is just unreadable, replace with switch
// statements for kernels. It's not that bad. Or at least clean up the generality stuff

ki_Mat Kernel::get_proxy_mat(std::vector<double> center, int num_points,
                             double r, const QuadTree * tree,
                             const std::vector<int>& box_inds) const {
  if (domain_dimension == 3) {
    return get_proxy_mat3d(center, r, tree, box_inds);
  }

  // each row is a pxy point, cols are box dofs
  double pxy_w = 2.0 * M_PI * r / num_points;
  std::vector<double>  pxy_p, pxy_n;

  for (int i = 0; i < num_points; i++) {
    double ang = 2 * M_PI * i * (1.0 / num_points);
    for (int k = 1; k < 2; k++) {   // modify this for annulus proxy
      double eps = (k - 1) * 1;
      pxy_p.push_back(center[0] + (r + eps) * cos(ang));
      pxy_p.push_back(center[1] + (r + eps) * sin(ang));
      pxy_n.push_back(cos(ang));
      pxy_n.push_back(sin(ang));
    }
  }
  int lda = solution_dimension * pxy_p.size();
  ki_Mat ret(lda, box_inds.size());
  for (int j = 0; j < box_inds.size(); j++) {
    int box_ind = box_inds[j];
    int j_point_index = box_ind / solution_dimension;
    int j_points_vec_index = j_point_index * domain_dimension;
    double bp[2] = {boundary_points_[j_points_vec_index],
                    boundary_points_[j_points_vec_index + 1]
                   };
    double bn[2] =  {boundary_normals(j_points_vec_index),
                     boundary_normals(j_points_vec_index + 1)
                    };
    double bw =  boundary_weights(j_point_index);
    for (int i = 0; i < pxy_p.size(); i += 2) {
      double pp[2] = {pxy_p[i], pxy_p[i + 1]};
      double pn[2] = {pxy_n[i], pxy_n[i + 1]};
      double r[2] = {pp[0] - bp[0], pp[1] - bp[1]};
      switch (solution_dimension) {
        case 1:
          // box to pxy
          one_d_kern((i / 2) + lda * j, &ret, r[0], r[1],  pn[0], pn[1],
                     bn[0], bn[1], bw,  0);
          // pxy to box
          one_d_kern((pxy_p.size() / 2) + (i / 2) + lda * j, &ret,
                     -r[0], -r[1],  bn[0], bn[1], pn[0], pn[1], pxy_w, 0);
          break;
        case 2:
          for (int pxy_parity = 0; pxy_parity < 2; pxy_parity++) {
            // box to pxy
            two_d_kern(i + pxy_parity + lda * j, pxy_parity, box_ind % 2,
                       &ret, r[0], r[1],  pn[0], pn[1], bn[0], bn[1], bw, 0);
            // pxy to box
            two_d_kern(pxy_p.size() + i + pxy_parity + lda * j,
                       box_ind % 2, pxy_parity, &ret, -r[0], -r[1],
                       bn[0], bn[1], pn[0], pn[1], pxy_w, 0);
          }
          break;
      }
    }
  }
  return ret;
}


ki_Mat Kernel::get_proxy_mat3d(std::vector<double> center,
                               double r, const QuadTree * tree,
                               const std::vector<int>& box_inds) const {

  // each row is a pxy point, cols are box dofs
  // r *= 2;
  std::vector<double> pxy_p, pxy_n, pxy_w;
  for (int i = 0; i < pxy_thetas.size(); i++) {
    double theta = pxy_thetas[i];
    for (int j = 0; j < pxy_phis.size(); j++) {   // modify this for annulus proxy
      double phi = pxy_phis[j];
      pxy_p.push_back(center[0] + r * sin(phi) * cos(theta));
      pxy_p.push_back(center[1] + r * sin(phi) * sin(theta));
      pxy_p.push_back(center[2] + r * cos(phi));
      pxy_n.push_back(sin(phi) * cos(theta));
      pxy_n.push_back(sin(phi) * sin(theta));
      pxy_n.push_back(cos(phi));
      pxy_w.push_back(pow(r, 2) *pxy_phi_weights[j]*pxy_theta_weights[i]);
    }
  }

  int lda = 2 * solution_dimension * pxy_w.size();

  ki_Mat ret(lda, box_inds.size());
  for (int j = 0; j < box_inds.size(); j++) {
    int box_ind = box_inds[j];
    int j_point_index = box_ind / solution_dimension;
    int j_points_vec_index = j_point_index * domain_dimension;

    double bp[3] = {boundary_points_[j_points_vec_index],
                    boundary_points_[j_points_vec_index + 1],
                    boundary_points_[j_points_vec_index + 2]
                   };
    double bn[3] = {boundary_normals(j_points_vec_index),
                    boundary_normals(j_points_vec_index + 1),
                    boundary_normals(j_points_vec_index + 2)
                   };

    double bw =  boundary_weights(j_point_index);
    for (int i = 0; i < pxy_p.size(); i += domain_dimension) {
      double pp[3] = {pxy_p[i],
                      pxy_p[i + 1],
                      pxy_p[i + 2]
                     };
      double pn[3] = {pxy_n[i],
                      pxy_n[i + 1],
                      pxy_n[i + 2]
                     };

      double r[3] = {pp[0] - bp[0], pp[1] - bp[1], pp[2] - bp[2]};
      if (solution_dimension == 1) {
        three_d_laplace((i / domain_dimension) + lda * j, &ret,
                        r[0], r[1], r[2], 0.,
                        bn[0], bn[1], bn[2],
                        bw);
        three_d_laplace((pxy_p.size() / domain_dimension) + (i / domain_dimension) + lda
                        * j, &ret,
                        -r[0], -r[1], -r[2], 0.,
                        pn[0], pn[1], pn[2],
                        pxy_w[i / domain_dimension]);
      } else {
        for (int pxy_parity = 0; pxy_parity < 3; pxy_parity++) {
          // box to pxy
          three_d_stokes(i + pxy_parity + lda * j, pxy_parity, box_ind % 3,
                         &ret, r[0], r[1], r[2], ki_Mat(3, 3), pn[0], pn[1], pn[2],
                         bn[0], bn[1], bn[2], bw);

          // pxy to box
          three_d_stokes(pxy_p.size() + i + pxy_parity + lda * j,
                         box_ind % 3, pxy_parity, &ret, -r[0], -r[1], -r[2], ki_Mat(3, 3),
                         bn[0], bn[1], bn[2], pn[0], pn[1], pn[2], pxy_w[i / domain_dimension]);
        }
      }

    }
  }
  return ret;
}


ki_Mat Kernel::forward() const {
  int tgt = solution_dimension * (domain_points.size() / domain_dimension);
  int src = solution_dimension * (boundary_points_.size() / domain_dimension);
  std::vector<int> tgt_inds(tgt), src_inds(src);
  for (int i = 0; i < tgt; i++) tgt_inds[i] = i;
  for (int j = 0; j < src; j++) src_inds[j] = j;
  return (*this)(tgt_inds, src_inds, true, true);
}

}  // namespace kern_interp
