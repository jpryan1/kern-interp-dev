// Copyright 2019 John Paul Ryan
#ifndef KERN_INTERP_SKEL_FACTORIZATION_SKEL_FACTORIZATION_H_
#define KERN_INTERP_SKEL_FACTORIZATION_SKEL_FACTORIZATION_H_

#include <atomic>
#include <set>
#include <vector>
#include "kern_interp/ki_mat.h"
#include "kern_interp/quadtree/quadtree.h"
#include "kern_interp/kernel/kernel.h"

#define MIN_DOFS_TO_COMPRESS 64
#define NODE_CAP INFINITY
#define LEVEL_CAP INFINITY

namespace kern_interp {

// Perhaps putting these in a class isn't best, maybe just a namespace.
// Or maybe make everything static.
class SkelFactorization {
 public:

  double id_tol;
  int fact_threads;

  SkelFactorization() {}
  SkelFactorization(double id_tol, int fact_threads);
  ~SkelFactorization() {}
  void decouple(const Kernel& K, QuadTreeNode* node);
  void decouple(const Kernel& kernel, MidLevelNode* node);
  int id_compress(const Kernel& K, const QuadTree* tree,
                  QuadTreeNode* node);
  int id_compress(const Kernel& kernel, const QuadTree* tree,
                  MidLevelNode* node) ;

  void skeletonize(const Kernel& K, QuadTree* tree);

  void make_id_mat(const Kernel& K, ki_Mat* pxy, const QuadTree* tree,
                   const QuadTreeNode* node);
  ki_Mat make_proxy_mat(const Kernel& kernel, double cntr_x,
                        double cntr_y, double r,
                        const QuadTree* tree,
                        const std::vector<int>& box_indices);

  void apply_sweep_matrix(const ki_Mat& mat, ki_Mat* vec,
                          const std::vector<int>& a,
                          const std::vector<int>& b,
                          bool transpose) const;
  void apply_diag_matrix(const ki_Mat& mat, ki_Mat* vec,
                         const std::vector<int>& range) const;
  void apply_diag_inv_matrix(const ki_Mat& mat,
                             const std::vector<lapack_int>& piv, ki_Mat* vec,
                             const std::vector<int>& range) const;

  void solve(const QuadTree& quadtree, ki_Mat* x, const ki_Mat& b) const;
  void multiply_connected_solve(const QuadTree& quadtree, ki_Mat* x,
                                ki_Mat* alpha, const ki_Mat& b) const;
};

}  // namespace kern_interp

#endif  // KERN_INTERP_SKEL_FACTORIZATION_SKEL_FACTORIZATION_H_
