// Copyright 2019 John Paul Ryan
#ifndef KERN_INTERP_KI_MAT_H_
#define KERN_INTERP_KI_MAT_H_

#include <cblas.h>
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>
#include <vector>
#include <string>

#define NORMAL CblasNoTrans
#define TRANSPOSE CblasTrans

namespace kern_interp {


struct ki_Mat {
  // storage is column major, so by default lda is the height.
  double *mat;

  int lda_, height_, width_;
  ki_Mat();
  ~ki_Mat();

  ki_Mat(int h, int w);

  // Copy constructor
  ki_Mat(const ki_Mat &o);
  // Copy assignment
  ki_Mat& operator=(const ki_Mat& copy);
  // Move constructor
  ki_Mat(ki_Mat&& move);
  // Move assignment
  ki_Mat& operator=(ki_Mat&& move);

  double get(int i, int j) const;
  void set(int i, int j, double a);
  void addset(int i, int j, double a);
  void set_submatrix(const std::vector<int>& I_,
                     const std::vector<int>& J_, const ki_Mat& A,
                     bool transpose_A = false);
  void set_submatrix(int row_s, int row_e, int col_s,
                     int col_e,
                     const ki_Mat& A, bool transpose_A = false,
                     bool timing = false);
  void set_submatrix(const std::vector<int>& I_, int col_s,
                     int col_e,
                     const ki_Mat& A, bool transpose_A = false);
  void set_submatrix(int row_s, int row_e,
                     const std::vector<int>& J_,
                     const ki_Mat& A, bool transpose_A = false);

  void transpose_into(ki_Mat* transpose) const;
  void eye(int n);
  ki_Mat transpose() const;

  int height() const;
  int width() const;

  ki_Mat& operator-=(const ki_Mat& o);
  ki_Mat& operator+=(const ki_Mat& o);
  ki_Mat& operator*=(double o);
  ki_Mat& operator/=(double o);

  ki_Mat operator-() const;
  ki_Mat operator-(const ki_Mat& o) const;
  ki_Mat operator+(const ki_Mat& o) const;
  ki_Mat operator*(const ki_Mat& o) const;
  ki_Mat operator*(double o) const;
  ki_Mat operator/(double o) const;

  ki_Mat operator()(const std::vector<int>& I_,
                    const std::vector<int>& J_) const;
  ki_Mat operator()(const std::vector<int>& I_,
                    int col_s, int col_e) const;
  ki_Mat operator()(int row_s, int row_e,
                    const std::vector<int>& J_) const;
  ki_Mat operator()(int row_s, int row_e, int col_s,
                    int col_e) const;


  double one_norm() const;
  double vec_two_norm() const;
  double frob_norm() const;
  double max_entry_magnitude() const;

  void LU_factorize(ki_Mat* K_LU, std::vector<lapack_int>* piv) const;
  void left_multiply_inverse(const ki_Mat& K, ki_Mat* U) const;
  void right_multiply_inverse(const ki_Mat& K, ki_Mat* L) const;
  void left_multiply_inverse(const ki_Mat& K,
                             const std::vector<lapack_int>& piv,
                             ki_Mat* U) const;
  void right_multiply_inverse(const ki_Mat& K,
                              const std::vector<lapack_int>& piv,
                              ki_Mat* L) const;

  int id(std::vector<int>* p, ki_Mat* Z, double tol) const;
  std::vector<double> real_eigenvalues();
  double condition_number() const;

  static ki_Mat rand_vec(int height);
};

}  // namespace kern_interp

#endif  // KERN_INTERP_KI_MAT_H_
