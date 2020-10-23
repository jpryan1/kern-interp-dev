// Copyright 2019 John Paul Ryan
#include <string.h>
#include <omp.h>
#include <string>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include "kern_interp/ki_mat.h"

namespace kern_interp {



double ki_Mat::condition_number() const {
  // TODO(John) this should be done without using LU factorization if possible?
  ki_Mat X_copy = *this;
  std::vector<lapack_int> ipiv(height_);
  memset(&ipiv[0], 0, height_ * sizeof(lapack_int));
  double one_norm =  X_copy.one_norm();
  LAPACKE_dgetrf(LAPACK_COL_MAJOR, X_copy.height_, X_copy.width_, X_copy.mat,
                 X_copy.lda_, &ipiv[0]);
  double rcond;
  LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', X_copy.height(), X_copy.mat,
                 X_copy.lda_, one_norm, &rcond);
  return 1.0 / rcond;
}



ki_Mat::ki_Mat() {
  mat       = NULL;
  lda_      = 0;
  height_   = 0;
  width_    = 0;
}


ki_Mat::~ki_Mat() {
  if (mat) delete[] mat;
}

// Copy constructor
ki_Mat::ki_Mat(const ki_Mat &copy) {
  lda_ = copy.lda_;
  height_ = copy.height_;
  width_ = copy.width_;
  mat = new double[height_ * width_];
  memcpy(mat, copy.mat, width_ * height_ * sizeof(double));
}

// Copy assignment
ki_Mat& ki_Mat::operator=(const ki_Mat& copy) {
  if (height_ != copy.height_ || width_ != copy.width_ || lda_ != copy.lda_) {
    if (mat) delete[] mat;
    lda_ = copy.lda_;
    height_ = copy.height_;
    width_ = copy.width_;
    mat = new double[height_ * width_];
  }
  memcpy(mat, copy.mat, width_ * height_ * sizeof(double));
  return *this;
}


// Move constructor
ki_Mat::ki_Mat(ki_Mat&& move) {
  lda_ = move.lda_;
  height_ = move.height_;
  width_ = move.width_;
  mat = move.mat;
  move.mat = nullptr;
}


// Move assignment
ki_Mat& ki_Mat::operator= (ki_Mat&& move) {
  if (mat) delete[] mat;
  if (height_ != move.height_ || width_ != move.width_ || lda_ != move.lda_) {
    lda_ = move.lda_;
    height_ = move.height_;
    width_ = move.width_;
  }
  mat = move.mat;
  move.mat = nullptr;
  return *this;
}


ki_Mat::ki_Mat(int h, int w) {
  lda_      = h;
  height_   = h;
  width_    = w;
  mat       = new double[height_ * width_];
  memset(mat, 0, height_ * width_ * sizeof(double));
}


double ki_Mat::get(int i, int j) const {

  assert(i < height_ && j < width_ && mat != NULL);
  return mat[i + lda_ * j];
}


void ki_Mat::set(int i, int j, double a) {
  assert(i < height_ && j < width_ && mat != NULL);
  mat[i + lda_ * j] = a;
}


void ki_Mat::addset(int i, int j, double a) {
  assert(i < height_ && j < width_ && mat != NULL);
  mat[i + lda_ * j] += a;
}


void ki_Mat::set_submatrix(const std::vector<int>& I_,
                           const std::vector<int>& J_,
                           const ki_Mat& A, bool transpose_A) {
  if (transpose_A) {
    assert(I_.size() == A.width_ && J_.size() == A.height_);
    for (int i = 0; i < I_.size(); i++) {
      for (int j = 0; j < J_.size(); j++) {
        set(I_[i], J_[j], A.get(j, i));
      }
    }
  } else {
    assert(I_.size() == A.height_ && J_.size() == A.width_);
    for (int i = 0; i < I_.size(); i++) {
      for (int j = 0; j < J_.size(); j++) {
        set(I_[i], J_[j], A.get(i, j));
      }
    }
  }
}


void ki_Mat::set_submatrix(int row_s, int row_e,
                           int col_s, int col_e,
                           const ki_Mat& A, bool transpose_A, bool timing) {
  if (transpose_A) {
    for (int i = 0; i < row_e - row_s; i++) {
      for (int j = 0; j < col_e - col_s; j++) {
        set(i + row_s, j + col_s, A.get(j, i));
      }
    }
    assert(row_e - row_s == A.width_ && col_e - col_s == A.height_);
  } else {
    assert(row_e - row_s == A.height_ && col_e - col_s == A.width_);
    for (int j = 0; j < col_e - col_s; j++) {
      memcpy(&(mat[row_s + lda_ * (j + col_s)]), &(A.mat[A.lda_ * j]),
             (row_e - row_s)*sizeof(double));
    }
  }
}


void ki_Mat::set_submatrix(const std::vector<int>& I_,
                           int col_s, int col_e,
                           const ki_Mat& A, bool transpose_A) {
  if (transpose_A) {
    assert(I_.size() == A.width_ &&  col_e - col_s  == A.height_);
    for (int i = 0; i < I_.size(); i++) {
      for (int j = 0; j < col_e - col_s; j++) {
        set(I_[i], j + col_s, A.get(j, i));
      }
    }
  } else {
    assert(I_.size() == A.height_ &&  col_e - col_s  == A.width_);
    for (int i = 0; i < I_.size(); i++) {
      for (int j = 0; j < col_e - col_s; j++) {
        set(I_[i], j + col_s, A.get(i, j));
      }
    }
  }
}


void ki_Mat::set_submatrix(int row_s, int row_e,
                           const std::vector<int>& J_,
                           const ki_Mat& A, bool transpose_A) {
  if (transpose_A) {
    assert(row_e - row_s == A.width_ && J_.size() == A.height_);
    for (int i = 0; i < row_e - row_s; i++) {
      for (int j = 0; j < J_.size(); j++) {
        set(i + row_s, J_[j], A.get(j, i));
      }
    }
  } else {
    assert(row_e - row_s == A.height_ && J_.size() == A.width_);
    for (int j = 0; j < J_.size(); j++) {
      memcpy(&(mat[row_s + lda_ * J_[j] ]), &(A.mat[A.lda_ * j]),
             (row_e - row_s)*sizeof(double));
    }
  }
}


ki_Mat ki_Mat::operator()(int row_s, int row_e,
                          int col_s, int col_e) const {
  ki_Mat submatrix(row_e - row_s, col_e - col_s);
  for (int i = 0; i < row_e - row_s; i++) {
    for (int j = 0; j < col_e - col_s; j++) {
      submatrix.set(i, j, this->get(i + row_s, j + col_s));
    }
  }
  return submatrix;
}


void ki_Mat::transpose_into(ki_Mat* transpose) const {
  if (height_ != transpose->width_ || width_ != transpose->height_) {
    if (transpose->mat) delete[] transpose->mat;
    transpose->lda_    = width_;
    transpose->height_ = width_;
    transpose->width_  = height_;
    transpose->mat     = new double[height_ * width_];
  }
  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      transpose->mat[j + i * width_] = mat[i + lda_ * j];
    }
  }
}


void ki_Mat::eye(int n) {
  if (width_ != n || height_ != n || lda_ != n) {
    if (mat) delete[] mat;
    lda_    = n;
    height_ = n;
    width_  = n;
    mat     = new double[height_ * width_];
  }
  for (int i = 0; i < n; i++) {
    set(i, i, 1.0);
  }
}


int ki_Mat::height() const {
  return height_;
}


int ki_Mat::width() const {
  return width_;
}


ki_Mat& ki_Mat::operator-=(const ki_Mat& o) {
  assert(o.height_ == height_ && o.width_ == width_);

  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      mat[i + lda_ * j] =  mat[i + lda_ * j] - o. mat[i + lda_ * j];
    }
  }
  return *this;
}


ki_Mat& ki_Mat::operator+=(const ki_Mat& o) {
  assert(o.height_ == height_ && o.width_ == width_);
  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      mat[i + lda_ * j] =  mat[i + lda_ * j] + o.mat[i + lda_ * j];
    }
  }
  return *this;
}


ki_Mat& ki_Mat::operator*=(double o) {
  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      mat[i + lda_ * j] =  mat[i + lda_ * j] * o;
    }
  }
  return *this;
}

ki_Mat& ki_Mat::operator/=(double o) {
  assert(std::abs(o) > 1e-8 && "Error: divide matrix by 0.");
  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      mat[i + lda_ * j] =  mat[i + lda_ * j] / o;
    }
  }
  return *this;
}


ki_Mat ki_Mat::operator-() const {
  ki_Mat result(height_, width_);
  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      result.set(i, j, -this->get(i, j));
    }
  }
  return result;
}


ki_Mat ki_Mat::operator-(const ki_Mat& o) const {
  assert(o.height_ == height_ && o.width_ == width_);
  ki_Mat result(height_, width_);
  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      result.set(i, j, this->get(i, j) - o.get(i, j));
    }
  }
  return result;
}


ki_Mat ki_Mat::operator+(const ki_Mat& o) const {
  assert(o.height_ == height_ && o.width_ == width_);
  ki_Mat sum(height_, width_);
  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      sum.set(i, j, this->get(i, j) + o.get(i, j));
    }
  }
  return sum;
}


ki_Mat ki_Mat::operator*(const ki_Mat& o) const {
  ki_Mat result(height_, o.width_);

  cblas_dgemm(CblasColMajor, NORMAL, NORMAL, height_, o.width_,
              width_, 1., mat, height_, o.mat, o.height_, 0.,
              result.mat, height_);
  return result;
}


ki_Mat ki_Mat::operator*(double o) const {
  ki_Mat result(height_, width_);
  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      result.set(i, j, this->get(i, j) *o);
    }
  }
  return result;
}


ki_Mat ki_Mat::operator/(double o) const {
  assert(std::abs(o) > 1e-8 && "Error: divide matrix by 0.");
  ki_Mat result(height_, width_);
  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      result.set(i, j, this->get(i, j) / o);
    }
  }
  return result;
}


ki_Mat ki_Mat::operator()(const std::vector<int>& I_,
                          const std::vector<int>& J_) const {
  ki_Mat ret(I_.size(), J_.size());
  int olda_ = I_.size();
  for (int i = 0; i < I_.size(); i++) {
    for (int j = 0; j < J_.size(); j++) {
      assert(I_[i] < height() && J_[j] < width());
      ret.mat[i + olda_ * j] = get(I_[i], J_[j]);
    }
  }
  return ret;
}


ki_Mat ki_Mat::operator()(const std::vector<int>& I_,
                          int col_s, int col_e) const {
  ki_Mat ret(I_.size(), col_e - col_s);
  int olda_ = I_.size();
  for (int i = 0; i < I_.size(); i++) {
    for (int j = 0; j < col_e - col_s; j++) {
      assert(I_[i] < height() && col_s + j < width());
      ret.mat[i + olda_ * j] = get(I_[i], col_s + j);
    }
  }
  return ret;
}


ki_Mat ki_Mat::operator()(int row_s, int row_e,
                          const std::vector<int>& J_) const {
  ki_Mat ret(row_e - row_s, J_.size());
  int olda_ = row_e - row_s;
  for (int i = 0; i < row_e - row_s; i++) {
    for (int j = 0; j < J_.size(); j++) {
      assert(row_s + i < height() && J_[j] < width());
      ret.mat[i + olda_ * j] = get(row_s + i, J_[j]);
    }
  }
  return ret;
}


double ki_Mat::frob_norm() const {
  double sum = 0;
  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      sum += pow(get(i, j), 2);
    }
  }
  return sqrt(sum);
}

double ki_Mat::max_entry_magnitude() const {
  double max = 0;
  for (int i = 0; i < height_; i++) {
    for (int j = 0; j < width_; j++) {
      max = std::max(max, std::abs(get(i, j)));
    }
  }
  return max;
}


double ki_Mat::one_norm() const {
  double top = 0;
  for (int j = 0; j < width_; j++) {
    double sum = 0;
    for (int i = 0; i < height_; i++) {
      sum += fabs(get(i, j));
    }
    if (sum > top) {
      top = sum;
    }
  }
  return top;
}


double ki_Mat::vec_two_norm() const {
  assert(width() == 1);
  return frob_norm();
}


void ki_Mat::LU_factorize(ki_Mat* K_LU, std::vector<lapack_int>* piv) const {
  *K_LU = *this;
  *piv = std::vector<lapack_int>(height_);

  int status = LAPACKE_dgetrf(LAPACK_COL_MAJOR, K_LU->height_, K_LU->width_,
                              K_LU->mat,
                              K_LU->lda_, &(*piv)[0]);
  assert(status == 0);
}

void ki_Mat::left_multiply_inverse(const ki_Mat& K, ki_Mat* U) const {
  // X^-1K = U
  // aka, XU = K

  ki_Mat X_copy;
  *U = K;
  std::vector<lapack_int> ipiv;
  LU_factorize(&X_copy, &ipiv);

  int status = LAPACKE_dgetrs(LAPACK_COL_MAJOR , 'N' , X_copy.height_ ,
                              U->width_ , X_copy.mat , X_copy.lda_ ,
                              &ipiv[0] , U->mat, U->lda_);
  assert(status == 0);
}


void ki_Mat::right_multiply_inverse(const ki_Mat& K, ki_Mat* L) const {
  ki_Mat K_copy = K.transpose();
  std::vector<lapack_int> ipiv;
  // KX^-1 = L
  // aka X_T L^T = K^T
  ki_Mat X_copy;
  LU_factorize(&X_copy, &ipiv);

  int err2 = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'T', X_copy.height_,
                            K_copy.width_, X_copy.mat, X_copy.lda_, &ipiv[0],
                            K_copy.mat, K_copy.lda_);
  assert(err2 == 0);
  *L = K_copy.transpose();
}


void ki_Mat::left_multiply_inverse(const ki_Mat& K,
                                   const std::vector<lapack_int>& piv,
                                   ki_Mat* U) const {
  // X^-1K = U
  // aka, XU = K
  ki_Mat X_copy = *this;
  *U = K;

  int status = LAPACKE_dgetrs(LAPACK_COL_MAJOR , 'N' , this->height_ ,
                              U->width_ , this->mat , this->lda_ ,
                              &piv[0] , U->mat, U->lda_);
  assert(status == 0);
}


void ki_Mat::right_multiply_inverse(const ki_Mat& K,
                                    const std::vector<lapack_int>& piv,
                                    ki_Mat* L) const {
  ki_Mat K_copy = K.transpose();
  // KX^-1 = L
  // aka X_T L^T = K^T
  int err2 = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'T', this->height_,
                            K_copy.width_, this->mat, this->lda_, &piv[0],
                            K_copy.mat, K_copy.lda_);
  assert(err2 == 0);
  *L = K_copy.transpose();
}

ki_Mat ki_Mat::transpose() const {
  ki_Mat transpose(width(), height());
  transpose_into(&transpose);
  return transpose;
}


// Performs interpolative decomposition, and returns number of skeleton columns.
// Takes double /tol/, tolerance factorfor error in CPQR factorization.
// Populates /p/ with permutation, Z with linear transformation.
int ki_Mat::id(std::vector<int>* p, ki_Mat* Z, double tol) const {
  ki_Mat cpy = *this;
  assert(height() != 0 &&  width() != 0);
  std::vector<lapack_int> pvt(width_);
  memset(&pvt[0], 0, width_ * sizeof(lapack_int));

  // /tau/ will contain an output from dgeqp3 that we don't need.
  std::vector<double> tau(width_);

  int info1 = LAPACKE_dgeqp3(CblasColMajor, height_, width_, cpy.mat, lda_,
                             &pvt[0], &tau[0]);

  assert(info1 == 0);
  assert(width_ < height_
         && "Width>=Height in ID, maybe this should happen, but probably a bug");
  int skel = 0;
  double thresh = fabs(tol * cpy.get(0, 0));
  int smaller = std::min(height_, width_);
  for (int i = 1; i < smaller; i++) {
    // check if R_{i,i} / R_{0,0} < tol
    // std::cout << i << " out of " << width_ << std::endl;
    if (fabs(cpy.get(i, i)) < thresh) {
      skel = i;
      break;
    }
  }
  if (skel == 0) {
    // no compression to be done :/
    return 0;
  }
  *p = std::vector<int>(width_);

  for (int i = 0; i < width_; i++) {
    (*p)[i] = pvt[i] - 1;
  }

  int redund = width_ - skel;
  // set Z to be R_11^-1 R_12. Note 'U' (above diagonal) part of cp.mat
  // is the R matrix from dgeqp3.
  int info2 = LAPACKE_dtrtrs(CblasColMajor, 'U', 'N', 'N', skel, redund,
                             cpy.mat, cpy.lda_, cpy.mat + cpy.lda_ * skel,
                             cpy.lda_);
  assert(info2 == 0);
  *Z = cpy(0, skel, skel, skel + redund);
  return skel;
}


std::vector<double> ki_Mat::real_eigenvalues() {
  std::vector<double> eigvs;
  double *eigs = new double[width_];
  double *imags = new double[width_];
  int info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', width_, mat,
                           width_, eigs, imags,
                           nullptr, 1, nullptr, 1);
  assert(info == 0);
  for (int i = 0; i < width_; i++) {
    if (fabs(imags[i]) < 1e-14) {
      eigvs.push_back(eigs[i]);
    }
  }
  delete[] imags;
  delete[] eigs;
  return eigvs;
}

ki_Mat ki_Mat::rand_vec(int height) {
  ki_Mat random(height, 1);
  for (int i = 0; i < height; i++) {
    random.set(i, 0, rand() / (0. + RAND_MAX));
  }
  return random;
}


}  // namespace kern_interp
