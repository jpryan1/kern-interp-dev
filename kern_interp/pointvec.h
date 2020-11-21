// Copyright 2019 John Paul Ryan
#ifndef KERN_INTERP_POINTVEC_H_
#define KERN_INTERP_POINTVEC_H_

#include <vector>


namespace kern_interp {


struct PointVec {
  std::vector<double> a;

  PointVec() {}
  PointVec(double m, double n);
  PointVec(double x, double y, double z);
  explicit PointVec(const std::vector<double>& b);

  double norm() const;
  double dot(const PointVec& o) const;
  PointVec operator-(const PointVec &o) const;
  PointVec operator+(const PointVec &o) const;
  PointVec operator*(const double d) const;
};

}  // namespace kern_interp

#endif  // KERN_INTERP_POINTVEC_H_
