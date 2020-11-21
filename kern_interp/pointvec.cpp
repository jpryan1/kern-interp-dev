// Copyright 2019 John Paul Ryan
#include <string>
#include <cmath>
#include "kern_interp/pointvec.h"

namespace kern_interp {

PointVec::PointVec(double m, double n) {
  a.push_back(m);
  a.push_back(n);
}


PointVec::PointVec(double x, double y, double z) {
  a.push_back(x);
  a.push_back(y);
  a.push_back(z);
}


PointVec::PointVec(const std::vector<double>& b) {
  a = b;
}


double PointVec::norm() const {
  double tot = 0.;
  for (int i = 0; i < a.size(); i++) tot += pow(a[i], 2);
  return sqrt(tot);
}


double PointVec::dot(const PointVec& o) const {
  double tot = 0.;
  for (int i = 0; i < a.size(); i++) tot += a[i] * o.a[i];
  return tot;
}


PointVec PointVec::operator-(const PointVec &o) const {
  std::vector<double> out;
  for (int i = 0; i < a.size(); i++) out.push_back(a[i] - o.a[i]);
  return PointVec(out);
}

PointVec PointVec::operator+(const PointVec &o) const {
  std::vector<double> out;
  for (int i = 0; i < a.size(); i++) out.push_back(a[i] + o.a[i]);
  return PointVec(out);
}


PointVec PointVec::operator*(const double d) const {
  std::vector<double> out;
  for (int i = 0; i < a.size(); i++) out.push_back(a[i] * d);
  return PointVec(out);
}

}  // namespace kern_interp
