#include "clamp.hh"

#include <algorithm>
#include <cassert>
#include <limits>

using namespace Geometry;

double PeriodicCurve::getKnot(int i) const {
  if (i >= (int)n + 2)
    return knots.back() + (knots[i-n-1] - knots.front());
  if (i >= 0)
    return knots[i];
  return knots.front() - (knots.back() - knots[n+i+1]);
}

Point3D PeriodicCurve::getCP(int i) const {
  return cp[((i % ((int)n + 1)) + (n + 1)) % ((int)n + 1)];
}

void PeriodicCurve::setCP(int i, const Point3D &p) {
  cp[((i % ((int)n + 1)) + (n + 1)) % ((int)n + 1)] = p;
}

size_t PeriodicCurve::findSpan(double u, size_t *multi) const
{
  auto range = std::equal_range(knots.begin(), knots.end(), u);
  if (multi)
    *multi = range.second - range.first;
  return (range.second - knots.begin()) - 1;
}

void PeriodicCurve::basisFunctionDerivatives(size_t i, double u, size_t d, DoubleMatrix &der) const
{
  der.clear(); der.resize(d + 1);
  DoubleVector left(p + 1), right(p + 1), a[2];
  a[0].resize(p + 1); a[1].resize(p + 1);
  DoubleMatrix ndu(p + 1);
  ndu[0].resize(p + 1); ndu[0][0] = 1.0;
  for (size_t j = 1; j <= p; ++j) {
    ndu[j].resize(p + 1);
    left[j] = u - getKnot(i+1-j);
    right[j] = getKnot(i+j) - u;
    double saved = 0.0;
    for (size_t r = 0; r < j; ++r) {
      // lower triangle
      ndu[j][r] = right[r+1] + left[j-r];
      double tmp = ndu[r][j-1] / ndu[j][r];
      // upper triangle
      ndu[r][j] = saved + tmp * right[r+1];
      saved = tmp * left[j-r];
    }
    ndu[j][j] = saved;
  }
  for (size_t j = 0; j <= p; ++j)
    der[0].push_back(ndu[j][p]);
  for (size_t r = 0; r <= p; ++r) {
    size_t s1 = 0, s2 = 1;
    a[0][0] = 1.0;
    for (size_t k = 1; k <= d; ++k) {
      double dd = 0.0;
      int rk = r - k;
      int pk = p - k;
      if (r >= k) {
        a[s2][0] = a[s1][0] / ndu[pk+1][rk];
        dd = a[s2][0] * ndu[rk][pk];
      }
      size_t j1 = rk >= -1 ? 1 : -rk;
      size_t j2 = (int)r - 1 <= pk ? k - 1 : p - r;
      for (size_t j = j1; j <= j2; ++j) {
        a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j];
        dd += a[s2][j] * ndu[rk+j][pk];
      }
      if (r <= (size_t)pk) {
        a[s2][k] = -a[s1][k-1] / ndu[pk+1][r];
        dd += a[s2][k] * ndu[r][pk];
      }
      der[k].push_back(dd);
      std::swap(s1, s2);
    }
  }
  size_t r = p;
  for (size_t k = 1; k <= d; ++k) {
    for (size_t j = 0; j <= p; ++j)
      der[k][j] *= r;
    r *= p - k;
  }
}

Point3D PeriodicCurve::derivatives(double u, size_t d, VectorVector &der) const
{
  if (u >= knots.back()) // >= instead of == because of numerical problems
    u = knots.front();
  size_t du = std::min(d, p);
  der.clear();
  size_t span = findSpan(u);
  DoubleMatrix nder; basisFunctionDerivatives(span, u, du, nder);
  for (size_t k = 0; k <= du; ++k) {
    der.emplace_back(0.0, 0.0, 0.0);
    for (size_t j = 0; j <= p; ++j)
      der[k] += getCP(span-p+j) * nder[k][j];
  }
  for (size_t k = p + 1; k <= d; ++k)
    der.emplace_back(0.0, 0.0, 0.0);
  return der[0];
}

PeriodicCurve PeriodicCurve::insertKnot(double u, size_t r) const {
  size_t s;
  size_t k = findSpan(u, &s);
  if (s >= p)
    return *this;
  r = std::min(r, p - s);

  PeriodicCurve result;
  result.p = p;
  result.n = n + r;
  result.knots.reserve(knots.size() + r);
  std::copy_n(knots.begin(), k + 1, std::back_inserter(result.knots));
  std::fill_n(std::back_inserter(result.knots), r, u);
  std::copy(knots.begin() + k + 1, knots.end(), std::back_inserter(result.knots));

  result.cp.resize(cp.size() + r);
  // Interesting interval: [k-p+1, k-s)
  if (k + 1 >= p)
    for (size_t i = 0; i < k + 1 - p; ++i)
      result.setCP(i, getCP(i));
  for (size_t i = k - s; i < cp.size(); ++i)
    result.setCP(i + r, getCP(i));

  PointVector tmp; tmp.reserve(p - s + 1);
  for (size_t i = 0; i < p - s + 1; ++i)
    tmp.push_back(getCP(i + k - p));

  int L = k - p + 1;
  for (size_t j = 1; j <= r; ++j, ++L) {
    for (size_t i = 0; i <= p - j - s; ++i) {
      double alpha = (u - getKnot(L+i)) / (getKnot(i+k+1) - getKnot(L+i));
      tmp[i] = tmp[i+1] * alpha + tmp[i] * (1.0 - alpha);
    }
    result.setCP(L, tmp[0]);
    result.setCP(k + r - j - s, tmp[p-j-s]);
  }
  if (p > s + r + 1)
    for (size_t i = 0; i < p - s - 1 - r; ++i)
      result.setCP(L + i, tmp[i+1]);

  return result;
}

BSCurve PeriodicCurve::clamp(double u) const {
  auto c = insertKnot(u, p);
  int span = c.findSpan(u);

  DoubleVector k;
  std::fill_n(std::back_inserter(k), p, u);
  for (size_t i = 0; i < c.knots.size(); ++i)
    k.push_back(c.getKnot(span + i));
  auto last = k.back();
  std::fill_n(std::back_inserter(k), p, last);

  PointVector cpts;
  for (size_t i = 0; i < k.size() - p - 1; ++i)
    cpts.push_back(c.getCP(span - p + i));

  BSCurve result(p, k, cpts);
  result.normalize();
  return result;
}
