#include <geometry.hh>

struct PeriodicCurve
{
  size_t p;                     // degree
  size_t n;                     // n + 1 = cp.size()
  Geometry::DoubleVector knots; // n + 2 values
  Geometry::PointVector cp;     // n + 1 points

  double getKnot(int i) const;
  Geometry::Point3D getCP(int i) const;
  void setCP(int i, const Geometry::Point3D &p);

  size_t findSpan(double u, size_t *multi = nullptr) const;
  void basisFunctionDerivatives(size_t i, double u, size_t d, Geometry::DoubleMatrix &der) const;
  Geometry::Point3D derivatives(double u, size_t d, Geometry::VectorVector &der) const;
  PeriodicCurve insertKnot(double u, size_t r) const;

  Geometry::BSCurve clamp(double u) const;
};
