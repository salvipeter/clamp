#include <fstream>

#include "clamp.hh"

using namespace Geometry;

void writeBSCurve(const BSCurve &pc, std::string filename, size_t resolution = 100) {
  std::ofstream f(filename);
  for (size_t i = 0; i < resolution; ++i) {
    double u = (double)i / (resolution - 1);
    u = pc.basis().knots().front() + u * (pc.basis().knots().back() - pc.basis().knots().front());
    f << "v " << pc.eval(u) << std::endl;
  }
  f << 'l';
  for (size_t i = 1; i <= resolution; ++i)
    f << ' ' << i;
  f << std::endl;
}

void writeCurve(const PeriodicCurve &pc, std::string filename, size_t resolution = 100) {
  std::ofstream f(filename);
  for (size_t i = 0; i < resolution; ++i) {
    double u = (double)i / resolution;
    u = pc.knots.front() + u * (pc.knots.back() - pc.knots.front());
    VectorVector der;
    f << "v " << pc.derivatives(u, 0, der) << std::endl;
  }
  f << 'l';
  for (size_t i = 1; i <= resolution; ++i)
    f << ' ' << i;
  f << ' ' << 1 << std::endl;
}

void writePoly(const PointVector &pv, std::string filename) {
  std::ofstream f(filename);
  for (const auto &p : pv)
    f << "v " << p << std::endl;
  f << 'l';
  for (size_t i = 1; i <= pv.size(); ++i)
    f << ' ' << i;
  f << std::endl;
}

void showBSCurve(const BSCurve &c) {
  std::cout << "Deg: " << c.basis().degree() << std::endl;
  std::cout << "Knots:";
  for (auto k : c.basis().knots())
    std::cout << ' ' << k;
  std::cout << std::endl;
  std::cout << "Points:";
  for (const auto &p : c.controlPoints())
    std::cout << " (" << p[0] << ',' << p[1] << ')';
  std::cout << std::endl;
}

void showCurve(const PeriodicCurve &pc) {
  std::cout << "Deg: " << pc.p << std::endl;
  std::cout << "Knots:";
  for (auto k : pc.knots)
    std::cout << ' ' << k;
  std::cout << std::endl;
  std::cout << "Points:";
  for (const auto &p : pc.cp)
    std::cout << " (" << p[0] << ',' << p[1] << ')';
  std::cout << std::endl;
}

int main() {
  PeriodicCurve pc;
  pc.p = 3;
  pc.n = 3;
  pc.knots = { 1, 2, 5, 8, 15 };
  pc.cp = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } };
  double u = 1;

  auto clamped = pc.clamp(u);
  writeBSCurve(clamped, "/tmp/clamped.obj");
  writePoly(clamped.controlPoints(), "/tmp/clamped-net.obj");
  showBSCurve(clamped);

  writeCurve(pc, "/tmp/pc1.obj");
  writePoly(pc.cp, "/tmp/net1.obj");
  showCurve(pc);
  pc = pc.insertKnot(u, pc.p);
  writeCurve(pc, "/tmp/pc2.obj");
  writePoly(pc.cp, "/tmp/net2.obj");
  showCurve(pc);
}
