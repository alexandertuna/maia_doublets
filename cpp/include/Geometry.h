#pragma once

// Geometry shared by the stages: lines, circles, and the s-z line fit

#include "Constants.h"
#include "Objects.h"

#include <cmath>
#include <vector>

namespace doublets {

// Into [-pi, pi)
inline double wrapPhi(double phi) {
  return phi - 2 * PI * std::floor((phi + PI) / (2 * PI));
}

inline float etaFromRZ(float r, float z) {
  return -std::log(std::tan(std::atan2(r, z) / 2));
}

// Slices are computed in float, like Python does with its float32 columns,
// so that objects near a slice edge land in the same slice
inline int phiSlice(float phi, const SliceGrid& grid) {
  return std::floor((phi + float(DETECTOR_MAX_PHI)) / float(2 * DETECTOR_MAX_PHI) * grid.n_phi);
}

inline int etaSlice(float eta, const SliceGrid& grid) {
  return std::floor((eta + float(DETECTOR_MAX_ETA)) / float(2 * DETECTOR_MAX_ETA) * grid.n_eta);
}

// xy: distance from the origin to the line through two points = |cross product| / length.
// In double: the cross product cancels two terms of ~1e6 mm^2, too much for a float.
inline double distanceToOrigin(double x0, double y0, double x1, double y1) {
  return std::abs(x0 * y1 - x1 * y0) / std::hypot(x1 - x0, y1 - y0);
}

// rz: z of the line through two points, at r = 0
inline double zAtOrigin(double r0, double z0, double r1, double z1) {
  return z0 - r0 * (z1 - z0) / (r1 - r0);
}

struct Circle {
  double x, y, radius;
  bool ok;  // false if the three points are on a line
};

// The circle through three hits. Like the Python code, it uses each hit's stored r for r^2.
Circle circleThrough(const Hit& a, const Hit& b, const Hit& c);

// Sum of squared distances of the given hits from the circle
double chi2XY(const std::vector<Hit>& hits, const int* indices, int n, const Circle& circle);

// Fit z = z0 + tanlambda * s, with s the arc length along the circle from the first hit,
// and return the sum of squared z residuals
double chi2SZ(const std::vector<Hit>& hits, const int* indices, int n, const Circle& circle);

}  // namespace doublets
