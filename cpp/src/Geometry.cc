#include "Geometry.h"

namespace doublets {

Circle circleThrough(const Hit& a, const Hit& b, const Hit& c) {
  const double xa = a.x, ya = a.y, ra2 = double(a.r) * a.r;
  const double xb = b.x, yb = b.y, rb2 = double(b.r) * b.r;
  const double xc = c.x, yc = c.y, rc2 = double(c.r) * c.r;

  Circle circle;
  const double d = 2 * (xa * (yb - yc) + xb * (yc - ya) + xc * (ya - yb));
  circle.ok = d != 0;
  circle.x = (ra2 * (yb - yc) + rb2 * (yc - ya) + rc2 * (ya - yb)) / d;
  circle.y = (ra2 * (xc - xb) + rb2 * (xa - xc) + rc2 * (xb - xa)) / d;
  circle.radius = std::hypot(xa - circle.x, ya - circle.y);
  return circle;
}

double chi2XY(const std::vector<Hit>& hits, const int* indices, int n, const Circle& circle) {
  if (!circle.ok) {
    return n * BAD_CHI2;
  }
  double chi2 = 0;
  for (int i = 0; i < n; ++i) {
    const Hit& h = hits[indices[i]];
    const double diff = std::hypot(h.x - circle.x, h.y - circle.y) - circle.radius;
    chi2 += diff * diff;
  }
  return chi2;
}

double chi2SZ(const std::vector<Hit>& hits, const int* indices, int n, const Circle& circle) {
  if (!circle.ok) {
    return BAD_CHI2;
  }

  // arc length of each hit from the first one
  if (n > N_LAYERS_IN_T8) {  // the most hits of any object
    return BAD_CHI2;
  }
  double s[N_LAYERS_IN_T8], z[N_LAYERS_IN_T8];
  const Hit& first = hits[indices[0]];
  const double phi_first = std::atan2(first.y - circle.y, first.x - circle.x);
  for (int i = 0; i < n; ++i) {
    const Hit& h = hits[indices[i]];
    const double phi = std::atan2(h.y - circle.y, h.x - circle.x);
    s[i] = circle.radius * wrapPhi(phi - phi_first);
    z[i] = h.z;
  }

  // least-squares line through (s, z)
  double s_mean = 0, z_mean = 0;
  for (int i = 0; i < n; ++i) {
    s_mean += s[i] / n;
    z_mean += z[i] / n;
  }
  double s_ss = 0, s_sz = 0;
  for (int i = 0; i < n; ++i) {
    s_ss += (s[i] - s_mean) * (s[i] - s_mean);
    s_sz += (s[i] - s_mean) * (z[i] - z_mean);
  }
  const double tanlambda = s_sz / s_ss;
  const double z0 = z_mean - tanlambda * s_mean;

  double chi2 = 0;
  for (int i = 0; i < n; ++i) {
    const double resid = z[i] - (z0 + tanlambda * s[i]);
    chi2 += resid * resid;
  }
  return chi2;
}

}  // namespace doublets
