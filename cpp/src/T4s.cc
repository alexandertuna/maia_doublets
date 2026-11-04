#include "Geometry.h"
#include "Pairing.h"
#include "Stages.h"

namespace doublets {

std::vector<T4> makeT4s(const std::vector<Hit>& hits,
                        const std::vector<T2>& t2s,
                        const Calib& calib,
                        bool apply_cuts) {
  std::vector<T4> t4s;

  for (auto [lower_gdl, upper_gdl] : {std::pair(0, 2), std::pair(4, 6)}) {
    const std::vector<int> lowers = okInGDoubleLayer(t2s, lower_gdl);
    const std::vector<int> uppers = okInGDoubleLayer(t2s, upper_gdl);
    if (lowers.empty() || uppers.empty()) {
      continue;
    }
    const float dz_cut = calib.get("t4_dz", lower_gdl);
    const float dr_cut = calib.get("t4_dr", lower_gdl);
    const float chi2_xy_cut = calib.get("t4_chi2_xy", lower_gdl);
    const float chi2_sz_cut = calib.get("t4_chi2_sz", lower_gdl);

    forEachNeighbouringPair(t2s, lowers, uppers, T4_GRID.n_phi, [&](int lower, int upper) {
      const T2& lo = t2s[lower];
      const T2& up = t2s[upper];

      T4 t4;
      t4.lower = lower;
      t4.upper = upper;
      for (int i = 0; i < N_LAYERS_IN_T2; ++i) {
        t4.hits[i] = lo.hits[i];
        t4.hits[i + N_LAYERS_IN_T2] = up.hits[i];
      }
      t4.gdoublelayer = lo.gdoublelayer;

      // With cuts, each is applied as soon as it is known, cheapest first: most pairs fail.
      t4.dphi = wrapPhi(up.phi - lo.phi);
      if (apply_cuts && !(std::abs(t4.dphi) < PI / 2)) {
        return;
      }

      // the lines through both T2s, extrapolated to the origin
      t4.dz = zAtOrigin(lo.r, lo.z, up.r, up.z);
      if (apply_cuts && !(std::abs(t4.dz) < dz_cut)) {
        return;
      }
      t4.dr = distanceToOrigin(lo.x, lo.y, up.x, up.y);
      if (apply_cuts && !(std::abs(t4.dr) < dr_cut)) {
        return;
      }

      // circle through the first, middle, and last hits, and how far all 8 are from it
      const Circle circle = circleThrough(hits[t4.hits[0]], hits[t4.hits[4]], hits[t4.hits[7]]);
      t4.chi2_xy = chi2XY(hits, t4.hits, N_LAYERS_IN_T4, circle);
      if (apply_cuts && !(t4.chi2_xy < chi2_xy_cut)) {
        return;
      }
      t4.chi2_sz = chi2SZ(hits, t4.hits, N_LAYERS_IN_T4, circle);
      if (apply_cuts && !(t4.chi2_sz < chi2_sz_cut)) {
        return;
      }

      // position (eta is the average of the T2s', as in Python)
      t4.x = (lo.x + up.x) / 2;
      t4.y = (lo.y + up.y) / 2;
      t4.z = (lo.z + up.z) / 2;
      t4.r = (lo.r + up.r) / 2;
      t4.phi = std::atan2(t4.y, t4.x);
      t4.eta = (lo.eta + up.eta) / 2;
      t4.phi_slice = phiSlice(t4.phi, T8_GRID);
      t4.eta_slice = etaSlice(t4.eta, T8_GRID);

      t4.ok = std::abs(t4.dphi) < PI / 2 &&
              std::abs(t4.dz) < dz_cut &&
              std::abs(t4.dr) < dr_cut &&
              t4.chi2_xy < chi2_xy_cut &&
              t4.chi2_sz < chi2_sz_cut;
      if (t4.ok || !apply_cuts) {
        t4s.push_back(t4);
      }
    });
  }
  return t4s;
}

}  // namespace doublets
