#include "Geometry.h"
#include "Pairing.h"
#include "Stages.h"

namespace doublets {

std::vector<T2> makeT2s(const std::vector<Hit>& hits,
                        const std::vector<MD>& mds,
                        const Calib& calib,
                        bool apply_cuts) {
  std::vector<T2> t2s;

  for (auto [lower_gdl, upper_gdl] : {std::pair(0, 1),
                                      std::pair(2, 3),
                                      std::pair(4, 5),
                                      std::pair(6, 7)}) {
    const std::vector<int> lowers = okInGDoubleLayer(mds, lower_gdl);
    const std::vector<int> uppers = okInGDoubleLayer(mds, upper_gdl);
    if (lowers.empty() || uppers.empty()) {
      continue;
    }
    const MD& first = mds[lowers.front()];
    const float dz_cut = calib.get("t2_dz", first.system, first.doublelayer);
    const float dr_cut = calib.get("t2_dr", first.system, first.doublelayer);
    const float chi2_xy_cut = calib.get("t2_chi2_xy", first.system, first.doublelayer);
    const float chi2_sz_cut = calib.get("t2_chi2_sz", first.system, first.doublelayer);

    forEachNeighbouringPair(mds, lowers, uppers, T2_GRID.n_phi, [&](int lower, int upper) {
      const MD& lo = mds[lower];
      const MD& up = mds[upper];

      T2 t2;
      t2.lower = lower;
      t2.upper = upper;
      t2.hits[0] = lo.lower;
      t2.hits[1] = lo.upper;
      t2.hits[2] = up.lower;
      t2.hits[3] = up.upper;
      t2.system = lo.system;
      t2.doublelayer = lo.doublelayer;
      t2.gdoublelayer = lo.gdoublelayer;

      // With cuts, each is applied as soon as it is known, cheapest first: most pairs fail.

      // the lines through both MDs, extrapolated to the origin
      t2.dz = zAtOrigin(lo.r, lo.z, up.r, up.z);
      if (apply_cuts && !(std::abs(t2.dz) < dz_cut)) {
        return;
      }
      t2.dr = distanceToOrigin(lo.x, lo.y, up.x, up.y);
      if (apply_cuts && !(std::abs(t2.dr) < dr_cut)) {
        return;
      }
      t2.dphi = wrapPhi(up.phi - lo.phi);
      if (apply_cuts && !(std::abs(t2.dphi) < PI / 2)) {
        return;
      }

      // circle through the first three hits, and how far the fourth is from it
      const Circle circle = circleThrough(hits[t2.hits[0]], hits[t2.hits[1]], hits[t2.hits[2]]);
      t2.chi2_xy = chi2XY(hits, &t2.hits[3], 1, circle);
      if (apply_cuts && !(t2.chi2_xy < chi2_xy_cut)) {
        return;
      }
      t2.chi2_sz = chi2SZ(hits, t2.hits, N_LAYERS_IN_T2, circle);
      if (apply_cuts && !(t2.chi2_sz < chi2_sz_cut)) {
        return;
      }

      // position
      t2.x = (lo.x + up.x) / 2;
      t2.y = (lo.y + up.y) / 2;
      t2.z = (lo.z + up.z) / 2;
      t2.r = (lo.r + up.r) / 2;
      t2.phi = std::atan2(t2.y, t2.x);
      t2.eta = etaFromRZ(t2.r, t2.z);
      t2.phi_slice = phiSlice(t2.phi, T4_GRID);
      t2.eta_slice = etaSlice(t2.eta, T4_GRID);

      t2.ok = std::abs(t2.dz) < dz_cut &&
              std::abs(t2.dr) < dr_cut &&
              std::abs(t2.dphi) < PI / 2 &&
              t2.chi2_xy < chi2_xy_cut &&
              t2.chi2_sz < chi2_sz_cut;
      if (t2.ok || !apply_cuts) {
        t2s.push_back(t2);
      }
    });
  }
  return t2s;
}

}  // namespace doublets
