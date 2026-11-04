#include "Geometry.h"
#include "Pairing.h"
#include "Stages.h"

namespace doublets {

std::vector<T8> makeT8s(const std::vector<Hit>& hits,
                        const std::vector<T4>& t4s,
                        const Calib& calib,
                        bool apply_cuts) {
  std::vector<T8> t8s;

  const int lower_gdl = 0;  // IT
  const int upper_gdl = 4;  // OT
  const std::vector<int> lowers = okInGDoubleLayer(t4s, lower_gdl);
  const std::vector<int> uppers = okInGDoubleLayer(t4s, upper_gdl);
  const float dz_cut = calib.get("t8_dz", lower_gdl);
  const float dr_cut = calib.get("t8_dr", lower_gdl);
  const float chi2_xy_cut = calib.get("t8_chi2_xy", lower_gdl);
  const float chi2_sz_cut = calib.get("t8_chi2_sz", lower_gdl);

  forEachNeighbouringPair(t4s, lowers, uppers, T8_GRID.n_phi, [&](int lower, int upper) {
    const T4& lo = t4s[lower];
    const T4& up = t4s[upper];
    T8 t8;
    t8.lower = lower;
    t8.upper = upper;
    for (int i = 0; i < N_LAYERS_IN_T4; ++i) {
      t8.hits[i] = lo.hits[i];
      t8.hits[i + N_LAYERS_IN_T4] = up.hits[i];
    }
    t8.gdoublelayer = lo.gdoublelayer;

    // the lines through both T4s, extrapolated to the origin
    t8.dr = distanceToOrigin(lo.x, lo.y, up.x, up.y);
    t8.dz = zAtOrigin(lo.r, lo.z, up.r, up.z);

    // circle through the first, middle, and last hits, and how far the other 13 are from it
    const int first = 0, middle = N_LAYERS_IN_T4, last = N_LAYERS_IN_T8 - 1;
    const Circle circle = circleThrough(hits[t8.hits[first]],
                                        hits[t8.hits[middle]],
                                        hits[t8.hits[last]]);
    constexpr int N_OTHERS = N_LAYERS_IN_T8 - 3;
    int others[N_OTHERS];
    int n_others = 0;
    for (int i = 0; i < N_LAYERS_IN_T8; ++i) {
      if (i != first && i != middle && i != last) {
        others[n_others++] = t8.hits[i];
      }
    }
    t8.chi2_xy = chi2XY(hits, others, N_OTHERS, circle);
    t8.chi2_sz = chi2SZ(hits, t8.hits, N_LAYERS_IN_T8, circle);

    t8.ok = std::abs(t8.dz) < dz_cut &&
            std::abs(t8.dr) < dr_cut &&
            t8.chi2_xy < chi2_xy_cut &&
            t8.chi2_sz < chi2_sz_cut;
    if (t8.ok || !apply_cuts) {
      t8s.push_back(t8);
    }
  });
  return t8s;
}

}  // namespace doublets
