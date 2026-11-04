#pragma once

// Constants shared by the stages, as in maia_doublets/constants.py

#include <cmath>

namespace doublets {

constexpr double PI = M_PI;
constexpr double SPEED_OF_LIGHT = 299.792458;  // mm/ns
constexpr double MAGNETIC_FIELD = 5.0;         // T
constexpr double BAD_CHI2 = 1e6;               // when the circle fit fails

// The barrel tracker systems (the cellID's system field) and their layers
constexpr int INNER_TRACKER_BARREL = 3;
constexpr int OUTER_TRACKER_BARREL = 5;
constexpr int N_LAYERS_PER_SYSTEM = 8;  // effective layers, i.e. 4 double layers
constexpr int N_DOUBLELAYERS_PER_SYSTEM = N_LAYERS_PER_SYSTEM / 2;
constexpr int N_GDOUBLELAYERS = 2 * N_DOUBLELAYERS_PER_SYSTEM;  // IT and OT: IT is 0-3, OT is 4-7

// Hits per object
constexpr int N_LAYERS_IN_T2 = 4;
constexpr int N_LAYERS_IN_T4 = 8;
constexpr int N_LAYERS_IN_T8 = 16;

// (phi, eta) grids for pairing objects into the next stage.
// MDs are binned on the T2 grid, T2s on the T4 grid, and T4s on the T8 grid.
// The slice widths are multiplied by 2 to speed up the processing (fewer slices).
constexpr double DETECTOR_MAX_PHI = PI;
constexpr double DETECTOR_MAX_ETA = 2.5;  // must include all background hits
constexpr double MAX_T2_DPHI = 0.10 * 2;
constexpr double MAX_T2_DETA = 0.01 * 2;
constexpr double MAX_T4_DPHI = 0.25 * 2;
constexpr double MAX_T4_DETA = 0.025 * 2;
constexpr double MAX_T8_DPHI = 0.25 * 2;
constexpr double MAX_T8_DETA = 0.025 * 2;

struct SliceGrid {
  int n_phi;
  int n_eta;
};
constexpr SliceGrid T2_GRID{int(2 * DETECTOR_MAX_PHI / MAX_T2_DPHI), int(2 * DETECTOR_MAX_ETA / MAX_T2_DETA)};
constexpr SliceGrid T4_GRID{int(2 * DETECTOR_MAX_PHI / MAX_T4_DPHI), int(2 * DETECTOR_MAX_ETA / MAX_T4_DETA)};
constexpr SliceGrid T8_GRID{int(2 * DETECTOR_MAX_PHI / MAX_T8_DPHI), int(2 * DETECTOR_MAX_ETA / MAX_T8_DETA)};

}  // namespace doublets
