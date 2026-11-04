#pragma once

// The objects of the track finding chain: hits -> mini-doublets (MDs) -> T2 -> T4 -> T8.
// Conventions follow the Python code (see the repo README):
//   doublelayer  = layer / 2, within a system
//   gdoublelayer = global doublelayer: IT is 0-3, OT is 4-7
//   objects are labelled by their innermost doublelayer

#include "Constants.h"

namespace EVENT {
class TrackerHit;
}

namespace doublets {

struct Hit {
  float x, y, z, r;
  int system;  // INNER_TRACKER_BARREL or OUTER_TRACKER_BARREL
  int layer;   // effective layer within the system, 0-7
  int module;  // phi module
  int sensor;  // z sensor
  EVENT::TrackerHit* lcio;

  int doublelayer() const { return layer / 2; }
  int sublayer() const { return layer % 2; }  // 0 = inner sensor of the double layer
  int gdoublelayer() const {
    return doublelayer() + (system == OUTER_TRACKER_BARREL ? N_DOUBLELAYERS_PER_SYSTEM : 0);
  }
};

// Each object keeps only what later stages need. Python also computes features for
// plots (charge, pt, angles), which the C++ code leaves out.
//
// Every object has:
//   dr: xy distance of closest approach to the origin of the line through its two halves
//   dz: rz z at r = 0 of the line through its two halves
//   x, y, z, r, phi, eta: position, the average of its two halves
//   phi_slice, eta_slice: (phi, eta) grid used to pair it into the next stage
//   ok: passed this stage's cuts

struct MiniDoublet {
  int lower, upper;  // indices into the event's hit vector
  int system, doublelayer, gdoublelayer, module, sensor;
  float dr, dz;
  float x, y, z, r, phi, eta;
  int phi_slice, eta_slice;
  bool ok;
};

// Two MDs in neighbouring double layers (4 layers)
struct T2 {
  int lower, upper;  // indices into the event's MD vector
  int hits[N_LAYERS_IN_T2];  // indices into the event's hit vector, inside out
  int system, doublelayer, gdoublelayer;
  float dr, dz;
  float x, y, z, r, phi, eta;
  int phi_slice, eta_slice;
  float dphi, chi2_xy, chi2_sz;
  bool ok;
};

// Two T2s (8 layers)
struct T4 {
  int lower, upper;  // indices into the event's T2 vector
  int hits[N_LAYERS_IN_T4];
  int gdoublelayer;
  float dr, dz;
  float x, y, z, r, phi, eta;
  int phi_slice, eta_slice;
  float dphi, chi2_xy, chi2_sz;
  bool ok;
};

// Two T4s (16 layers): one in the IT and one in the OT
struct T8 {
  int lower, upper;  // indices into the event's T4 vector
  int hits[N_LAYERS_IN_T8];
  int gdoublelayer;
  float dr, dz;
  float chi2_xy, chi2_sz;
  bool ok;
};

}  // namespace doublets
