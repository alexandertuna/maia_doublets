#pragma once

#include "Calib.h"
#include "Objects.h"

#include <vector>

namespace doublets {

// Make mini-doublets: pairs of hits in the same (system, doublelayer, module, sensor),
// one on the inner sensor and one on the outer sensor. See maia_doublets/md.py.
//
// With apply_cuts, only MDs passing the dr and dz cuts are returned. Without,
// every pair is returned with its `ok` flag set (for validation against Python).
std::vector<MiniDoublet> makeMiniDoublets(const std::vector<Hit>& hits,
                                          const Calib& calib,
                                          bool apply_cuts);

}  // namespace doublets
