#pragma once

#include "Calib.h"
#include "Objects.h"

#include <vector>

namespace doublets {

// Each stage pairs objects of the previous stage that passed its cuts, in neighbouring
// double layers and neighbouring (phi, eta) slices.
//
// With apply_cuts, only objects passing this stage's cuts are returned. Without, every
// candidate is returned with its `ok` flag set (for validation against Python).

// T2: MDs in gdoublelayers (0, 1), (2, 3), (4, 5), (6, 7). See maia_doublets/t2.py.
std::vector<T2> makeT2s(const std::vector<Hit>& hits,
                        const std::vector<MD>& mds,
                        const Calib& calib,
                        bool apply_cuts);

// T4: T2s in gdoublelayers (0, 2) and (4, 6). See maia_doublets/t4.py.
std::vector<T4> makeT4s(const std::vector<Hit>& hits,
                        const std::vector<T2>& t2s,
                        const Calib& calib,
                        bool apply_cuts);

// T8: T4s in gdoublelayers (0, 4), i.e. IT and OT. See maia_doublets/t8.py.
std::vector<T8> makeT8s(const std::vector<Hit>& hits,
                        const std::vector<T4>& t4s,
                        const Calib& calib,
                        bool apply_cuts);

}  // namespace doublets
