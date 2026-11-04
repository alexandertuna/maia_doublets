#pragma once

#include <map>
#include <utility>
#include <vector>

namespace doublets {

// Call pair(lower, upper) for every lower and upper object in the same or a neighbouring
// (phi, eta) slice. This is the Python merge on (phi_slice, eta_slice) with the upper
// slices shifted by -1, 0, +1: phi wraps around, eta doesn't.
//
// `lowers` and `uppers` are indices into `objects`, which have phi_slice and eta_slice.
template <typename Object, typename Function>
void forEachNeighbouringPair(const std::vector<Object>& objects, const std::vector<int>& lowers,
                             const std::vector<int>& uppers, int n_phi_slices, Function pair) {
  // File each upper object under the 9 slices it can pair with
  std::map<std::pair<int, int>, std::vector<int>> uppers_by_slice;
  for (int upper : uppers) {
    for (int phi_shift : {-1, 0, 1}) {
      for (int eta_shift : {-1, 0, 1}) {
        const int phi = ((objects[upper].phi_slice + phi_shift) % n_phi_slices + n_phi_slices) % n_phi_slices;
        const int eta = objects[upper].eta_slice + eta_shift;
        uppers_by_slice[{phi, eta}].push_back(upper);
      }
    }
  }

  for (int lower : lowers) {
    auto found = uppers_by_slice.find({objects[lower].phi_slice, objects[lower].eta_slice});
    if (found == uppers_by_slice.end()) {
      continue;
    }
    for (int upper : found->second) {
      pair(lower, upper);
    }
  }
}

// Indices of the objects that passed the previous stage's cuts, in one gdoublelayer
template <typename Object>
std::vector<int> okInGDoubleLayer(const std::vector<Object>& objects, int gdoublelayer) {
  std::vector<int> indices;
  for (int i = 0; i < int(objects.size()); ++i) {
    if (objects[i].ok && objects[i].gdoublelayer == gdoublelayer) {
      indices.push_back(i);
    }
  }
  return indices;
}

}  // namespace doublets
