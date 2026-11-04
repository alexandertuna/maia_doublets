#include "MiniDoublets.h"

#include "Geometry.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <utility>

namespace doublets {

namespace {

// A hit's position, copied into double-module order so that each double module is contiguous
// in memory
struct SortedHit {
  float x, y, z, r;
  int sublayer;
  int index;  // into the event's hit vector
};

// A double module: the inner and outer sensors of one (system, doublelayer, module, sensor),
// like a CMS Phase-2 stacked module. sorted[begin, middle) are the hits on the lower (inner)
// sensor and sorted[middle, end) on the upper (outer) one, each sorted by z.
struct DoubleModule {
  int begin, middle, end;
};

struct DoubleModules {
  std::vector<SortedHit> sorted;
  std::vector<DoubleModule> double_modules;
};

// system:3 | doublelayer:2 | module:11 | sensor:8: the cellID's fields with the double layer in
// place of the layer, so both sensors of a double module share a key
uint32_t doubleModuleKey(const Hit& h) {
  return (uint32_t(h.system & 0x7) << 21) |
         (uint32_t(h.doublelayer() & 0x3) << 19) |
         (uint32_t(h.module & 0x7ff) << 8) |
         (uint32_t(h.sensor & 0xff) << 0);
}

// One pass of a radix sort on (double-module key << 32 | hit index): reorder `in` into `out` by
// the 12-bit digit of the double-module key at `shift`, keeping the order of equal digits
// (a stable counting sort)
constexpr int DIGIT_BITS = 12;
void radixPass(const std::vector<uint64_t>& in, int shift, std::vector<uint64_t>& out) {
  std::vector<int> start((1 << DIGIT_BITS) + 1, 0);
  for (uint64_t item : in) {
    ++start[((item >> (32 + shift)) & ((1 << DIGIT_BITS) - 1)) + 1];
  }
  for (size_t d = 1; d < start.size(); ++d) {
    start[d] += start[d - 1];
  }
  for (uint64_t item : in) {
    out[start[(item >> (32 + shift)) & ((1 << DIGIT_BITS) - 1)]++] = item;
  }
}

// Group hits by double module with a radix sort on the 24-bit double-module key (two 12-bit
// passes, so the cost scales with the number of hits, not with the 2^24 possible keys), then
// sort each double module (~1000 hits) by sublayer and z, with the hit index breaking ties.
// Each hit is one 64-bit item, double-module key << 32 | hit index, so both passes read memory
// in order.
DoubleModules groupHitsIntoDoubleModules(const std::vector<Hit>& hits) {
  std::vector<uint64_t> items(hits.size()), buffer(hits.size());
  for (size_t i = 0; i < hits.size(); ++i) {
    items[i] = (uint64_t(doubleModuleKey(hits[i])) << 32) | uint64_t(i);
  }
  radixPass(items, 0, buffer);
  radixPass(buffer, DIGIT_BITS, items);

  DoubleModules out;
  out.sorted.resize(hits.size());
  for (size_t j = 0; j < items.size(); ++j) {
    const int index = int(items[j] & 0xffffffff);
    const Hit& h = hits[index];
    out.sorted[j] = {h.x, h.y, h.z, h.r, h.sublayer(), index};
  }

  for (size_t begin = 0; begin < items.size();) {
    size_t end = begin + 1;
    while (end < items.size() && (items[end] >> 32) == (items[begin] >> 32)) {
      ++end;
    }
    auto first = out.sorted.begin() + begin;
    auto last = out.sorted.begin() + end;
    std::sort(first, last, [](const SortedHit& a, const SortedHit& b) {
      return std::tuple(a.sublayer, a.z, a.index) < std::tuple(b.sublayer, b.z, b.index);
    });
    auto on_lower_sensor = [](const SortedHit& h) { return h.sublayer == 0; };
    auto upper = std::partition_point(first, last, on_lower_sensor);
    out.double_modules.push_back({int(begin), int(upper - out.sorted.begin()), int(end)});
    begin = end;
  }
  return out;
}

// Fill everything about an MD except dr, dz, and the ok flag
MiniDoublet buildMiniDoublet(const std::vector<Hit>& hits, int lower, int upper) {
  const Hit& lo = hits[lower];
  const Hit& up = hits[upper];
  MiniDoublet md;
  md.lower = lower;
  md.upper = upper;
  md.system = lo.system;
  md.doublelayer = lo.doublelayer();
  md.gdoublelayer = lo.gdoublelayer();
  md.module = lo.module;
  md.sensor = lo.sensor;

  // position
  md.x = (lo.x + up.x) / 2;
  md.y = (lo.y + up.y) / 2;
  md.z = (lo.z + up.z) / 2;
  md.r = (lo.r + up.r) / 2;
  md.phi = std::atan2(md.y, md.x);
  md.eta = etaFromRZ(md.r, md.z);
  md.phi_slice = phiSlice(md.phi, T2_GRID);
  md.eta_slice = etaSlice(md.eta, T2_GRID);
  return md;
}

}  // namespace

std::vector<MiniDoublet> makeMiniDoublets(const std::vector<Hit>& hits,
                                          const Calib& calib,
                                          bool apply_cuts) {
  std::vector<MiniDoublet> mds;

  const DoubleModules grouped = groupHitsIntoDoubleModules(hits);
  const std::vector<SortedHit>& sorted = grouped.sorted;

  for (const DoubleModule& double_module : grouped.double_modules) {
    if (double_module.begin == double_module.middle || double_module.middle == double_module.end) {
      continue;
    }
    const Hit& first = hits[sorted[double_module.begin].index];
    const float dr_cut = calib.get("md_dr", first.system, first.doublelayer());
    const float dz_cut = calib.get("md_dz", first.system, first.doublelayer());

    auto add = [&](int lower, int upper) {
      const SortedHit& lo = sorted[lower];
      const SortedHit& up = sorted[upper];

      // the lines through both hits, extrapolated to the origin.
      // With cuts, each is applied as soon as it is known: at 100% BIB only ~1 in 400-1000
      // candidates passes.
      const float slope_rz = (up.z - lo.z) / (up.r - lo.r);
      const float dz = lo.z - lo.r * slope_rz;
      if (apply_cuts && !(std::abs(dz) < dz_cut)) {
        return;
      }
      const float dr = distanceToOrigin(lo.x, lo.y, up.x, up.y);
      if (apply_cuts && !(std::abs(dr) < dr_cut)) {
        return;
      }

      MiniDoublet md = buildMiniDoublet(hits, lo.index, up.index);
      md.dr = dr;
      md.dz = dz;
      md.ok = std::abs(md.dr) < dr_cut && std::abs(md.dz) < dz_cut;
      mds.push_back(md);
    };

    // Without cuts, keep every pair
    if (!apply_cuts) {
      for (int lower = double_module.begin; lower < double_module.middle; ++lower) {
        for (int upper = double_module.middle; upper < double_module.end; ++upper) {
          add(lower, upper);
        }
      }
      continue;
    }

    // With cuts, only look at upper hits in the z window that the dz cut allows.
    // For one pair, dz = (z_lo * r_up - z_up * r_lo) / (r_up - r_lo), so |dz| < D means
    //   z_up in [z_lo * r_up - D * (r_up - r_lo), z_lo * r_up + D * (r_up - r_lo)] / r_lo.
    // The ends are linear in r_up, so evaluating them at the smallest and largest r_up
    // in the double module gives a window containing every pair that passes.
    float r_up_min = sorted[double_module.middle].r;
    float r_up_max = r_up_min;
    for (int upper = double_module.middle; upper < double_module.end; ++upper) {
      r_up_min = std::min(r_up_min, sorted[upper].r);
      r_up_max = std::max(r_up_max, sorted[upper].r);
    }

    for (int lower = double_module.begin; lower < double_module.middle; ++lower) {
      const SortedHit& lo = sorted[lower];

      // The window assumes the upper sensor is outside the lower one. If not, try every pair.
      if (r_up_min <= lo.r) {
        for (int upper = double_module.middle; upper < double_module.end; ++upper) {
          add(lower, upper);
        }
        continue;
      }

      auto z_window = [&](float r_up) {
        const float d = dz_cut * (r_up - lo.r);
        return std::pair((lo.z * r_up - d) / lo.r, (lo.z * r_up + d) / lo.r);
      };
      const auto [lo_min, hi_min] = z_window(r_up_min);
      const auto [lo_max, hi_max] = z_window(r_up_max);
      const float z_min = std::min(lo_min, lo_max);
      const float z_max = std::max(hi_min, hi_max);

      // upper hits are sorted by z
      const SortedHit* upper_begin = &sorted[double_module.middle];
      const SortedHit* upper_end = upper_begin + (double_module.end - double_module.middle);
      auto below = [](const SortedHit& upper, float z) { return upper.z < z; };
      const SortedHit* it = std::lower_bound(upper_begin, upper_end, z_min, below);
      for (; it != upper_end && it->z <= z_max; ++it) {
        add(lower, int(it - sorted.data()));
      }
    }
  }
  return mds;
}

}  // namespace doublets
