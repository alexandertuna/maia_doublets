#pragma once

#include <nlohmann/json.hpp>

#include <string>

namespace doublets {

// Cut thresholds from maia_doublets/calibs/{geo}_{sim|digi_XXum}.json.
// md_* and t2_* are keyed by (system, doublelayer); t4_* and t8_* by gdoublelayer.
// A missing entry gives 0, so nothing passes that cut (as in calib.py).
class Calib {
 public:
  explicit Calib(const std::string& path);

  float get(const std::string& feature, int system, int doublelayer) const;
  float get(const std::string& feature, int gdoublelayer) const;

 private:
  nlohmann::json _json;
};

}  // namespace doublets
