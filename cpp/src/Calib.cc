#include "Calib.h"

#include <fstream>
#include <stdexcept>

namespace doublets {

Calib::Calib(const std::string& path) {
  std::ifstream file(path);
  if (!file) {
    throw std::runtime_error("Cannot open calibration file " + path);
  }
  _json = nlohmann::json::parse(file);
}

float Calib::get(const std::string& feature, int system, int doublelayer) const {
  if (!_json.contains(feature)) {
    return 0;
  }
  const auto& by_system = _json[feature];
  const std::string sy = std::to_string(system);
  const std::string dl = std::to_string(doublelayer);
  if (!by_system.contains(sy) || !by_system[sy].contains(dl)) {
    return 0;
  }
  return by_system[sy][dl].get<float>();
}

float Calib::get(const std::string& feature, int gdoublelayer) const {
  if (!_json.contains(feature)) {
    return 0;
  }
  const std::string gdl = std::to_string(gdoublelayer);
  if (!_json[feature].contains(gdl)) {
    return 0;
  }
  return _json[feature][gdl].get<float>();
}

}  // namespace doublets
