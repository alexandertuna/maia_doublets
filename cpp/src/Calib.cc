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
