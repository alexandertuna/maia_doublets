#pragma once

#include "Calib.h"
#include "Objects.h"

#include <marlin/Processor.h>

#include <array>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// Marlin processor running the maia_doublets track finding on IT/OT barrel hits:
// hits -> mini-doublets (MDs) -> T2 -> T4 -> T8.
// Hits outside [TimeMin, TimeMax] (the hit time as stored, which the digitiser has
// already corrected for time of flight) are dropped when reading them.
//
// The CPU time of each stage is printed in end(), next to the object counts.
class DoubleLayerTrackingProc : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() override { return new DoubleLayerTrackingProc; }

  DoubleLayerTrackingProc();
  DoubleLayerTrackingProc(const DoubleLayerTrackingProc&) = delete;
  DoubleLayerTrackingProc& operator=(const DoubleLayerTrackingProc&) = delete;

  void init() override;
  void processEvent(LCEvent* evt) override;
  void end() override;

 private:
  std::vector<doublets::Hit> readHits(LCEvent* evt) const;

  // CPU time and object counts of one stage, summed over the job
  struct StageSummary {
    double cpu = 0;
    std::array<long, doublets::N_GDOUBLELAYERS> all{};  // by gdoublelayer
    std::array<long, doublets::N_GDOUBLELAYERS> ok{};   // ... that pass the cuts

    template <typename Object>
    void count(const std::vector<Object>& objects);
  };
  void printSummary(const std::string& name, const StageSummary& summary) const;

  // parameters
  std::vector<std::string> _hitCollections;
  std::string _calibFile;
  float _timeMin = -std::numeric_limits<float>::infinity();  // ns
  float _timeMax = std::numeric_limits<float>::infinity();   // ns
  bool _applyCuts = true;
  std::string _dumpDir;

  std::unique_ptr<doublets::Calib> _calib;

  // validation dumps, one CSV per stage
  std::ofstream _dumpMDs, _dumpT2s, _dumpT4s, _dumpT8s;

  int _nEvents = 0;
  long _nHits = 0;
  double _cpuReadHits = 0;
  StageSummary _mds, _t2s, _t4s, _t8s;
};
