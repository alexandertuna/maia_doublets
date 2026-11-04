#include "DoubleLayerTrackingProc.h"

#include "MiniDoublets.h"
#include "Stages.h"

#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHit.h>

#include <cmath>
#include <ctime>
#include <iomanip>

using namespace doublets;

DoubleLayerTrackingProc aDoubleLayerTrackingProc;

namespace {

// CPU time of this process, in seconds
double cpuSeconds() {
  return double(std::clock()) / CLOCKS_PER_SEC;
}

// The v06/v07 IT/OT barrel cellID, as in slcio.py decode_cellid:
// system:5, side:-2, sensor:8, layer:5, module:11 (the "effective" layout)
int decodeSystem(unsigned int cellid) { return cellid & 0x1f; }
int decodeSensor(unsigned int cellid) { return (cellid >> 7) & 0xff; }
int decodeLayer(unsigned int cellid) { return (cellid >> 15) & 0x1f; }
int decodeModule(unsigned int cellid) { return (cellid >> 20) & 0x7ff; }

// Validation dumps: one row per object, with the positions of its hits (to match objects
// with Python's), a few features, and the ok flag
void writeHeader(std::ofstream& out, int n_hits, const std::string& features) {
  out << "i_event,gdoublelayer";
  for (int i = 0; i < n_hits; ++i) {
    out << ",x_" << i << ",y_" << i << ",z_" << i;
  }
  out << "," << features << ",ok\n";
  out << std::setprecision(9);  // enough digits to round-trip a float
}

void writeHits(std::ofstream& out, int i_event, int gdoublelayer, const std::vector<Hit>& hits, const int* indices,
               int n_hits) {
  out << i_event << "," << gdoublelayer;
  for (int i = 0; i < n_hits; ++i) {
    const Hit& h = hits[indices[i]];
    out << "," << h.x << "," << h.y << "," << h.z;
  }
}

}  // namespace

template <typename Object>
void DoubleLayerTrackingProc::StageSummary::count(const std::vector<Object>& objects) {
  for (const Object& object : objects) {
    all[object.gdoublelayer] += 1;
    ok[object.gdoublelayer] += object.ok;
  }
}

DoubleLayerTrackingProc::DoubleLayerTrackingProc() : marlin::Processor("DoubleLayerTrackingProc") {
  _description = "maia_doublets track finding on IT/OT barrel double layers";

  registerInputCollections(LCIO::TRACKERHITPLANE, "TrackerHitCollectionNames", "IT/OT barrel hit collections",
                           _hitCollections, std::vector<std::string>{"IBTrackerHits", "OBTrackerHits"});
  registerProcessorParameter("CalibFile", "Cut thresholds, maia_doublets/calibs/*.json", _calibFile,
                             std::string(""));
  registerProcessorParameter("TimeMin", "Drop hits with time below this [ns]", _timeMin, _timeMin);
  registerProcessorParameter("TimeMax", "Drop hits with time above this [ns]", _timeMax, _timeMax);
  registerProcessorParameter("ApplyCuts", "Keep only objects that pass the cuts", _applyCuts, true);
  registerProcessorParameter("DumpDir", "If set, write every object to {mds,t2s,t4s,t8s}.csv here (for validation)",
                             _dumpDir, std::string(""));
}

void DoubleLayerTrackingProc::init() {
  printParameters();
  _calib = std::make_unique<Calib>(_calibFile);
  if (!_dumpDir.empty()) {
    _dumpMDs.open(_dumpDir + "/mds.csv");
    _dumpT2s.open(_dumpDir + "/t2s.csv");
    _dumpT4s.open(_dumpDir + "/t4s.csv");
    _dumpT8s.open(_dumpDir + "/t8s.csv");
    writeHeader(_dumpMDs, N_LAYERS_IN_MD, "dr,dz");
    writeHeader(_dumpT2s, N_LAYERS_IN_T2, "dr,dz,chi2_xy,chi2_sz");
    writeHeader(_dumpT4s, N_LAYERS_IN_T4, "dr,dz,chi2_xy,chi2_sz");
    writeHeader(_dumpT8s, N_LAYERS_IN_T8, "dr,dz,chi2_xy,chi2_sz");
  }
}

std::vector<Hit> DoubleLayerTrackingProc::readHits(LCEvent* evt) const {
  std::vector<Hit> hits;
  size_t n_hits = 0;
  for (const std::string& name : _hitCollections) {
    n_hits += evt->getCollection(name)->getNumberOfElements();
  }
  hits.reserve(n_hits);
  for (const std::string& name : _hitCollections) {
    LCCollection* col = evt->getCollection(name);
    for (int i = 0; i < col->getNumberOfElements(); ++i) {
      auto* lcio = static_cast<EVENT::TrackerHit*>(col->getElementAt(i));
      if (lcio->getTime() < _timeMin || lcio->getTime() > _timeMax) {
        continue;
      }
      const unsigned int cellid = lcio->getCellID0();
      const double* pos = lcio->getPosition();

      Hit hit;
      hit.system = decodeSystem(cellid);
      hit.layer = decodeLayer(cellid);
      hit.module = decodeModule(cellid);
      hit.sensor = decodeSensor(cellid);
      if (hit.system != INNER_TRACKER_BARREL && hit.system != OUTER_TRACKER_BARREL) {
        continue;
      }
      if (hit.layer >= N_LAYERS_PER_SYSTEM) {
        continue;
      }
      // Like Python: r from the double-precision position, then everything as float
      hit.x = pos[0];
      hit.y = pos[1];
      hit.z = pos[2];
      hit.r = std::hypot(pos[0], pos[1]);
      hit.lcio = lcio;
      hits.push_back(hit);
    }
  }
  return hits;
}

void DoubleLayerTrackingProc::processEvent(LCEvent* evt) {
  double start = cpuSeconds();
  const std::vector<Hit> hits = readHits(evt);
  _cpuReadHits += cpuSeconds() - start;

  start = cpuSeconds();
  const std::vector<MiniDoublet> mds = makeMiniDoublets(hits, *_calib, _applyCuts);
  _mds.cpu += cpuSeconds() - start;

  start = cpuSeconds();
  const std::vector<T2> t2s = makeT2s(hits, mds, *_calib, _applyCuts);
  _t2s.cpu += cpuSeconds() - start;

  start = cpuSeconds();
  const std::vector<T4> t4s = makeT4s(hits, t2s, *_calib, _applyCuts);
  _t4s.cpu += cpuSeconds() - start;

  start = cpuSeconds();
  const std::vector<T8> t8s = makeT8s(hits, t4s, *_calib, _applyCuts);
  _t8s.cpu += cpuSeconds() - start;

  streamlog_out(MESSAGE) << "event " << _nEvents << ": " << hits.size() << " hits, " << mds.size() << " MDs, "
                         << t2s.size() << " T2s, " << t4s.size() << " T4s, " << t8s.size() << " T8s" << std::endl;

  if (!_dumpDir.empty()) {
    for (const MiniDoublet& md : mds) {
      const int indices[N_LAYERS_IN_MD] = {md.lower, md.upper};
      writeHits(_dumpMDs, _nEvents, md.gdoublelayer, hits, indices, N_LAYERS_IN_MD);
      _dumpMDs << "," << md.dr << "," << md.dz << "," << md.ok << "\n";
    }
    for (const T2& t2 : t2s) {
      writeHits(_dumpT2s, _nEvents, t2.gdoublelayer, hits, t2.hits, N_LAYERS_IN_T2);
      _dumpT2s << "," << t2.dr << "," << t2.dz << "," << t2.chi2_xy << "," << t2.chi2_sz << "," << t2.ok << "\n";
    }
    for (const T4& t4 : t4s) {
      writeHits(_dumpT4s, _nEvents, t4.gdoublelayer, hits, t4.hits, N_LAYERS_IN_T4);
      _dumpT4s << "," << t4.dr << "," << t4.dz << "," << t4.chi2_xy << "," << t4.chi2_sz << "," << t4.ok << "\n";
    }
    for (const T8& t8 : t8s) {
      writeHits(_dumpT8s, _nEvents, t8.gdoublelayer, hits, t8.hits, N_LAYERS_IN_T8);
      _dumpT8s << "," << t8.dr << "," << t8.dz << "," << t8.chi2_xy << "," << t8.chi2_sz << "," << t8.ok << "\n";
    }
  }

  _nHits += hits.size();
  _mds.count(mds);
  _t2s.count(t2s);
  _t4s.count(t4s);
  _t8s.count(t8s);
  ++_nEvents;
}

void DoubleLayerTrackingProc::printSummary(const std::string& name, const StageSummary& summary) const {
  long all = 0, ok = 0;
  for (int gdl = 0; gdl < N_GDOUBLELAYERS; ++gdl) {
    all += summary.all[gdl];
    ok += summary.ok[gdl];
  }
  streamlog_out(MESSAGE) << "  " << name << ": " << summary.cpu << " s CPU, " << all << " made, " << ok
                         << " passing cuts" << std::endl;
  for (int gdl = 0; gdl < N_GDOUBLELAYERS; ++gdl) {
    if (summary.all[gdl] > 0) {
      streamlog_out(MESSAGE) << "    gdoublelayer " << gdl << ": " << summary.all[gdl] << " made, " << summary.ok[gdl]
                             << " passing cuts" << std::endl;
    }
  }
}

void DoubleLayerTrackingProc::end() {
  for (std::ofstream* dump : {&_dumpMDs, &_dumpT2s, &_dumpT4s, &_dumpT8s}) {
    if (dump->is_open()) {
      dump->close();
    }
  }
  streamlog_out(MESSAGE) << "Summary over " << _nEvents << " events:" << std::endl
                         << "  hits: " << _cpuReadHits << " s CPU, " << _nHits << " read" << std::endl;
  printSummary("MDs", _mds);
  printSummary("T2s", _t2s);
  printSummary("T4s", _t4s);
  printSummary("T8s", _t8s);
}
