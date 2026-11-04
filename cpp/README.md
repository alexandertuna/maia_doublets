# C++ port

A C++ port of the maia_doublets track finding, for measuring how fast a C++
implementation runs. It's a Marlin processor (`DoubleLayerTrackingProc`) run by `k4run`,
with Gaudi's timing auditors.

The Python code stays the reference. The C++ code follows it step by step and is
validated against it object by object.

| Stage | C++ | Python | Validated |
|---|---|---|---|
| Hits | `DoubleLayerTrackingProc::readHits` | `slcio.py` | yes |
| MDs | `MDs.cc` | `md.py` | yes |
| T2 | `T2s.cc` | `t2.py` | yes |
| T4 | `T4s.cc` | `t4.py` | yes |
| T8 | `T8s.cc` | `t8.py` | yes |

Validation (10 µm, no time window, as in Python):
- v07, muon gun, file 300 (1000 events), without cuts: every MD, T2, T4, and T8 matches
  Python's one-to-one (same event and hits), with the same `ok` flag.
- v07, BIB, 10% file 0: MD counts per double layer match, with and without cuts.
- v07, BIB, 20% file 0: MDs, T2s, and T4s passing cuts match, per double layer.
- v06, muon gun, file 300 (1000 events), and BIB, 10% file 0, with cuts: the numbers of
  MDs, and of T2s, T4s, and T8s per double layer, passing cuts match.

## Files

- `include/Objects.h`: hits, MDs, T2, T4, T8
- `include/Constants.h`: constants from `constants.py`, including the (phi, eta) slice grids
- `include/Geometry.h`, `src/Geometry.cc`: lines to the origin, circles, the s-z fit
- `include/Pairing.h`: pairs objects in neighbouring (phi, eta) slices (the Python merges)
- `include/Calib.h`, `src/Calib.cc`: reads the cut thresholds from `maia_doublets/calibs/*.json`
- `include/MDs.h`, `src/MDs.cc`: MD making
- `include/Stages.h`, `src/T2s.cc`, `src/T4s.cc`, `src/T8s.cc`: T2, T4, T8 making
- `include/DoubleLayerTrackingProc.h`, `src/DoubleLayerTrackingProc.cc`: the Marlin processor
- `steer_doublets.py`: k4run steering
- `build.sh`, `env.sh`, `run.sh`: build, set up, run

## Usage

Inside the container (see the top-level README):

```bash
./cpp/build.sh                  # build into cpp/build, install into cpp/install
./cpp/run.sh                    # one v06 10% BIB event, MAIA-like time window
MAIA_DATASET=muonGun_pT_2p0_2p1 MAIA_FILE=300 MAIA_EVENTS=10 ./cpp/run.sh
MAIA_DETECTOR=v07 ./cpp/run.sh  # the same with v07
```

`run.sh` takes `MAIA_DETECTOR` (`v06`, the default, or `v07`), `MAIA_RESOLUTION` (e.g. `10um`), `MAIA_DATASET`,
`MAIA_FILE`, `MAIA_EVENTS`, and `MAIA_TIME_WINDOW` (`maia`, `wide`, or `none`) and writes
`output/${MAIA_DETECTOR}_${MAIA_DATASET}_digi_${MAIA_RESOLUTION}/doublets_${MAIA_TIME_WINDOW}_${MAIA_FILE}.log`.
The steering file takes `--detector` and `--resolution` and reads the cuts from
`maia_doublets/calibs/{detector}_digi_{resolution}.json` (or from `--calibFile`).
The end of the log has the CPU time and object counts of each stage.
The processor drops hits outside the time window as it reads them (`TimeMin`, `TimeMax`),
so this is part of the `hits` time. `--timeWindow none` keeps every hit, like the Python code.

For validation, `k4run cpp/steer_doublets.py --noCuts --dumpDir DIR ...` keeps
objects that fail the cuts and writes every object to `DIR/{mds,t2s,t4s,t8s}.csv`.

## Differences from Python

- Objects keep only what later stages use. Python also computes features for plots
  (charge, pt, angles, truth matching), which the C++ code leaves out.
- dr (distance of closest approach in xy) is computed in double precision. In float32
  the cross product loses ~0.1 mm for OT hits; Python's slope/intercept form agrees
  with the double result to < 1 µm.
- With cuts, MDs are made only from upper hits in the z window the dz cut allows
  (like `merge_binned` in `md.py`), found by binary search on hits sorted by z.
- Without cuts, Python's MD "all" count with `--fast-mds` counts only pairs in the
  z window, so compare it with the C++ count with cuts, not without.
