# CLAUDE.md

## Working agreement (read this first)

Alex owns this repo and this physics result. These rules take priority over speed:

1. **Ask before editing anything in the repo.** Describe the change (file, what, why) and wait for a yes. This includes code, `run/run.sh`, calibration JSONs, paper scripts, and docs. One approval covers one change, not later ones.
2. **The scratchpad is fair game.** Exploratory scripts, quick checks, and throwaway plots go in the session scratchpad directory, never in the repo. Give the full absolute path of anything written there.
3. **Show exactly what you ran.** For every command or script you run, show the exact command line (and the script contents, or its path, if it is new), plus the key output. Alex should be able to copy-paste it and reproduce the result.
4. **Keep Alex in the loop and explain your reasoning.** Say what you are about to do and why before doing it; say what you found and what you think it means afterwards. Flag assumptions, uncertainties, and anything surprising in the numbers. Explain physics or code reasoning rather than just giving conclusions; Alex enjoys learning.
5. **Never touch outputs and inputs silently.** Don't overwrite or delete anything under `output/`, `run/`, or the sample directories, and don't regenerate calibration constants (`--calibrate` rewrites `maia_doublets/calibs/*.json`) without asking.
6. **No git commits or pushes** unless Alex asks.

## What this repo is

A study of whether a MAIA-like muon collider detector with **double-layers in the tracker** (CMS Phase-2 style "stacked" modules) enables geometric, LST-style track finding in the presence of beam-induced background (BIB).

- MAIA paper: https://arxiv.org/abs/2502.00181
- ICHEP 2026 talk (preliminary result): https://indico.cern.ch/event/1522800/contributions/6982863/attachments/3321557/5948113/Tuna-2026-07-28-lstmucol-03.pdf
- Datasets are produced with https://github.com/alexandertuna/maia_datasets
- Geometries (k4geo fork):
  - **v06**: grouped double-layers in IT and OT: https://github.com/alexandertuna/k4geo/blob/KITP_10TeV_ITOT_doublets_v06/MuColl/MAIA/compact/MAIA_v0/
  - **v07**: equally-spaced double-layers in IT and OT: https://github.com/alexandertuna/k4geo/tree/KITP_10TeV_ITOT_doublets_v07/MuColl/MAIA/compact/MAIA_v0
  - Older versions (v01, v04, v05) still appear in the code and calibrations.
- `docs/lst-mucol-paper-outline.md`: outline of the planned journal paper (detector-design oriented, algorithm held fixed as a probe).

## Environment

Everything runs inside the MuColl apptainer image (see `README.md`):

```bash
apptainer run /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v2.9.8-amd64
setup_mucoll
source run/env_ubuntu24/bin/activate   # venv with `pip install -e .`
```

The package installs a CLI entry point `maia_doublets` (`maia_doublets/main.py:main`).

## Pipeline

`main.py` steers a chain in which each object is built only from objects that passed the previous step's cuts:

hits (`slcio.py`, `HitMaker`) → mini-doublets, MDs (`md.py`) → T2 (`t2.py`, 2 double-layers / 4 layers) → T4 (`t4.py`, 8 layers) → T8 (`t8.py`, 16 layers) → plots (`plot.py`) and cutflow (`--cutflow`, ndjson).

- `datasets.py`: hard-coded sample paths keyed by dataset name and `(geo, "sim")` or `(geo, "digi", smear)`. `-i/--input` overrides with a glob.
- `calib.py` + `calibs/{geo}_{sim|digi_XXum}.json`: cut thresholds (3σ signal intervals) measured with `--calibrate` on the muon-gun sample.
- `constants.py`: detector system IDs, collection names, units, physics constants (B = 5 T).
- Signal vs background is inferred from input filenames (`muonGun`/`pionGun` → signal). Background runs apply all cuts by default; signal runs apply only the cuts requested with `--cut-mds/--cut-t2s/--cut-t4s/--cut-t8s`. `--no-cuts` disables all.
- Intermediate DataFrames can be cached with `--write-*`/`--read-*` pickles; outputs conventionally go to `output/{geo}_{dataset}_digi_{smear}/`.
- `run/run.sh` is a scratch log of commands Alex has run (mostly commented out). It is a good reference for typical invocations.
- `maia_doublets/demos/`: standalone studies and checks. `paper/`: figure scripts for the paper (`python fig-*.py` → `fig-*.pdf`), sharing `style.py`.

Example (efficiency with all cuts, v07, 10 µm smearing):

```bash
maia_doublets --geo v07 --dataset muonGun_pT_2p0_2p1 --digi --smear 10um --cut-mds --cut-t2s --cut-t4s --cut-t8s --plot
```

## Code conventions (from README)

- **dl** (double-layer) = layer // 2 within a system; starts at 0 in both IT and OT.
- **gdl** (global double-layer) ignores the system: IT is 0–3, OT starts at 4.
- Objects are labeled by their **innermost** double-layer (a T2 spanning dl 0 and 1 has dl 0).
- **`i_mcp`**: truth-matched MCParticle index; `0xffffffff` (`NO_MCP`) means the hits do not share one truth particle.
- **`xx_ok`**: the object passed all selection cuts (e.g. `t4_ok`).
- Code style: pandas DataFrames throughout, `logging` for output, type hints on function signatures. Match the surrounding code.
