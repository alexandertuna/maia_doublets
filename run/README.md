# Run directory

This is how I run the code:

```bash
# Go to the run/ directory
cd run/

# Create a virtual environment (once)
python -m venv env_ubuntu24

# Set up the MuColl stack and activate the environment (on each login)
source env.sh

# pip install the repo
pip install -e ..

# Run! One script per task; defaults can be changed with environment variables
./run_plot_signal.sh
MAIA_DETECTOR=v07 MAIA_RESOLUTION=50um ./run_efficiency.sh
```

| Script | What it does | Default `MAIA_DATASET` |
|---|---|---|
| `run_calibrate.sh` | Measure cut thresholds, rewriting `maia_doublets/calibs/*.json` | `muonGun_pT_2p0_2p1` |
| `run_plot_signal.sh` | Write pickles, cutflow, and plots | `muonGun_pT_2p0_2p1` |
| `run_plot_background.sh` | Write pickles, cutflow, and plots (with `--fast-mds`) | `neutrinoGun` |
| `run_mds.sh` | Write uncut MDs, two layers at a time | `neutrinoGun10` |
| `run_efficiency.sh` | Apply all cuts and make plots | `muonGun_pT_2p0_2p1` |

Other environment variables: `MAIA_DETECTOR` (default `v06`) and `MAIA_RESOLUTION` (default `10um`).

Pickles and cutflows go to `output/${MAIA_DETECTOR}_${MAIA_DATASET}_digi_${MAIA_RESOLUTION}/`.
