This directory contains scripts which make plots for the paper.

A plot-making environment can be set up like:

```bash
apptainer run /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v2.9.8-amd64
Apptainer> setup_mucoll
```

A plot can be made like:

```bash
python fig-scaling.py
```

This will produce a plot named `fig-scaling.pdf`.

`style.py` provides a common style for plots.
