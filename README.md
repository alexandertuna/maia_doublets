# maia_doublets

Set up the analysis environment like:

```bash
apptainer run /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-alma9:v2.9.8-amd64
apptainer> setup_mucoll
```

# Code conventions

There are a few conventions used throughout the code which can be helpful to know.

## Double-layer index

Double-layer (dl) refers to layer mod 2, and is used frequently throughout the code.

Double-layer starts from zero in the inner tracker and zero in the outer tracker.

## Global double-layer index

Global double-layer (gdl) refers to the double-layer index when ignoring the parent detector.

Put another way: global double-layer starts from zero in the inner tracker and four in the outer tracker.

## Objects are referenced by their innermost double-layer

For example, a T2 which includes hits from double-layer 0 and 1 will have a "double-layer" of 0.

The implication is a T2 with double-layer 0 naturally includes double-layer 1, since it's a T2.

## `xx_ok`

`xx_ok` is a flag which evaluates if an object has passed all the selection criteria.

For example, if `t4_ok` is True, that T4 has passed all the selection cuts.

## T2 are only built from MDs which pass all the MD cuts

The same is true for T4s built from T2s, and T8s built from T4s.

Keep this in mind when applying the `--cut-xx` flags.

# To-do items

This is a non-exhaustive list of things which I should address:

- Mixing "doublet" with "md"
- Mixing "linesegment and "ls" with "t2"
- Mixing global double layer with system, layer, double-layer
- Using double-layers at all? Is that silly?
