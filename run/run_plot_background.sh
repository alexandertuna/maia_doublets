#!/usr/bin/env bash
#
# Write pickles, cutflow, and plots for a background sample.
#
# Run like:
# > ./run_plot_background.sh
# > MAIA_DETECTOR=v07 MAIA_RESOLUTION=50um ./run_plot_background.sh
#
set -eo pipefail
cd "$(dirname "$0")"
source env.sh

GEO=${MAIA_DETECTOR:-v06}
SMEAR=${MAIA_RESOLUTION:-10um}
DATASET=${MAIA_DATASET:-neutrinoGun}
PKL_DIR=../output/${GEO}_${DATASET}_digi_${SMEAR}
mkdir -p ${PKL_DIR}

set -x
maia_doublets \
  --geo ${GEO} \
  --dataset ${DATASET} \
  --digi \
  --smear ${SMEAR} \
  --cutflow ${PKL_DIR}/cutflow.ndjson \
  --write-mcps ${PKL_DIR}/mcps.pkl \
  --write-hits ${PKL_DIR}/hits.pkl \
  --write-mds ${PKL_DIR}/mds.pkl \
  --write-t2s ${PKL_DIR}/t2s.pkl \
  --write-t4s ${PKL_DIR}/t4s.pkl \
  --write-t8s ${PKL_DIR}/t8s.pkl \
  --fast-mds --plot
