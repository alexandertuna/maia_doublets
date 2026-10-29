#!/usr/bin/env bash
#
# Find the overall efficiency, applying all cuts.
#
# Run like:
# > ./run_efficiency.sh
# > MAIA_DETECTOR=v07 MAIA_RESOLUTION=50um ./run_efficiency.sh
#
set -eo pipefail
cd "$(dirname "$0")"
source env.sh

GEO=${MAIA_DETECTOR:-v06}
SMEAR=${MAIA_RESOLUTION:-10um}
DATASET=${MAIA_DATASET:-muonGun_pT_2p0_2p1}

set -x
maia_doublets \
  --geo ${GEO} \
  --dataset ${DATASET} \
  --digi \
  --smear ${SMEAR} \
  --cut-mds \
  --cut-t2s \
  --cut-t4s \
  --cut-t8s \
  --plot
