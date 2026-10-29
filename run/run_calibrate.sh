#!/usr/bin/env bash
#
# Measure cut thresholds (3sigma intervals) and write them to
# maia_doublets/calibs/${MAIA_DETECTOR}_digi_${MAIA_RESOLUTION}.json.
# The jsons are git controlled, so `git diff` shows whether they changed.
#
# Run like:
# > ./run_calibrate.sh
# > MAIA_DETECTOR=v07 MAIA_RESOLUTION=50um ./run_calibrate.sh
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
  --calibrate
