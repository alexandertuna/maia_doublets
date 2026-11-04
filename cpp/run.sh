#!/usr/bin/env bash
#
# Run the C++ maia_doublets track finding on one v06 or v07 sample file, for timing.
# Per-algorithm timing and per-stage object counts are printed at the end of the log.
#
# Run inside the container, like:
# > ./cpp/run.sh
# > MAIA_DATASET=muonGun_pT_2p0_2p1 MAIA_FILE=300 MAIA_EVENTS=10 ./cpp/run.sh
# > MAIA_TIME_WINDOW=wide ./cpp/run.sh
# > MAIA_DETECTOR=v07 ./cpp/run.sh
#
set -eo pipefail
cd "$(dirname "$0")"

source ../run/env.sh
source env.sh

GEO=${MAIA_DETECTOR:-v06}
SMEAR=${MAIA_RESOLUTION:-10um}
DATASET=${MAIA_DATASET:-neutrinoGun}
FILE=${MAIA_FILE:-0}
EVENTS=${MAIA_EVENTS:-1}
WINDOW=${MAIA_TIME_WINDOW:-none}

# File names differ between datasets (e.g. neutrinoGun10 has neutrinoGun_digi_0.slcio)
DATA_DIR=/ceph/users/atuna/work/maia/maia_datasets/samples/${GEO}/${DATASET}/${SMEAR}
INPUT=$(ls ${DATA_DIR}/*_digi_${FILE}.slcio)

OUTPUT_DIR=../output/${GEO}_${DATASET}_digi_${SMEAR}
NAME=doublets_${WINDOW}_${FILE}
mkdir -p ${OUTPUT_DIR}

set -x
k4run steer_doublets.py \
  -n ${EVENTS} \
  --inputFile ${INPUT} \
  --timeWindow ${WINDOW} \
  --detector ${GEO} \
  --resolution ${SMEAR} \
  2>&1 | tee ${OUTPUT_DIR}/${NAME}.log
