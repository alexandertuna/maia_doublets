#!/usr/bin/env bash
#
# Write MDs to disk without cuts,
# without stressing the memory by writing all layers at once.
#
# Run like:
# > ./run_mds.sh
# > MAIA_DETECTOR=v07 MAIA_RESOLUTION=50um ./run_mds.sh
#
set -eo pipefail
cd "$(dirname "$0")"
source env.sh

GEO=${MAIA_DETECTOR:-v06}
SMEAR=${MAIA_RESOLUTION:-10um}
DATASET=${MAIA_DATASET:-neutrinoGun10}
PKL_DIR=../output/${GEO}_${DATASET}_digi_${SMEAR}
mkdir -p ${PKL_DIR}

CMD="maia_doublets \
  --geo ${GEO} \
  --dataset ${DATASET} \
  --digi \
  --smear ${SMEAR} \
  --no-cuts \
  --stop-after-mds"

set -x
${CMD} --write-mds ${PKL_DIR}/mds_ITB0_ITB1.pkl --cutflow ${PKL_DIR}/cutflow_ITB0_ITB1.ndjson --layers ITB0 ITB1
${CMD} --write-mds ${PKL_DIR}/mds_ITB2_ITB3.pkl --cutflow ${PKL_DIR}/cutflow_ITB2_ITB3.ndjson --layers ITB2 ITB3
${CMD} --write-mds ${PKL_DIR}/mds_ITB4_ITB5.pkl --cutflow ${PKL_DIR}/cutflow_ITB4_ITB5.ndjson --layers ITB4 ITB5
${CMD} --write-mds ${PKL_DIR}/mds_ITB6_ITB7.pkl --cutflow ${PKL_DIR}/cutflow_ITB6_ITB7.ndjson --layers ITB6 ITB7
${CMD} --write-mds ${PKL_DIR}/mds_OTB0_OTB1.pkl --cutflow ${PKL_DIR}/cutflow_OTB0_OTB1.ndjson --layers OTB0 OTB1
${CMD} --write-mds ${PKL_DIR}/mds_OTB2_OTB3.pkl --cutflow ${PKL_DIR}/cutflow_OTB2_OTB3.ndjson --layers OTB2 OTB3
${CMD} --write-mds ${PKL_DIR}/mds_OTB4_OTB5.pkl --cutflow ${PKL_DIR}/cutflow_OTB4_OTB5.ndjson --layers OTB4 OTB5
${CMD} --write-mds ${PKL_DIR}/mds_OTB6_OTB7.pkl --cutflow ${PKL_DIR}/cutflow_OTB6_OTB7.ndjson --layers OTB6 OTB7
