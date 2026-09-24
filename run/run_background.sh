#
# Run me like:
# > source run_background.sh
#
GEO="v06"
SMEAR="10um"
DATASET="neutrinoGun60"

DATA_DIR="/ceph/users/atuna/work/maia/maia_datasets/samples/${GEO}/${DATASET}/${SMEAR}"
OUTPUT_DIR=$(dirname "${BASH_SOURCE[0]}")/../output
PKL_DIR=${OUTPUT_DIR}/${GEO}_${DATASET}_digi_${SMEAR}

mkdir -p ${PKL_DIR}
echo "Output directory: ${PKL_DIR}"

for IT in $(seq 0 9); do

    INPUT=${DATA_DIR}/${DATASET}_digi_${IT}.slcio

    echo "Running ${INPUT} ..."
    maia_doublets \
    -i ${INPUT} \
    --geo ${GEO} \
    --digi \
    --smear ${SMEAR} \
    --fast-mds \
    --dataset ${DATASET} \
    --cutflow ${PKL_DIR}/cutflow_${IT}.ndjson \
    --write-mcps ${PKL_DIR}/mcps_${IT}.pkl \
    --write-hits ${PKL_DIR}/hits_${IT}.pkl \
    --write-mds ${PKL_DIR}/mds_${IT}.pkl \
    --write-t2s ${PKL_DIR}/t2s_${IT}.pkl \
    --write-t4s ${PKL_DIR}/t4s_${IT}.pkl \
    --write-t8s ${PKL_DIR}/t8s_${IT}.pkl \
    &> ${PKL_DIR}/log_${IT}.txt

done
