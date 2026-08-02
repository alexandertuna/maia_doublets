GEO="v01"
SMEAR="10um"
DATA_DIR="/ceph/users/atuna/work/maia/maia_noodling/samples/${GEO}/neutrinoGun/${SMEAR}"
# neutrinoGun_digi_3.slcio

for IT in $(seq 50 99); do

    OUTDIR=${GEO}_background100_digi_${SMEAR}

    maia_doublets \
    -i ${DATA_DIR}/neutrinoGun_digi_${IT}.slcio \
    --geo ${GEO} \
    --digi \
    --smear ${SMEAR} \
    --fast-mds \
    --cutflow ${OUTDIR}/cutflow_${IT}.ndjson \
    --write-mcps ${OUTDIR}/mcps_${IT}.pkl \
    --write-hits ${OUTDIR}/hits_${IT}.pkl \
    --write-mds ${OUTDIR}/mds_${IT}.pkl \
    --write-t2s ${OUTDIR}/t2s_${IT}.pkl \
    --write-t4s ${OUTDIR}/t4s_${IT}.pkl \
    --write-t8s ${OUTDIR}/t8s_${IT}.pkl \
    &> ${OUTDIR}/log_${IT}.txt

done
