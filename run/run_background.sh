GEO="v05"
SMEAR="10um"
DATA_DIR="/ceph/users/atuna/work/maia/maia_noodling/samples/${GEO}/neutrinoGun/${SMEAR}"
# neutrinoGun_digi_3.slcio

for IT in $(seq 10 19); do

    OUTDIR=${GEO}_background100_digi_${SMEAR}

    maia_doublets \
    -i ${DATA_DIR}/neutrinoGun_digi_${IT}.slcio \
    --geo ${GEO} \
    --digi \
    --smear ${SMEAR} \
    --fast-mds \
    --cutflow ${OUTDIR}/cutflow_${IT}.ndjson \
    --write-mcps ${OUTDIR}/mcps_${IT}.pkl \
    --write-simhits ${OUTDIR}/simhits_${IT}.pkl \
    --write-mds ${OUTDIR}/mds_${IT}.pkl \
    &> ${OUTDIR}/log_${IT}.txt

done
