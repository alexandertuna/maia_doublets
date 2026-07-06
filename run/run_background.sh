GEO="v01"
SMEAR="10um"
DATA_DIR="/ceph/users/atuna/work/maia/maia_noodling/samples/${GEO}/neutrinoGun/${SMEAR}"
# neutrinoGun_digi_3.slcio

for IT in $(seq 0 9); do

    maia_doublets \
    -i ${DATA_DIR}/neutrinoGun_digi_${IT}.slcio \
    --geo ${GEO} \
    --digi \
    --smear ${SMEAR} \
    --fast-mds \
    --write-mcps ${GEO}_background100_digi_${SMEAR}/mcps_${IT}.pkl \
    --write-simhits ${GEO}_background100_digi_${SMEAR}/simhits_${IT}.pkl \
    --write-mds ${GEO}_background100_digi_${SMEAR}/mds_${IT}.pkl \
    &> ${GEO}_background100_digi_${SMEAR}/log_${IT}.txt

done
