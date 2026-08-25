#
# Run like:
# > source run.sh
#
OUTPUT_DIR=$(dirname "${BASH_SOURCE[0]}")/../output
mkdir -p ${OUTPUT_DIR}

# How to calibrate cut thresholds (3sigma intervals)
# maia_doublets --geo v01 --signal --digi --smear 10um --calibrate
# maia_doublets --geo v05 --signal --digi --smear 10um --calibrate
# maia_doublets --geo v06 --signal --digi --smear 10um --calibrate
# maia_doublets --geo v07 --signal --digi --smear 10um --calibrate

# How to make plots
# GEO="v07"
# SMEAR="10um"
# PKL_DIR=${OUTPUT_DIR}/${GEO}_signal_digi_${SMEAR}
# mkdir -p ${PKL_DIR}
# maia_doublets \
#     --geo ${GEO} \
#     --signal \
#     --digi \
#     --smear ${SMEAR} \
#     --write-hits ${PKL_DIR}/hits.pkl \
#     --write-mcps ${PKL_DIR}/mcps.pkl \
#     --write-mds ${PKL_DIR}/mds.pkl \
#     --write-t2s ${PKL_DIR}/t2s.pkl \
#     --write-t4s ${PKL_DIR}/t4s.pkl \
#     --write-t8s ${PKL_DIR}/t8s.pkl \
#     --plot

# How to write background MDs to disk without cuts
#   and without stressing the memory by writing all layers at once
GEO="v07"
SMEAR="10um"
PKL_DIR=${OUTPUT_DIR}/${GEO}_neutrinoGun10_digi_${SMEAR}
mkdir -p ${PKL_DIR}
CMD="maia_doublets \
    --geo ${GEO} \
    --neutrinoGun10 \
    --digi \
    --smear ${SMEAR} \
    --no-cuts \
    --stop-after-mds"
# time ${CMD} --write-mds ${PKL_DIR}/mds_ITB0_ITB1.pkl --layers ITB0 ITB1
# time ${CMD} --write-mds ${PKL_DIR}/mds_ITB2_ITB3.pkl --layers ITB2 ITB3
# time ${CMD} --write-mds ${PKL_DIR}/mds_ITB4_ITB5.pkl --layers ITB4 ITB5
# time ${CMD} --write-mds ${PKL_DIR}/mds_ITB6_ITB7.pkl --layers ITB6 ITB7
# time ${CMD} --write-mds ${PKL_DIR}/mds_OTB0_OTB1.pkl --layers OTB0 OTB1
# time ${CMD} --write-mds ${PKL_DIR}/mds_OTB2_OTB3.pkl --layers OTB2 OTB3
# time ${CMD} --write-mds ${PKL_DIR}/mds_OTB4_OTB5.pkl --layers OTB4 OTB5
# time ${CMD} --write-mds ${PKL_DIR}/mds_OTB6_OTB7.pkl --layers OTB6 OTB7

# How to find the overall efficiency
# maia_doublets --geo v05 --signal --digi --smear 10um --cut-mds --cut-t2s --cut-t4s --cut-t8s --plot
# maia_doublets --geo v06 --signal --digi --smear 10um --cut-mds --cut-t2s --cut-t4s --cut-t8s --plot
# maia_doublets --geo v07 --signal --digi --smear 10um --cut-mds --cut-t2s --cut-t4s --cut-t8s --plot

# How to plot the 0-10 GeV pionGun sample
# maia_doublets --geo v01 -i "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/pionGun_pT_0_10/10um/pionGun_pT_0_10_digi_3*" --digi --smear 10um --plot
# maia_doublets --geo v01 -i "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/pionGun_pT_0_10/10um/pionGun_pT_0_10_digi_3*" --digi --smear 10um --cut-mds --cut-t2s --cut-t4s --cut-t8s --plot

# How to plot the 0-10 GeV muonGun sample
# maia_doublets --geo v01 -i "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_0_10/10um/muonGun_pT_0_10_digi_3*" --digi --smear 10um --cut-mds --cut-t2s --cut-t4s --cut-t8s --plot

# How to pickle the 0-10 GeV muonGun sample
# GEO="v01"
# SMEAR="10um"
# PKL_DIR=${OUTPUT_DIR}/${GEO}_muonGun_pT_0_10_digi_${SMEAR}
# maia_doublets \
#     --geo ${GEO} \
#     -i "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_0_10/10um/muonGun_pT_0_10_digi_3*" \
#     --digi \
#     --smear ${SMEAR} \
#     --write-hits ${PKL_DIR}/hits.pkl \
#     --write-mcps ${PKL_DIR}/mcps.pkl \
#     --write-mds ${PKL_DIR}/mds.pkl \
#     --write-t2s ${PKL_DIR}/t2s.pkl \
#     --write-t4s ${PKL_DIR}/t4s.pkl \
#     --write-t8s ${PKL_DIR}/t8s.pkl \
#     --cut-mds \
#     --cut-t2s \
#     --cut-t4s \
#     --cut-t8s


# How to run background neutrinoGun
GEO="v06"
SMEAR="10um"
PKL_DIR=${OUTPUT_DIR}/${GEO}_neutrinoGun_digi_${SMEAR}
mkdir -p ${PKL_DIR}
maia_doublets \
  --geo ${GEO} \
  --neutrinoGun \
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

# How to make background dz plot for OTB L01 only
# GEO="v01"
# SMEAR="10um"
# maia_doublets \
#   --geo ${GEO} \
#   --digi \
#   --smear ${SMEAR} \
#   --layers OTB0 OTB1 \
#   --cutflow "" \
#   --plot \
#   --signal
# #   --background10

# How to make the scatter plot demo
# time python ../maia_doublets/demos/demo_scatter2d.py \
#     --hits v01_neutrinoGun_digi_10um/hits_0.pkl \
#     --mds v01_neutrinoGun_digi_10um/mds_0.pkl \
#     --t2s v01_neutrinoGun_digi_10um/t2s_0.pkl \
#     --t4s v01_neutrinoGun_digi_10um/t4s_0.pkl \
#     --t8s v01_neutrinoGun_digi_10um/t8s_0.pkl \
#     --signal v01_signal_digi_10um/t8s.pkl \
#     --output tmp.png
