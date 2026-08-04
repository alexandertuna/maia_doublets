#
# Run like:
# > source run.sh
#
OUTPUT_DIR=$(dirname "${BASH_SOURCE[0]}")/../output
mkdir -p ${OUTPUT_DIR}

# How to calibrate cut thresholds (3sigma intervals)
# maia_doublets --geo v01 --signal --digi --smear 10um --calibrate

# How to make plots
GEO="v01"
SMEAR="10um"
PKL_DIR=${OUTPUT_DIR}/${GEO}_signal_digi_${SMEAR}
maia_doublets \
    --geo ${GEO} \
    --signal \
    --digi \
    --smear ${SMEAR} \
    --write-hits ${PKL_DIR}/hits.pkl \
    --write-mcps ${PKL_DIR}/mcps.pkl \
    --write-mds ${PKL_DIR}/mds.pkl \
    --write-t2s ${PKL_DIR}/t2s.pkl \
    --write-t4s ${PKL_DIR}/t4s.pkl \
    --write-t8s ${PKL_DIR}/t8s.pkl \
    --cut-mds \
    --cut-t2s \
    --cut-t4s \
    --cut-t8s \
    --plot

# How to find the overall efficiency
# maia_doublets --geo v01 --signal --digi --smear 10um --cut-mds --cut-t2s --cut-t4s --cut-t8s --plot

# How to plot the 0-10 GeV pionGun sample
# maia_doublets --geo v01 -i "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/pionGun_pT_0_10/10um/pionGun_pT_0_10_digi_3*" --digi --smear 10um --plot
# maia_doublets --geo v01 -i "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/pionGun_pT_0_10/10um/pionGun_pT_0_10_digi_3*" --digi --smear 10um --cut-mds --cut-t2s --cut-t4s --cut-t8s --plot

# How to plot the 0-10 GeV muonGun sample
# maia_doublets --geo v01 -i "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_0_10/10um/muonGun_pT_0_10_digi_3*" --digi --smear 10um --cut-mds --cut-t2s --cut-t4s --cut-t8s --plot

# How to pickle the 0-10 GeV muonGun sample
# GEO="v01"
# SMEAR="10um"
# OUTDIR=${GEO}_muonGun_pT_0_10_digi_${SMEAR}
# maia_doublets \
#     --geo ${GEO} \
#     -i "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_0_10/10um/muonGun_pT_0_10_digi_3*" \
#     --digi \
#     --smear ${SMEAR} \
#     --write-hits ${OUTDIR}/hits.pkl \
#     --write-mcps ${OUTDIR}/mcps.pkl \
#     --write-mds ${OUTDIR}/mds.pkl \
#     --write-t2s ${OUTDIR}/t2s.pkl \
#     --write-t4s ${OUTDIR}/t4s.pkl \
#     --write-t8s ${OUTDIR}/t8s.pkl \
#     --cut-mds \
#     --cut-t2s \
#     --cut-t4s \
#     --cut-t8s


# How to run background neutrinoGun
# GEO="v01"
# SMEAR="10um"
# OUTDIR=${GEO}_background100_digi_${SMEAR}
# maia_doublets \
#   --geo ${GEO} \
#   --background100 \
#   --digi \
#   --smear ${SMEAR} \
#   --fast-mds \
#   --cutflow ${OUTDIR}/cutflow.ndjson \
#   --write-mcps ${OUTDIR}/mcps.pkl \
#   --write-hits ${OUTDIR}/hits.pkl \
#   --write-mds ${OUTDIR}/mds.pkl # --plot

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
#     --hits v01_background100_digi_10um/hits_0.pkl \
#     --mds v01_background100_digi_10um/mds_0.pkl \
#     --t2s v01_background100_digi_10um/t2s_0.pkl \
#     --t4s v01_background100_digi_10um/t4s_0.pkl \
#     --t8s v01_background100_digi_10um/t8s_0.pkl \
#     --signal v01_signal_digi_10um/t8s.pkl \
#     --output tmp.png
