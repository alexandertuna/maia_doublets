# How to calibrate cut thresholds (3sigma intervals)
# maia_doublets --geo v01 --signal --digi --smear 10um --calibrate

# How to make plots
GEO="v01"
SMEAR="10um"
OUTDIR=${GEO}_signal_digi_${SMEAR}
maia_doublets \
    --geo ${GEO} \
    --signal \
    --digi \
    --smear ${SMEAR} \
    --write-simhits ${OUTDIR}/hits.pkl \
    --write-mcps ${OUTDIR}/mcps.pkl \
    --write-mds ${OUTDIR}/mds.pkl \
    --write-t2s ${OUTDIR}/t2s.pkl \
    --write-t4s ${OUTDIR}/t4s.pkl \
    --write-t8s ${OUTDIR}/t8s.pkl \
    --plot

# How to find the overall efficiency
# maia_doublets --geo v01 --signal --digi --smear 10um --cut-mds --cut-t2s --cut-t4s --cut-t8s --plot

# How to plot the 0-10 GeV muonGun sample
# maia_doublets --geo v01 -i "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_0_10/10um/muonGun_pT_0_10_digi_3*" --digi --smear 10um --cut-mds --cut-t2s --cut-t4s --cut-t8s --plot

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
#   --write-simhits ${OUTDIR}/simhits.pkl \
#   --write-mds ${OUTDIR}/mds.pkl # --plot

