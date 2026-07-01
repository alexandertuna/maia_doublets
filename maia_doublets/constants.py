import numpy as np

MCP_PKL = "mcps.pkl"
SIMHIT_PKL = "simhits.pkl"
DOUBLET_PKL = "doublets.pkl"

MUON = 13
MUON_NEUTRINO = 14
PARTICLES_OF_INTEREST = [
    MUON,
    MUON_NEUTRINO,
]

SPEED_OF_LIGHT = 299.792458  # mm/ns

SIGNAL = "muonGun"
NO_MCP = np.uint32(0xffffffff)

# https://github.com/AIDASoft/podio/blob/master/include/podio/ObjectID.h
PODIO_NO_MCP = -2

ONE_MM = 1.0
ZERO_POINT_ZERO_ONE_MM = 0.01
ONE_GEV = 1.0
ONE_POINT_FIVE_GEV = 1.5

MM_TO_CM = 0.1
CM_TO_MM = 10.0
BYTE_TO_MB = 1e-6
MEV_TO_GEV = 1e-3

CODE = "/ceph/users/atuna/work/maia"
XML = f"{CODE}/k4geo/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml"


EPSILON = 1e-6
MCPARTICLE = "MCParticle"
MAGNETIC_FIELD = 5.0 # Tesla

BARREL_TRACKER_MAX_ETA = 0.65
BARREL_TRACKER_MAX_RADIUS = 1446.0

OUTSIDE_BOUNDS = 0
INSIDE_BOUNDS = 1
UNDEFINED_BOUNDS = 2
POSSIBLE_BOUNDS = [OUTSIDE_BOUNDS, INSIDE_BOUNDS, UNDEFINED_BOUNDS]
BOUNDS = {
    OUTSIDE_BOUNDS: "outside",
    INSIDE_BOUNDS: "inside",
    UNDEFINED_BOUNDS: "undefined",
}

MIN_COSTHETA = 0.0
MIN_SIMHIT_PT_FRACTION = 0.7
MAX_TIME = 3.0 # in ns

INNER_TRACKER_BARREL_COLLECTION = "InnerTrackerBarrelCollection"
OUTER_TRACKER_BARREL_COLLECTION = "OuterTrackerBarrelCollection"

INNER_TRACKER_BARREL_RELATIONS = "IBTrackerHitsRelations"
OUTER_TRACKER_BARREL_RELATIONS = "OBTrackerHitsRelations"

INNER_TRACKER_BARREL_HITS = "IBTrackerHits"
OUTER_TRACKER_BARREL_HITS = "OBTrackerHits"

NOT_USED = 0
VERTEX_TRACKER_BARREL = 1
VERTEX_TRACKER_ENDCAP = 2
INNER_TRACKER_BARREL = 3
INNER_TRACKER_ENDCAP = 4
OUTER_TRACKER_BARREL = 5
OUTER_TRACKER_ENDCAP = 6

N_LAYERS_IN_T4 = 8

CUT_MISSING = np.zeros((10, 10))

#
# np.arrays below are indexed by (system, doublelayer)
# i.e.:
#   [dl0, dl1, dl2, dl3], # NOT_USED
#   [dl0, dl1, dl2, dl3], # VERTEX_TRACKER_BARREL
#   [dl0, dl1, dl2, dl3], # VERTEX_TRACKER_ENDCAP
#   [dl0, dl1, dl2, dl3], # INNER_TRACKER_BARREL
#   [dl0, dl1, dl2, dl3], # INNER_TRACKER_ENDCAP
#   [dl0, dl1, dl2, dl3], # OUTER_TRACKER_BARREL
#   [dl0, dl1, dl2, dl3], # OUTER_TRACKER_ENDCAP
#

REQ_PASSTHROUGH = "no cuts"
REQ_RZ = "rz req."
REQ_XY = "xy req."
REQ_RZ_XY = "both req."
DOUBLET_REQS = [
    REQ_PASSTHROUGH,
    REQ_XY,
    REQ_RZ,
    REQ_RZ_XY,
]

T2_REQ_DR_POS = "dr req."
T2_REQ_DZ_POS = "dz req."
T2_REQ_XY_ANG = "dtheta xy req."
T2_REQ_XY_CHI2 = "xy chi2 req."
T2_REQ_RZ_ANG = "dtheta rz req."
T2_REQ_ALL = "all reqs"
T2_REQS = [
    REQ_PASSTHROUGH,
    T2_REQ_DR_POS,
    T2_REQ_DZ_POS,
    # T2_REQ_XY_ANG,
    T2_REQ_XY_CHI2,
    T2_REQ_RZ_ANG,
    T2_REQ_ALL,
]

NICKNAMES = {
    INNER_TRACKER_BARREL: "ITB",
    OUTER_TRACKER_BARREL: "OTB",
}
NICKNAME_TO_SYSTEM = {
    value: key for key, value in NICKNAMES.items()
}

# For converting system layer to global layer
LAYER_OFFSET = np.array([
    0, # NOT_USED
    0, # VERTEX_TRACKER_BARREL
    0, # VERTEX_TRACKER_ENDCAP
    0, # INNER_TRACKER_BARREL
    0, # INNER_TRACKER_ENDCAP
    8, # OUTER_TRACKER_BARREL
    0, # OUTER_TRACKER_ENDCAP
])

# Do not change this easily
# Its measured for 2 GeV muons
# We multiply by 2 to speed up the processing (fewer slices)
DETECTOR_MAX_PHI = np.pi
DETECTOR_MAX_ETA = 2.5 # Must include all background hits
MAX_T2_DPHI = np.array([
    2 * DETECTOR_MAX_PHI,
    2 * DETECTOR_MAX_PHI,
    2 * DETECTOR_MAX_PHI,
    0.10 * 2,
    2 * DETECTOR_MAX_PHI,
    0.10 * 2,
    2 * DETECTOR_MAX_PHI,
])
MAX_T2_DETA = np.array([
    2 * DETECTOR_MAX_ETA,
    2 * DETECTOR_MAX_ETA,
    2 * DETECTOR_MAX_ETA,
    0.01 * 2,
    2 * DETECTOR_MAX_ETA,
    0.01 * 2,
    2 * DETECTOR_MAX_ETA,
])
N_T2_PHI_SLICES = (2 * DETECTOR_MAX_PHI / MAX_T2_DPHI).astype(int)
N_T2_ETA_SLICES = (2 * DETECTOR_MAX_ETA / MAX_T2_DETA).astype(int)
MAX_T4_DPHI = 0.25 * 2
MAX_T4_DETA = 0.025 * 2
N_T4_PHI_SLICES = int(2 * DETECTOR_MAX_PHI / MAX_T4_DPHI)
N_T4_ETA_SLICES = int(2 * DETECTOR_MAX_ETA / MAX_T4_DETA)
