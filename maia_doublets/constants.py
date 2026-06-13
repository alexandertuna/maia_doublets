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

MD_DZ_CUT = {
    ("v01", "sim"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0.62, 1.02, 5.49, 7.37],
        [0, 0, 0, 0],
        [22, 29, 112, 137],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "00um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0.62, 1.02, 5.49, 7.37],
        [0, 0, 0, 0],
        [22, 29, 112, 137],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "10um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 13, 15],
        [0, 0, 0, 0],
        [31, 38, 120, 146],
        [0, 0, 0, 0],
    ]),
    ("v05", "sim"): np.array([
        22, # mm # doublelayer 0
        43, # mm # doublelayer 1
        79, # mm # doublelayer 2
        137, # mm # doublelayer 3
    ]),
}
MD_DR_CUT = {
    ("v01", "sim"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [7.2, 12.1, 102.0, 118.5],
        [0, 0, 0, 0],
        [260, 313, 718, 806],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "00um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [7.2, 12.1, 102.0, 118.5],
        [0, 0, 0, 0],
        [260, 313, 718, 806],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "10um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 108, 125],
        [0, 0, 0, 0],
        [267, 320, 723, 809],
        [0, 0, 0, 0],
    ]),
    ("v05", "sim"): np.array([
        260, # mm # doublelayer 0
        408, # mm # doublelayer 1
        589, # mm # doublelayer 2
        804, # mm # doublelayer 3
    ]),
}

T2_DZ_CUT = {
    ("v01", "sim"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0.6, 0, 5.8, 0],
        [0, 0, 0, 0],
        [24, 0, 120, 0],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "00um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0.6, 0, 5.8, 0],
        [0, 0, 0, 0],
        [24, 0, 120, 0],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "10um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 6.1, 0],
        [0, 0, 0, 0],
        [24, 0, 120, 0],
        [0, 0, 0, 0],
    ]),
    ("v05", "sim"): np.array([
        24, # mm # doublelayer 0
        0.0, # mm # doublelayer 1
        120, # mm # doublelayer 2
    ]),
}
T2_DR_CUT = {
    ("v01", "sim"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [8.8, 0, 107.9, 0],
        [0, 0, 0, 0],
        [281, 0, 757, 0],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "00um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [8.8, 0, 107.9, 0],
        [0, 0, 0, 0],
        [281, 0, 757, 0],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "10um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 108, 0],
        [0, 0, 0, 0],
        [282, 0, 753, 0],
        [0, 0, 0, 0],
    ]),
    ("v05", "sim"): np.array([
        281, # mm # doublelayer 0
        0.0, # mm # doublelayer 1
        757, # mm # doublelayer 2
    ]),
}
T2_DTHETA_RZ_CUT = {
    ("v01", "sim"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0.00370, 0, 0.00445, 0],
        [0, 0, 0, 0],
        [0.0075, 0.0, 0.0130, 0.0],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "00um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0.00370, 0, 0.00445, 0],
        [0, 0, 0, 0],
        [0.0075, 0.0, 0.0130, 0.0],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "10um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0.02738, 0],
        [0, 0, 0, 0],
        [0.02838, 0, 0.03067, 0],
        [0, 0, 0, 0],
    ]),
    ("v05", "sim"): np.array([
        0.0075, # doublelayer 0
        0.0000, # doublelayer 1
        0.0130, # doublelayer 2
    ]),
}
T2_CHI2_XY_CUT = {
    ("v01", "sim"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0.04, 0.0, 0.04, 0.0],
        [0, 0, 0, 0],
        [0.04, 0.0, 0.04, 0.0],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "00um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0.04, 0.0, 0.04, 0.0],
        [0, 0, 0, 0],
        [0.04, 0.0, 0.04, 0.0],
        [0, 0, 0, 0],
    ]),
    ("v01", "digi", "10um"): np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0.00431, 0],
        [0, 0, 0, 0],
        [0.00403, 0, 0.00347, 0],
        [0, 0, 0, 0],
    ]),
    ("v05", "sim"): np.array([
        0.040, # doublelayer 0
        0.000, # doublelayer 1
        0.040, # doublelayer 2
    ]),
}
T4_DZ_CUT = {
    ("v01", "sim"): np.array([
        54.0, # gdoublelayer 0
        54.0, # gdoublelayer 1
        54.0, # gdoublelayer 2
        54.0, # gdoublelayer 3
        54.0, # gdoublelayer 4
        54.0, # gdoublelayer 5
        54.0, # gdoublelayer 6
        54.0, # gdoublelayer 7
    ]),
    ("v01", "digi", "00um"): np.array([
        54.0, # gdoublelayer 0
    ]),
    ("v01", "digi", "10um"): np.array([
        0.0,
        0.0,
        11.9,
        0.0,
        53.2,
        0.0,
        0.0,
        0.0,
    ]),
    ("v05", "sim"): np.array([
        0.0, # doublelayer 0
    ]),
}
T4_DR_CUT = {
    ("v01", "sim"): np.array([
        459.0, # gdoublelayer 0
        459.0, # gdoublelayer 1
        459.0, # gdoublelayer 2
        459.0, # gdoublelayer 3
        459.0, # gdoublelayer 4
        459.0, # gdoublelayer 5
        459.0, # gdoublelayer 6
        459.0, # gdoublelayer 7
    ]),
    ("v01", "digi", "00um"): np.array([
        459.0, # gdoublelayer 0
    ]),
    ("v01", "digi", "10um"): np.array([
        0.0,
        0.0,
        174.7,
        0.0,
        459.0,
        0.0,
        0.0,
        0.0,
    ]),
    ("v05", "sim"): np.array([
        0.0, # gdoublelayer 0
    ]),
}
T4_DTHETA_RZ_CUT = {
    ("v01", "sim"): np.array([
        0.053, # gdoublelayer 0
        0.053, # gdoublelayer 1
        0.053, # gdoublelayer 2
        0.053, # gdoublelayer 3
        0.053, # gdoublelayer 4
        0.053, # gdoublelayer 5
        0.053, # gdoublelayer 6
        0.053, # gdoublelayer 7
    ]),
    ("v01", "digi", "00um"): np.array([
        0.053, # gdoublelayer 0
    ]),
    ("v01", "digi", "10um"): np.array([
        0.0,
        0.0,
        0.01834,
        0.0,
        0.05294,
        0.0,
        0.0,
        0.0,
    ]),
    ("v05", "sim"): np.array([
        0.0, # gdoublelayer 0
    ]),
}
T4_CHI2_XY_CUT = {
    ("v01", "sim"): np.array([
        0.352, # gdoublelayer 0
        0.352, # gdoublelayer 1
        0.352, # gdoublelayer 2
        0.352, # gdoublelayer 3
        0.352, # gdoublelayer 4
        0.352, # gdoublelayer 5
        0.352, # gdoublelayer 6
        0.352, # gdoublelayer 7
    ]),
    ("v01", "digi", "00um"): np.array([
        0.352, # gdoublelayer 0
    ]),
    ("v01", "digi", "10um"): np.array([
        0.0,
        0.0,
        0.09388,
        0.0,
        0.42751,
        0.0,
        0.0,
        0.0,
    ]),
    ("v05", "sim"): np.array([
        0.0, # gdoublelayer 0
    ]),
}

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
