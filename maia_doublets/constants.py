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

INNER_TRACKER_BARREL = 3
OUTER_TRACKER_BARREL = 5

MD_DZ_CUT = {
    ("v01", "sim"): np.array([
        22, # mm # doublelayer 0
        29, # mm # doublelayer 1
        112, # mm # doublelayer 2
        137, # mm # doublelayer 3
    ]),
    ("v01", "digi", "00um"): np.array([
        22, # mm # doublelayer 0
        29, # mm # doublelayer 1
        112, # mm # doublelayer 2
        137, # mm # doublelayer 3
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
        260, # mm # doublelayer 0
        313, # mm # doublelayer 1
        718, # mm # doublelayer 2
        806, # mm # doublelayer 3
    ]),
    ("v01", "digi", "00um"): np.array([
        260, # mm # doublelayer 0
        313, # mm # doublelayer 1
        718, # mm # doublelayer 2
        806, # mm # doublelayer 3
    ]),
    ("v05", "sim"): np.array([
        260, # mm # doublelayer 0
        408, # mm # doublelayer 1
        589, # mm # doublelayer 2
        804, # mm # doublelayer 3
    ]),
}

LS_DZ_CUT = {
    ("v01", "sim"): np.array([
        24, # mm # doublelayer 0
        0.0, # mm # doublelayer 1
        120, # mm # doublelayer 2
    ]),
    ("v01", "digi", "00um"): np.array([
        24, # mm # doublelayer 0
        0.0, # mm # doublelayer 1
        120, # mm # doublelayer 2
    ]),
    ("v05", "sim"): np.array([
        24, # mm # doublelayer 0
        0.0, # mm # doublelayer 1
        120, # mm # doublelayer 2
    ]),
}
LS_DR_CUT = {
    ("v01", "sim"): np.array([
        281, # mm # doublelayer 0
        0.0, # mm # doublelayer 1
        757, # mm # doublelayer 2
    ]),
    ("v01", "digi", "00um"): np.array([
        281, # mm # doublelayer 0
        0.0, # mm # doublelayer 1
        757, # mm # doublelayer 2
    ]),
    ("v05", "sim"): np.array([
        281, # mm # doublelayer 0
        0.0, # mm # doublelayer 1
        757, # mm # doublelayer 2
    ]),
}
LS_DTHETA_RZ_CUT = {
    ("v01", "sim"): np.array([
        0.0075, # doublelayer 0
        0.0000, # 0.0453, # doublelayer 1
        0.0130, # doublelayer 2
    ]),
    ("v01", "digi", "00um"): np.array([
        0.0075, # doublelayer 0
        0.0000, # 0.0453, # doublelayer 1
        0.0130, # doublelayer 2
    ]),
    ("v05", "sim"): np.array([
        0.0075, # doublelayer 0
        0.0000, # doublelayer 1
        0.0130, # doublelayer 2
    ]),
}
LS_DTHETA_XY_CUT = {
    ("v01", "sim"): np.array([
        0.070, # doublelayer 0
        0.000, # 0.396, # doublelayer 1
        0.078, # doublelayer 2
    ]),
    ("v01", "digi", "00um"): np.array([
        0.070, # doublelayer 0
        0.000, # 0.396, # doublelayer 1
        0.078, # doublelayer 2
    ]),
    ("v05", "sim"): np.array([
        0.000, # doublelayer 0
        0.000, # doublelayer 1
        0.000, # doublelayer 2
    ]),
}
LS_CHI2_XY_CUT = {
    ("v01", "sim"): np.array([
        0.040, # doublelayer 0
        0.000, # doublelayer 1
        0.040, # doublelayer 2
    ]),
    ("v01", "digi", "00um"): np.array([
        0.040, # doublelayer 0
        0.000, # doublelayer 1
        0.040, # doublelayer 2
    ]),
    ("v05", "sim"): np.array([
        0.040, # doublelayer 0
        0.000, # doublelayer 1
        0.040, # doublelayer 2
    ]),
}
T4_DZ_CUT = {
    ("v01", "sim"): np.array([
        54.0, # doublelayer 0
    ]),
    ("v01", "digi", "00um"): np.array([
        54.0, # doublelayer 0
    ]),
    ("v05", "sim"): np.array([
        0.0, # doublelayer 0
    ]),
}
T4_DR_CUT = {
    ("v01", "sim"): np.array([
        459.0, # doublelayer 0
    ]),
    ("v01", "digi", "00um"): np.array([
        459.0, # doublelayer 0
    ]),
    ("v05", "sim"): np.array([
        0.0, # doublelayer 0
    ]),
}
T4_DTHETA_RZ_CUT = {
    ("v01", "sim"): np.array([
        0.053, # doublelayer 0
    ]),
    ("v01", "digi", "00um"): np.array([
        0.053, # doublelayer 0
    ]),
    ("v05", "sim"): np.array([
        0.0, # doublelayer 0
    ]),
}
T4_CHI2_XY_CUT = {
    ("v01", "sim"): np.array([
        0.352, # doublelayer 0
    ]),
    ("v01", "digi", "00um"): np.array([
        0.352, # doublelayer 0
    ]),
    ("v05", "sim"): np.array([
        0.0, # doublelayer 0
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

LS_REQ_DR_POS = "dr req."
LS_REQ_DZ_POS = "dz req."
LS_REQ_XY_ANG = "dtheta xy req."
LS_REQ_XY_CHI2 = "xy chi2 req."
LS_REQ_RZ_ANG = "dtheta rz req."
LS_REQ_ALL = "all reqs"
LS_REQS = [
    REQ_PASSTHROUGH,
    LS_REQ_DR_POS,
    LS_REQ_DZ_POS,
    # LS_REQ_XY_ANG,
    LS_REQ_XY_CHI2,
    LS_REQ_RZ_ANG,
    LS_REQ_ALL,
]

NICKNAMES = {
    INNER_TRACKER_BARREL: "ITB",
    OUTER_TRACKER_BARREL: "OTB",
}

# Do not change this easily
# Its measured for 2 GeV muons
# We multiply by 2 to speed up the processing (fewer slices)
DETECTOR_MAX_PHI = np.pi
DETECTOR_MAX_ETA = 2.5 # Must include all background hits
MAX_LS_DPHI = 0.10 * 2
MAX_LS_DETA = 0.01 * 2
N_LS_PHI_SLICES = int(2 * DETECTOR_MAX_PHI / MAX_LS_DPHI)
N_LS_ETA_SLICES = int(2 * DETECTOR_MAX_ETA / MAX_LS_DETA)
MAX_T4_DPHI = 0.25 * 2
MAX_T4_DETA = 0.025 * 2
N_T4_PHI_SLICES = int(2 * DETECTOR_MAX_PHI / MAX_T4_DPHI)
N_T4_ETA_SLICES = int(2 * DETECTOR_MAX_ETA / MAX_T4_DETA)
