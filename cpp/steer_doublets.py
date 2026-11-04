# Run the C++ maia_doublets track finding on digitised v06 or v07 samples with k4run,
# with Gaudi's timing auditors.
#
# The cuts come from maia_doublets/calibs/{detector}_digi_{resolution}.json, unless --calibFile is given.
#
# The processor drops hits outside the time window as it reads them:
#   --timeWindow wide: keep the digitised window, [-5, 15] ns
#   --timeWindow maia: MAIA's IT/OT window, [-0.18, 0.3] ns
#   --timeWindow none: no time cut, as in the Python code
# Both sensors of each double layer are kept.
#
# Example:
# > k4run steer_doublets.py -n 10 --detector v06 --inputFile muonGun_pT_2p0_2p1_digi_300.slcio
from Configurables import ApplicationMgr, AuditorSvc
from Gaudi.Configuration import *

from Configurables import LcioEvent, EventDataSvc, MarlinProcessorWrapper

from k4FWCore.parseArgs import parser

import os

# Set by cpp/env.sh (inside k4run, __file__ is not this file)
CALIB_DIR = os.path.join(os.environ["MAIA_DOUBLETS_CPP_DIR"], "..", "maia_doublets", "calibs")

parser.add_argument("--inputFile", type=str, nargs="+", required=True, help="Input digitised slcio file(s)")
parser.add_argument("--timeWindow", type=str, default="maia", choices=["wide", "maia", "none"], help="Hit time window")
parser.add_argument("--detector", type=str, default="v06", choices=["v06", "v07"], help="Geometry, for the cut thresholds")
parser.add_argument("--resolution", type=str, default="10um", help="Hit smearing, for the cut thresholds (e.g. 00um, 10um)")
parser.add_argument("--calibFile", type=str, default="", help="Cut thresholds, instead of the ones for --detector and --resolution")
parser.add_argument("--noCuts", action="store_true", default=False, help="Keep objects that fail the cuts (for validation)")
parser.add_argument("--dumpDir", type=str, default="", help="Write every object to {mds,t2s,t4s,t8s}.csv here (for validation)")
the_args = parser.parse_args()

calib_file = the_args.calibFile or os.path.join(CALIB_DIR, f"{the_args.detector}_digi_{the_args.resolution}.json")
if not os.path.isfile(calib_file):
    raise FileNotFoundError(f"No calibration file {calib_file}")

TIME_WINDOWS = {
    # (IT/OT min, IT/OT max) [ns]
    "wide": ("-5.0", "15.0"),
    "maia": ("-0.18", "0.3"),
}

algList = []
evtsvc = EventDataSvc()

read = LcioEvent()
read.OutputLevel = INFO
read.Files = the_args.inputFile
algList.append(read)

DoubletTracking = MarlinProcessorWrapper("DoubletTracking")
DoubletTracking.OutputLevel = INFO
DoubletTracking.ProcessorType = "DoubleLayerTrackingProc"
DoubletTracking.Parameters = {
    "TrackerHitCollectionNames": ["IBTrackerHits", "OBTrackerHits"],
    "CalibFile": [calib_file],
    "ApplyCuts": ["false" if the_args.noCuts else "true"],
    "DumpDir": [the_args.dumpDir],
}
if the_args.timeWindow != "none":
    time_min, time_max = TIME_WINDOWS[the_args.timeWindow]
    DoubletTracking.Parameters["TimeMin"] = [time_min]
    DoubletTracking.Parameters["TimeMax"] = [time_max]
algList.append(DoubletTracking)

# Per-algorithm timing, printed at the end of the job:
#   AlgTimingAuditor: wall time (steady_clock)
#   ChronoAuditor:    CPU time (user and system) and wall time, via ChronoStatSvc
ApplicationMgr(TopAlg=algList,
               EvtSel='NONE',
               EvtMax=-1,
               ExtSvc=[evtsvc, AuditorSvc(Auditors=["AlgTimingAuditor", "ChronoAuditor"])],
               AuditAlgorithms=True,
               OutputLevel=INFO
               )
