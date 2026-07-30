"""
Print the area of our silicon tracker detectors.
"""

import numpy as np
import dd4hep
import DDRec
import logging
logger = logging.getLogger(__name__)


XML = "/ceph/users/atuna/work/maia/k4geo/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml"
DETECTOR_NAMES = [
    "InnerTrackerBarrel",
    "OuterTrackerBarrel",
]
CM_TO_M = 0.01
CM2_TO_M2 = CM_TO_M * CM_TO_M

def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    # initialize detector and surface manager
    detector = dd4hep.Detector.getInstance()
    detector.fromCompact(XML)
    surfman = DDRec.SurfaceManager(detector)
    dets = {name: detector.detector(name) for name in DETECTOR_NAMES}
    maps = {name: surfman.map(det.name()) for name, det in dets.items()}

    # initialize area counters
    area_total = {}
    area_per_layer = {}

    # measure area
    for name, surfmap in maps.items():
        logger.info(f"Measuring area of {len(surfmap)} surfaces for {name} ... ")
        area_total[name] = 0
        area_per_layer[name] = {}
        for i_surf, surf_pair in enumerate(surfmap):
            id = surf_pair.first
            surf = surf_pair.second
            system, side, layer, module, sensor = decode(id)
            if layer not in area_per_layer[name]:
                area_per_layer[name][layer] = 0
            length_u = surf.length_along_u() # cm
            length_v = surf.length_along_v() # cm
            area = length_u * length_v
            area_total[name] += area
            area_per_layer[name][layer] += area

    # announce
    for name in DETECTOR_NAMES:
        logger.info(f"Total area of {name}: {area_total[name]*CM2_TO_M2:.6f} m^2")
        for layer, area in area_per_layer[name].items():
            logger.info(f"  Layer {layer}: {area*CM2_TO_M2:.6f} m^2")


def decode(cellid: int):
    """
    This is valid for my geometry versions v01 and v05!
    """
    system = np.right_shift(cellid, 0) & 0b1_1111
    side = np.right_shift(cellid, 5) & 0b11
    layer = np.right_shift(cellid, 7) & 0b11_1111
    module = np.right_shift(cellid, 13) & 0b111_1111_1111
    sensor = np.right_shift(cellid, 24) & 0b1111_1111
    return system, side, layer, module, sensor

if __name__ == "__main__":
    main()
