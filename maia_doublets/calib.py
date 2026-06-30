import json
import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)

from maia_doublets.constants import NO_MCP
from maia_doublets.constants import MUON
from maia_doublets.constants import ONE_POINT_FIVE_GEV
from maia_doublets.constants import BARREL_TRACKER_MAX_ETA
from maia_doublets.constants import ZERO_POINT_ZERO_ONE_MM


class MDCalibrator:

    def __init__(self, doublets: pd.DataFrame, calib_json: str) -> None:
        self.doublets = doublets
        self.calib_json = calib_json
        self.percentile = 99.7
        self.features = [
            "doublet_dz",
            "doublet_dr",
        ]
        self.system = "doublet_system"
        self.doublelayer = "doublet_doublelayer"
        self.detectable = "doublet_detectable"
        self.groupby = [
            self.system,
            self.doublelayer,
        ]
        self.calib = {feature: {} for feature in self.features}
        logger.info(f"Calibrating MDs {self.features}")
        logger.info(f"len(doublets) = {len(doublets)}")
        logger.info(f"Systems: {self.doublets[self.system].unique()}")
        logger.info(f"Doublelayers: {self.doublets[self.doublelayer].unique()}")


    def calibrate(self, update_calibration: bool = True) -> None:
        mask = (
            self.doublets[self.detectable] &
            (self.doublets["i_mcp"] != NO_MCP)
        )
        for feature in self.features:
            for (cols, group) in self.doublets[mask].groupby(self.groupby):
                (system, doublelayer) = [str(col) for col in cols]
                if system not in self.calib[feature]:
                    self.calib[feature][system] = {}
                perc = np.percentile(np.abs(group[feature]), self.percentile)
                self.calib[feature][system][doublelayer] = perc
        if update_calibration:
            self.update_calibration()


    def update_calibration(self) -> None:
        calib_dict = read_calibration(self.calib_json)
        calib_dict = update_calibration(calib_dict, self.calib)
        write_calibration(calib_dict, self.calib_json)


def read_calibration(calib_json: str) -> dict:
    try:
        with open(calib_json, "r") as fi:
            calib_dict = json.load(fi)
    except FileNotFoundError:
        calib_dict = {}
    return calib_dict


def update_calibration(old_calib: dict, new_calib: dict) -> dict:
    for feature in new_calib:
        if feature not in old_calib:
            old_calib[feature] = {}
        for system, doublelayer_dict in new_calib[feature].items():
            if system not in old_calib[feature]:
                old_calib[feature][system] = {}
            for doublelayer, perc in doublelayer_dict.items():
                old_calib[feature][system][doublelayer] = perc
    return old_calib


def write_calibration(calib_dict: dict, calib_json: str) -> None:
    with open(calib_json, "w") as fo:
        json.dump(calib_dict, fo, indent=4)
