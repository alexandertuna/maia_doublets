import numpy as np
import pandas as pd
import time
import logging
logger = logging.getLogger(__name__)

from constants import MD_DZ_CUT, MD_DR_CUT
from constants import MAGNETIC_FIELD, SPEED_OF_LIGHT
from constants import BYTE_TO_MB, MEV_TO_GEV, NO_MCP
from constants import N_LS_PHI_SLICES, N_LS_ETA_SLICES, DETECTOR_MAX_ETA, DETECTOR_MAX_PHI

class DoubletMaker:


    def __init__(self, geometry_version: str, sim: bool, smear: str, signal: bool, cut_doublets: bool, simhits: pd.DataFrame):
        self.signal = signal
        self.cut_doublets = cut_doublets
        key = (geometry_version, "sim") if sim else (geometry_version, "digi", smear)
        self.MD_DZ_CUT = MD_DZ_CUT[key]
        self.MD_DR_CUT = MD_DR_CUT[key]
        self.df = self.make_doublets(simhits)


    def make_doublets(self, df: pd.DataFrame) -> pd.DataFrame:
        logger.info("Making doublets ...")

        groupby_cols = [
            "file",
        ]
        if not self.signal:
            groupby_cols += [
                "i_event", # the event
                "simhit_system", # the system (IT, OT)
                "simhit_layer_div_2", # the double layer
                "simhit_module", # the phi-module
                # "simhit_sensor", # the z-sensor
            ]

        doublet_cols = [
            "file",
            "i_event", # the event
            "simhit_system", # the system (IT, OT)
            "simhit_layer_div_2", # the double layer
            "simhit_module", # the phi-module
            "simhit_sensor", # the z-sensor
        ]

        def make_doublets_from_group(group: pd.DataFrame) -> pd.DataFrame:

            lower_mask = group["simhit_layer_mod_2"] == 0
            upper_mask = group["simhit_layer_mod_2"] == 1

            # inner join to find doublets
            doublets = pd.merge(
                group[lower_mask],
                group[upper_mask],
                on=doublet_cols,
                how="inner",
                suffixes=("_lower", "_upper"),
            )

            # doublet feature: xy, dr at point of closest approach to origin
            slope_xy = np.divide(doublets["simhit_y_upper"] - doublets["simhit_y_lower"],
                                 doublets["simhit_x_upper"] - doublets["simhit_x_lower"])
            intercept_xy = doublets["simhit_y_lower"] - slope_xy * doublets["simhit_x_lower"]
            doublets["doublet_dr"] = np.abs(intercept_xy) / np.sqrt(1 + slope_xy**2)

            # doublet feature: rz
            slope_rz = np.divide(doublets["simhit_z_upper"] - doublets["simhit_z_lower"],
                                 doublets["simhit_r_upper"] - doublets["simhit_r_lower"])
            doublets["doublet_dz"] = doublets["simhit_z_lower"] - doublets["simhit_r_lower"] * slope_rz
            doublets["doublet_theta_rz"] = np.arctan(slope_rz)

            # record some numbers
            cutflow = {"all": len(doublets)}
            mask = {}

            # record some cut results
            dl = doublets["simhit_layer_div_2"]
            mask["dr"] = np.abs(doublets["doublet_dr"]) < self.MD_DR_CUT[dl]
            mask["dz"] = np.abs(doublets["doublet_dz"]) < self.MD_DZ_CUT[dl]
            mask["and"] = mask["dr"] & mask["dz"]
            doublets["doublet_ok"] = mask["and"].astype(bool)

            # remove as desired
            if self.cut_doublets:
                for cut in mask.keys():
                    cutflow[cut] = np.sum(mask[cut])
                doublets = doublets[mask["and"]]

            # rename some columns
            rename = {
                "simhit_system": "doublet_system",
                "simhit_layer_div_2": "doublet_doublelayer",
                "simhit_sensor": "doublet_sensor",
                "simhit_module": "doublet_module",
            }
            doublets = doublets.rename(columns=rename)

            # doublet feature, xy dphi
            phi_local = np.arctan2(doublets["simhit_y_upper"] - doublets["simhit_y_lower"],
                                   doublets["simhit_x_upper"] - doublets["simhit_x_lower"])
            phi_global = np.arctan2((doublets["simhit_y_lower"] + doublets["simhit_y_upper"]) / 2.0,
                                    (doublets["simhit_x_lower"] + doublets["simhit_x_upper"]) / 2.0)
            doublets["doublet_dphi"] = phi_local - phi_global
            doublets["doublet_dphi"] = (doublets["doublet_dphi"] + np.pi) % (2 * np.pi) - np.pi
            doublets["doublet_theta_xy"] = phi_local

            # doublet features: position
            doublets["doublet_r"] = (doublets["simhit_r_lower"] + doublets["simhit_r_upper"]) / 2
            doublets["doublet_z"] = (doublets["simhit_z_lower"] + doublets["simhit_z_upper"]) / 2
            doublets["doublet_x"] = (doublets["simhit_x_lower"] + doublets["simhit_x_upper"]) / 2
            doublets["doublet_y"] = (doublets["simhit_y_lower"] + doublets["simhit_y_upper"]) / 2
            doublets["doublet_phi"] = np.arctan2(doublets["doublet_y"], doublets["doublet_x"])
            doublets["doublet_theta"] = np.arctan2(doublets["doublet_r"], doublets["doublet_z"])
            doublets["doublet_eta"] = -np.log(np.tan(doublets["doublet_theta"] / 2))
            doublets["doublet_phi_slice"] = np.floor((doublets["doublet_phi"] + DETECTOR_MAX_PHI) / (2 * DETECTOR_MAX_PHI) * N_LS_PHI_SLICES).astype(np.int16)
            doublets["doublet_eta_slice"] = np.floor((doublets["doublet_eta"] + DETECTOR_MAX_ETA) / (2 * DETECTOR_MAX_ETA) * N_LS_ETA_SLICES).astype(np.int16)

            # guess charge from dphi:
            # positively charged particles have negative dphi, and vice versa
            doublets["doublet_q"] = (-1*np.sign(doublets["doublet_dphi"])).astype(np.int8)

            # pass-through the simhit positions
            for coord in ["x", "y", "r"]:
                doublets[f"doublet_{coord}_0"] = doublets[f"simhit_{coord}_lower"]
                doublets[f"doublet_{coord}_1"] = doublets[f"simhit_{coord}_upper"]

            # doublet feature: radius of circle composed of the two hits and the origin. R = abc/4K
            # then get pt from R
            circle_a = doublets["simhit_r_lower"]
            circle_b = doublets["simhit_r_upper"]
            circle_c = np.sqrt((doublets["simhit_x_upper"] - doublets["simhit_x_lower"])**2 +
                               (doublets["simhit_y_upper"] - doublets["simhit_y_lower"])**2)
            circle_K = 0.5 * np.abs(doublets["simhit_x_lower"] * doublets["simhit_y_upper"] -
                                    doublets["simhit_x_upper"] * doublets["simhit_y_lower"])
            doublets["doublet_circle_radius"] = np.divide(circle_a * circle_b * circle_c, 4.0 * circle_K)
            doublets["doublet_pt"] = SPEED_OF_LIGHT * MAGNETIC_FIELD * doublets["doublet_circle_radius"] * 1e-6
            doublets["doublet_qoverpt"] = doublets["doublet_q"] / doublets["doublet_pt"]

            # doublet feature: truth info
            mcp_ok = doublets["i_mcp_lower"] == doublets["i_mcp_upper"]
            doublets["i_mcp"] = doublets["i_mcp_lower"].where(mcp_ok, NO_MCP)
            if self.signal:
                doublets["doublet_first_exit"] = doublets["simhit_first_exit_lower"] & doublets["simhit_first_exit_upper"]
                for attr in [
                    "mcp_pt",
                    "mcp_eta",
                    "mcp_phi",
                    "mcp_pdg",
                    "mcp_q",
                    "mcp_vertex_r",
                    "mcp_vertex_z",
                    "mcp_qoverpt",
                ]:
                    doublets[attr] = doublets[f"{attr}_lower"].where(mcp_ok, 0)

            # drop columns which arent used downstream
            dropcols = ["i_mcp_lower", "i_mcp_upper"]
            dropcols.extend([col for col in doublets.columns if col.startswith("simhit_")])
            dropcols.extend([col for col in doublets.columns if col.startswith("mcp_") and col.endswith("_lower")])
            dropcols.extend([col for col in doublets.columns if col.startswith("mcp_") and col.endswith("_upper")])
            doublets.drop(columns=dropcols, inplace=True)

            return doublets, cutflow

        # group loop
        groups = df.groupby(groupby_cols)
        all_doublets, all_cutflows = [], []

        for i_group, (cols, group) in enumerate(groups):

            doublets, cutflow = make_doublets_from_group(group)

            all_doublets.append(doublets)
            all_cutflows.append(cutflow)

            if (self.signal and i_group % 100 == 0) or (not self.signal and i_group % 4 == 0):
                length = len(doublets)
                size = doublets.memory_usage(deep=True).sum() * BYTE_TO_MB
                logger.info(f"Processed group {i_group}/{len(groups)}, doublet size = {size:.1f} MB, n(doublets) = {length} ...")

        # concatenate doublets and cutflows
        logger.info(f"Concatenating doublets ...")
        doublets = pd.concat(all_doublets, ignore_index=True)
        cutflow = pd.DataFrame(all_cutflows)
        for col in cutflow.columns:
            logger.info(f"Doublets cutflow, {col}: {cutflow[col].sum()}")
        if len(doublets) == 0:
            raise ValueError("No doublets found in the DataFrame")

        # announcements
        logger.info(f"Total doublets: {len(doublets)}")
        logger.info(f"Total doublets size: {doublets.memory_usage(deep=True).sum() * BYTE_TO_MB:.1f} MB")
        counts = doublets.groupby(["doublet_system",
                                   "doublet_doublelayer"]).size()
        for (system, doublelayer), total in counts.items():
            logger.info(f"n(doublets) for system {system}, doublelayer {doublelayer}: {total}")

        return doublets
