"""
This module defines the T4Maker class, which creates a T4 from two T2s.
All T2s in layers 0-3 are considered for combination with T2s in layers 4-7.
The T2s are combined if they satisfy goodness criteria.
To avoid filling the memory with all possible combinations, combination is done with a groupby approach.

"""

import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)

from constants import BYTE_TO_MB, NO_MCP
from constants import T4_DZ_CUT, T4_DR_CUT, T4_DTHETA_RZ_CUT, T4_CHI2_XY_CUT
from constants import N_T4_PHI_SLICES

class T4Maker:

    def __init__(self, geometry_version: str, sim: bool, smear: str, t2s: pd.DataFrame, signal: bool, cut_t4s: bool):
        self.df = None
        self.signal = signal
        self.cut_t4s = cut_t4s
        self.t2s = t2s.copy()
        memory = self.t2s.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Making T4s. T2 dataframe size: {memory:.2f} MB")

        key = (geometry_version, "sim") if sim else (geometry_version, "digi", smear)
        self.T4_DZ_CUT = T4_DZ_CUT[key]
        self.T4_DR_CUT = T4_DR_CUT[key]
        self.T4_DTHETA_RZ_CUT = T4_DTHETA_RZ_CUT[key]
        self.T4_CHI2_XY_CUT = T4_CHI2_XY_CUT[key]

        self.prep_t2s()
        self.filter_t2s()
        self.sort_t2s()
        self.make_t4s()


    def prep_t2s(self) -> None:
        logger.info("Preprocessing T2s for T4 making ...")
        self.t2s["ls_doublelayer_div_4"] = self.t2s["ls_doublelayer"] // 4
        self.t2s["ls_doublelayer_mod_4"] = self.t2s["ls_doublelayer"] % 4


    def filter_t2s(self):
        # only consider "good" t2s
        logger.info("Filtering T2s for T4s ...")
        self.t2s = self.t2s[ self.t2s["ls_ok"] ]
        memory = self.t2s.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Memory usage after filtering T2s: {memory:.1f} MB")


    def sort_t2s(self):
        # sort T2s intuitively
        logger.info("Sorting T2s ...")
        cols = [
            "file",
            "i_event",
            "ls_system",
            "ls_doublelayer",
        ]
        self.t2s = self.t2s.sort_values(by=cols).reset_index(drop=True)


    def make_t4s(self) -> None:
        logger.info("Making T4s ...")

        # for each system,
        # layers 0123 combined with layers 4567
        # equivalently, doublelayers 01 combined with 23 (doublelayer // 4 == 0)

        groupby_cols = [
            "ls_system",
            "ls_doublelayer_div_4",
        ] if self.signal else [
            "file",
            "i_event",
            "ls_system",
            "ls_doublelayer_div_4",
        ]

        if self.cut_t4s:
            subgroup_cols = [
                "ls_eta_slice",
                "ls_phi_slice",
            ] if self.signal else [
                "file",
                "i_event",
                "ls_eta_slice",
                "ls_phi_slice",
            ]
        else:
            subgroup_cols = [
                "file",
            ]

        # how to merge lower and upper T2s into T4s
        merge_keys = [
            "file",
            "i_event",
            "ls_system",
            "ls_doublelayer_div_4",
        ]

        def make_t4s_from_group(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:

            group_t4s, group_cutflows = [], []

            entire_lower = df[ df["ls_doublelayer_mod_4"] == 0 ]
            entire_upper = df[ df["ls_doublelayer_mod_4"] != 0 ]
            subgroups = entire_lower.groupby(subgroup_cols)
            n_subgroup = len(subgroups)

            for i_subgroup, (subcols, lower) in enumerate(subgroups):

                # if not cutting, then take all upper
                # else, only consider upper in the same (or neighbor) eta/phi slice as lower
                if not self.cut_t4s:
                    ok = np.ones(len(entire_upper), dtype=bool)
                else:
                    if self.signal:
                        [eta_slice, phi_slice] = subcols
                    else:
                        [_, _, eta_slice, phi_slice] = subcols

                    # get upper data for the same phi slice
                    eta_ok = np.array([eta_slice-1, eta_slice, eta_slice+1]).astype(np.int16)
                    phi_ok = np.array([phi_slice-1, phi_slice, phi_slice+1]).astype(np.int16)
                    phi_ok = phi_ok % N_T4_PHI_SLICES
                    ok = (
                        entire_upper["ls_eta_slice"].isin(eta_ok) &
                        entire_upper["ls_phi_slice"].isin(phi_ok)
                    )
                upper = entire_upper[ok]

                # get all combinations of lower and upper
                t4s = lower.merge(
                    upper,
                    on=merge_keys,
                    how="inner",
                    suffixes=("_lower", "_upper"),
                )
                logger.info(f"{i_subgroup+1} / {n_subgroup}. Lower: {len(lower)}, Upper: {len(upper)}, Combos: {len(t4s)}")

                # the doublelayer
                t4s["t4_doublelayer"] = t4s["ls_doublelayer_lower"]

                # rz projection
                slope_rz = np.divide(t4s["ls_z_upper"] - t4s["ls_z_lower"],
                                     t4s["ls_r_upper"] - t4s["ls_r_lower"])
                t4s["t4_dz"] = t4s["ls_z_lower"] - t4s["ls_r_lower"] * slope_rz

                # xy projection
                slope_xy = np.divide(t4s["ls_y_upper"] - t4s["ls_y_lower"],
                                     t4s["ls_x_upper"] - t4s["ls_x_lower"])
                intercept_xy = t4s["ls_y_lower"] - t4s["ls_x_lower"] * slope_xy
                t4s["t4_dr"] = np.abs(intercept_xy) / np.sqrt(1 + slope_xy**2)

                # assign truth info
                mcp_ok = t4s["i_mcp_lower"] == t4s["i_mcp_upper"]
                t4s["i_mcp"] = t4s["i_mcp_lower"].where(mcp_ok, NO_MCP)
                if self.signal:
                    t4s["t4_first_exit"] = t4s["ls_first_exit_lower"] & t4s["ls_first_exit_upper"]
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
                        t4s[attr] = t4s[f"{attr}_lower"].where(mcp_ok, 0)

                # pass-through the simhit positions
                for coord in ["x", "y", "r"]:
                    t4s[f"t4_{coord}_0"] = t4s[f"ls_{coord}_0_lower"]
                    t4s[f"t4_{coord}_1"] = t4s[f"ls_{coord}_1_lower"]
                    t4s[f"t4_{coord}_2"] = t4s[f"ls_{coord}_2_lower"]
                    t4s[f"t4_{coord}_3"] = t4s[f"ls_{coord}_3_lower"]
                    t4s[f"t4_{coord}_4"] = t4s[f"ls_{coord}_0_upper"]
                    t4s[f"t4_{coord}_5"] = t4s[f"ls_{coord}_1_upper"]
                    t4s[f"t4_{coord}_6"] = t4s[f"ls_{coord}_2_upper"]
                    t4s[f"t4_{coord}_7"] = t4s[f"ls_{coord}_3_upper"]

                # more features
                t4s["t4_deta"] = t4s["ls_eta_upper"] - t4s["ls_eta_lower"]
                t4s["t4_dphi"] = t4s["ls_phi_upper"] - t4s["ls_phi_lower"]
                t4s["t4_dphi"] = (t4s["t4_dphi"] + np.pi) % (2 * np.pi) - np.pi

                # angle differences (handle wraparound)
                t4s["t4_dtheta_rz"] = t4s["ls_theta_rz_upper"] - t4s["ls_theta_rz_lower"]
                t4s["t4_dtheta_rz"] = (t4s["t4_dtheta_rz"] + np.pi) % (2 * np.pi) - np.pi

                # find the circle (radius, x_center, y_center) formed from the first three hits
                BAD_CHI2 = 1e6
                i0, i1, i2 = 0, 4, 7
                ixs = [1, 2, 3, 5, 6]
                circle_d = 2 * (t4s[f"t4_x_{i0}"] * (t4s[f"t4_y_{i1}"] - t4s[f"t4_y_{i2}"]) +
                                t4s[f"t4_x_{i1}"] * (t4s[f"t4_y_{i2}"] - t4s[f"t4_y_{i0}"]) +
                                t4s[f"t4_x_{i2}"] * (t4s[f"t4_y_{i0}"] - t4s[f"t4_y_{i1}"]))
                circle_x = np.divide(t4s[f"t4_r_{i0}"]**2 * (t4s[f"t4_y_{i1}"] - t4s[f"t4_y_{i2}"]) +
                                     t4s[f"t4_r_{i1}"]**2 * (t4s[f"t4_y_{i2}"] - t4s[f"t4_y_{i0}"]) +
                                     t4s[f"t4_r_{i2}"]**2 * (t4s[f"t4_y_{i0}"] - t4s[f"t4_y_{i1}"]),
                                     circle_d)
                circle_y = np.divide(t4s[f"t4_r_{i0}"]**2 * (t4s[f"t4_x_{i2}"] - t4s[f"t4_x_{i1}"]) +
                                     t4s[f"t4_r_{i1}"]**2 * (t4s[f"t4_x_{i0}"] - t4s[f"t4_x_{i2}"]) +
                                     t4s[f"t4_r_{i2}"]**2 * (t4s[f"t4_x_{i1}"] - t4s[f"t4_x_{i0}"]),
                                     circle_d)
                circle_r = np.sqrt((t4s[f"t4_x_{i0}"] - circle_x)**2 + (t4s[f"t4_y_{i0}"] - circle_y)**2)
                circle_ok = circle_d != 0
                if np.any(~circle_ok):
                    logger.warning(f"Found {np.sum(~circle_ok)} invalid circles with circle_d = 0")
                for ix in ixs:
                    circle_diff = np.sqrt((t4s[f"t4_x_{ix}"] - circle_x)**2 + (t4s[f"t4_y_{ix}"] - circle_y)**2) - circle_r
                    t4s[f"t4_chi2_{i0}{i1}{i2}_vs_{ix}"] = np.where(circle_ok, circle_diff**2, BAD_CHI2)

                # calculate the average diff
                chi2cols = [f"t4_chi2_{i0}{i1}{i2}_vs_{ix}" for ix in ixs]
                t4s[f"t4_chi2_{i0}{i1}{i2}"] = t4s[chi2cols].sum(axis=1)

                # rename some things
                rename = {
                    "ls_system": "t4_system",
                    "ls_doublelayer_lower": "t4_doublelayer_lower",
                    "ls_doublelayer_upper": "t4_doublelayer_upper",
                    "ls_module_lower": "t4_module_lower",
                    "ls_module_upper": "t4_module_upper",
                    "ls_sensor_lower": "t4_sensor_lower",
                    "ls_sensor_upper": "t4_sensor_upper",
                    "ls_dr_lower": "t4_dr_lower",
                    "ls_dr_upper": "t4_dr_upper",
                    "ls_dz_lower": "t4_dz_lower",
                    "ls_dz_upper": "t4_dz_upper",
                    "ls_ok_lower": "t4_ls_ok_lower",
                    "ls_ok_upper": "t4_ls_ok_upper",
                }
                t4s = t4s.rename(columns=rename)

                # drop other cols
                dropcols = ["i_mcp_lower", "i_mcp_upper"]
                dropcols.extend([col for col in t4s.columns if col.startswith("simhit_")])
                dropcols.extend([col for col in t4s.columns if col.startswith("doublet_")])
                dropcols.extend([col for col in t4s.columns if col.startswith("ls_")])
                dropcols.extend([col for col in t4s.columns if col.startswith("t2_")])
                dropcols.extend([col for col in t4s.columns if col.startswith("mcp_") and col.endswith("_lower")])
                dropcols.extend([col for col in t4s.columns if col.startswith("mcp_") and col.endswith("_upper")])
                t4s.drop(columns=dropcols, errors="ignore", inplace=True)

                # record some numbers
                cutflow = {"all": len(t4s)}

                # record some cut results
                dl = t4s["t4_doublelayer"]
                t4s["t4_ok_dphi"] = np.abs(t4s["t4_dphi"]) < np.pi / 2.0
                t4s["t4_ok_dz"] = np.abs(t4s["t4_dz"]) < self.T4_DZ_CUT[dl]
                t4s["t4_ok_dr"] = np.abs(t4s["t4_dr"]) < self.T4_DR_CUT[dl]
                t4s["t4_ok_dthetarz"] = np.abs(t4s["t4_dtheta_rz"]) < self.T4_DTHETA_RZ_CUT[dl]
                t4s["t4_ok_chi2xy"] = np.abs(t4s[f"t4_chi2_{i0}{i1}{i2}"]) < self.T4_CHI2_XY_CUT[dl]
                t4s["t4_ok"] = (
                    t4s["t4_ok_dphi"] &
                    t4s["t4_ok_dz"] &
                    t4s["t4_ok_dr"] &
                    t4s["t4_ok_dthetarz"] &
                    t4s["t4_ok_chi2xy"] &
                    np.ones(len(t4s), dtype=bool)
                )

                # remove as desired
                if self.cut_t4s:
                    for cut in [col for col in t4s.columns if col.startswith("t4_ok")]:
                        cutflow[cut] = np.sum(t4s[cut])
                    t4s = t4s[t4s["t4_ok"]]

                # save
                if len(t4s) > 0:
                    group_t4s.append(t4s)
                group_cutflows.append(cutflow)

            return group_t4s, group_cutflows

        # groupby
        groups = self.t2s.groupby(groupby_cols)
        n_group = len(groups)
        all_t4s, all_cutflows = [], []

        # evaluate
        for i_group, (cols, df) in enumerate(groups):
            logger.info(f"Processing group {i_group+1} / {n_group} for T4s (n={len(df)}) ...")
            t4s, cutflow = make_t4s_from_group(df)
            all_t4s.extend(t4s)
            all_cutflows.extend(cutflow)

        # merge cutflow
        cutflow = pd.DataFrame(all_cutflows)
        for col in cutflow.columns:
            logger.info(f"T4s cutflow, {col}: {cutflow[col].sum()}")

        # merge dataframes
        logger.info(f"Merging {len(all_t4s)} groups of T4s ...")
        if len(all_t4s) > 0:
            self.df = pd.concat(all_t4s, ignore_index=True)
        else:
            self.df = pd.DataFrame()
            return

        # sort them
        sortby = [
            "file",
            "i_event",
            "i_mcp",
            "t4_doublelayer",
        ]
        self.df = self.df.sort_values(by=sortby).reset_index(drop=True)

        # announce memory
        memory = self.df.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Memory usage of T4s: {memory:.1f} MB")

