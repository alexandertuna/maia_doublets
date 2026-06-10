import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)

from constants import LS_DZ_CUT, LS_DR_CUT
from constants import LS_DTHETA_RZ_CUT, LS_DTHETA_XY_CUT, LS_CHI2_XY_CUT
from constants import BYTE_TO_MB, NO_MCP
from constants import N_LS_PHI_SLICES
from constants import DETECTOR_MAX_PHI, DETECTOR_MAX_ETA
from constants import N_T4_PHI_SLICES, N_T4_ETA_SLICES

class LineSegment:

    #
    # To make doublets, we do 2 groupbys:
    #  Layers 01, 23, ... grouped by doublet_doublelayer_mod_2
    #  Layers 12, 34, ... grouped by doublet_doublelayer_plus_1_mod_2
    #

    def __init__(self, geometry_version: str, sim: bool, smear: str, doublets: pd.DataFrame, signal: bool, cut_line_segments: bool):
        self.df = None
        self.signal = signal
        self.cut_line_segments = cut_line_segments
        self.lower_suffix = "lower"
        self.upper_suffix = "upper"
        memory = doublets.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Making linesegments with doublets memory {memory:.1f} MB ...")

        key = (geometry_version, "sim") if sim else (geometry_version, "digi", smear)
        self.LS_DZ_CUT = LS_DZ_CUT[key]
        self.LS_DR_CUT = LS_DR_CUT[key]
        self.LS_DTHETA_RZ_CUT = LS_DTHETA_RZ_CUT[key]
        self.LS_DTHETA_XY_CUT = LS_DTHETA_XY_CUT[key]
        self.LS_CHI2_XY_CUT = LS_CHI2_XY_CUT[key]

        self.doublets = doublets.copy()
        self.prep_doublets()
        self.filter_doublets()
        self.sort_doublets()
        self.make_linesegments()


    def prep_doublets(self):
        # add columns for 2 groupbys
        self.doublets["doublet_doublelayer_div_2"] = self.doublets["doublet_doublelayer"] // 2
        self.doublets["doublet_doublelayer_mod_2"] = self.doublets["doublet_doublelayer"] % 2
        self.doublets["doublet_doublelayer_plus_1"] = self.doublets["doublet_doublelayer"] + 1
        self.doublets["doublet_doublelayer_plus_1_div_2"] = self.doublets["doublet_doublelayer_plus_1"] // 2
        self.doublets["doublet_doublelayer_plus_1_mod_2"] = self.doublets["doublet_doublelayer_plus_1"] % 2

        # announce memory
        memory = self.doublets.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Memory usage after adding/removing columns: {memory:.1f} MB")


    def filter_doublets(self):
        # only consider "good" doublets
        logger.info("Filtering doublets for line segments ...")
        self.doublets = self.doublets[ self.doublets["doublet_ok"] ]
        memory = self.doublets.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Memory usage after filtering doublets: {memory:.1f} MB")


    def sort_doublets(self):
        # sort doublets by file, event, system, doublelayer, sensor, module
        logger.info("Sorting doublets for line segments ...")
        cols = [
            "file",
            "i_event",
            "doublet_system",
            "doublet_doublelayer",
            "doublet_sensor",
            "doublet_module",
        ]
        self.doublets = self.doublets.sort_values(by=cols).reset_index(drop=True)


    def make_linesegments(self):
        logger.info("Making line segments ...")

        # two loops:
        # layers 01, 23, 45, ... (even)
        # layers 12, 34, 56, ... (odd)
        even, odd = 0, 1

        # groupby scheme
        groupby_cols = {
            even: [
                "doublet_system",
                "doublet_doublelayer_div_2",
            ] if self.signal else [
                "file",
                "i_event",
                "doublet_system",
                "doublet_doublelayer_div_2",
            ],
            odd: [
                "doublet_system",
                "doublet_doublelayer_plus_1_div_2",
            ] if self.signal else [
                "file",
                "i_event",
                "doublet_system",
                "doublet_doublelayer_plus_1_div_2",
            ],
        }

        # sub-groupby scheme
        if self.cut_line_segments:
            subgroup_cols = [
                "doublet_eta_slice",
                "doublet_phi_slice",
            ] if self.signal else [
                "file",
                "i_event",
                "doublet_eta_slice",
                "doublet_phi_slice",
            ]
        else:
            subgroup_cols = [
                "file",
            ]

        # how to merge lower MD and upper MD
        keys = {
            even: [
                "file",
                "i_event",
                "doublet_system",
                "doublet_doublelayer_div_2",
            ],
            odd: [
                "file",
                "i_event",
                "doublet_system",
                "doublet_doublelayer_plus_1_div_2",
            ],
        }

        lower_vs_upper = {
            even: "doublet_doublelayer_mod_2",
            odd: "doublet_doublelayer_plus_1_mod_2",
        }

        all_cutflows = []
        all_linesegments = []

        for start in [
            even,
            # odd,
        ]:

            for i_group, (cols, df) in enumerate(self.doublets.groupby(groupby_cols[start])):

                group_cutflows = []

                # progress bar
                if (self.signal and i_group % 10 == 0) or (not self.signal):
                    n_group = len(self.doublets.groupby(groupby_cols[start]))
                    logger.info(f"Processing group {i_group+1} / {n_group} for line segments (n={len(df)}) ...")

                # get lower doublets and upper doublets
                entire_lower = df[ df[lower_vs_upper[start]] == 0 ]
                entire_upper = df[ df[lower_vs_upper[start]] != 0 ]

                # get lower data for each eta,phi slice
                subgroups = entire_lower.groupby(subgroup_cols)
                n_subgroup = len(subgroups)
                for i_subgroup, (subcols, lower) in enumerate(subgroups):

                    # if not cutting line segments, then take all upper doublets for this group
                    if not self.cut_line_segments:
                        ok = np.ones(len(entire_upper), dtype=bool)

                    # else, only consider upper in the same (or neighbor) eta/phi slice as lower
                    else:
                        if self.signal:
                            [eta_slice, phi_slice] = subcols
                        else:
                            [_, _, eta_slice, phi_slice] = subcols

                        # get upper data for the same phi slice
                        eta_ok = np.array([eta_slice-1, eta_slice, eta_slice+1]).astype(np.int16)
                        phi_ok = np.array([phi_slice-1, phi_slice, phi_slice+1]).astype(np.int16)
                        phi_ok = phi_ok % N_LS_PHI_SLICES
                        ok = (
                            entire_upper["doublet_eta_slice"].isin(eta_ok) &
                            entire_upper["doublet_phi_slice"].isin(phi_ok)
                        )

                    upper = entire_upper[ok]

                    # get all combinations of lower and upper
                    segments = lower.merge(
                        upper,
                        on=keys[start],
                        how="inner",
                        suffixes=("_lower", "_upper"),
                    )

                    # the doublelayer
                    segments["ls_doublelayer"] = segments["doublet_doublelayer_lower"]

                    # rz projection
                    slope_rz = np.divide(segments["doublet_z_upper"] - segments["doublet_z_lower"],
                                         segments["doublet_r_upper"] - segments["doublet_r_lower"])
                    segments["ls_dz"] = segments["doublet_z_lower"] - segments["doublet_r_lower"] * slope_rz
                    segments["ls_theta_rz"] = np.arctan(slope_rz)

                    # xy projection
                    slope_xy = np.divide(segments["doublet_y_upper"] - segments["doublet_y_lower"],
                                         segments["doublet_x_upper"] - segments["doublet_x_lower"])
                    intercept_xy = segments["doublet_y_lower"] - segments["doublet_x_lower"] * slope_xy
                    segments["ls_dr"] = np.abs(intercept_xy) / np.sqrt(1 + slope_xy**2)

                    # cut some segments? do this early to save computations
                    if self.cut_line_segments:
                        dl = segments["ls_doublelayer"]
                        segments["ls_ok_dz"] = np.abs(segments["ls_dz"]) < self.LS_DZ_CUT[dl]
                        segments["ls_ok_dr"] = np.abs(segments["ls_dr"]) < self.LS_DR_CUT[dl]
                        segments = segments[segments["ls_ok_dz"] & segments["ls_ok_dr"]]

                    # assign truth info
                    mcp_ok = segments["i_mcp_lower"] == segments["i_mcp_upper"]
                    segments["i_mcp"] = segments["i_mcp_lower"].where(mcp_ok, NO_MCP)
                    if self.signal:
                        segments["ls_first_exit"] = segments["doublet_first_exit_lower"] & segments["doublet_first_exit_upper"]
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
                            segments[attr] = segments[f"{attr}_lower"].where(mcp_ok, 0)

                    # features: position
                    segments["ls_r"] = (segments["doublet_r_lower"] + segments["doublet_r_upper"]) / 2
                    segments["ls_z"] = (segments["doublet_z_lower"] + segments["doublet_z_upper"]) / 2
                    segments["ls_x"] = (segments["doublet_x_lower"] + segments["doublet_x_upper"]) / 2
                    segments["ls_y"] = (segments["doublet_y_lower"] + segments["doublet_y_upper"]) / 2
                    segments["ls_phi"] = np.arctan2(segments["ls_y"], segments["ls_x"])
                    segments["ls_theta"] = np.arctan2(segments["ls_r"], segments["ls_z"])
                    segments["ls_eta"] = -np.log(np.tan(segments["ls_theta"] / 2))
                    segments["ls_phi_slice"] = np.floor((segments["ls_phi"] + DETECTOR_MAX_PHI) / (2 * DETECTOR_MAX_PHI) * N_T4_PHI_SLICES).astype(np.int16)
                    segments["ls_eta_slice"] = np.floor((segments["ls_eta"] + DETECTOR_MAX_ETA) / (2 * DETECTOR_MAX_ETA) * N_T4_ETA_SLICES).astype(np.int16)

                    # assign more features
                    segments["ls_ddr"] = segments["doublet_dr_upper"] - segments["doublet_dr_lower"]
                    segments["ls_ddz"] = segments["doublet_dz_upper"] - segments["doublet_dz_lower"]
                    segments["ls_deta"] = segments["doublet_eta_upper"] - segments["doublet_eta_lower"]
                    segments["ls_dphi"] = segments["doublet_phi_upper"] - segments["doublet_phi_lower"]
                    segments["ls_dphi"] = (segments["ls_dphi"] + np.pi) % (2 * np.pi) - np.pi
                    segments["ls_dqoverpt"] = segments["doublet_qoverpt_upper"] - segments["doublet_qoverpt_lower"]
                    segments["ls_doublelayer_div_4"] = segments["doublet_doublelayer_lower"] // 4
                    segments["ls_doublelayer_mod_4"] = segments["doublet_doublelayer_lower"] % 4
                    segments["ls_doublelayer_even"] = (segments["doublet_doublelayer_lower"] % 2 == 0).astype(bool)

                    # pass-through the simhit positions
                    for coord in ["x", "y", "r"]:
                        segments[f"ls_{coord}_0"] = segments[f"doublet_{coord}_0_lower"]
                        segments[f"ls_{coord}_1"] = segments[f"doublet_{coord}_1_lower"]
                        segments[f"ls_{coord}_2"] = segments[f"doublet_{coord}_0_upper"]
                        segments[f"ls_{coord}_3"] = segments[f"doublet_{coord}_1_upper"]

                    # angle differences (handle wraparound)
                    segments["ls_dtheta_rz"] = segments["doublet_theta_rz_upper"] - segments["doublet_theta_rz_lower"]
                    segments["ls_dtheta_xy"] = segments["doublet_theta_xy_upper"] - segments["doublet_theta_xy_lower"]
                    segments["ls_dtheta_rz"] = (segments["ls_dtheta_rz"] + np.pi) % (2 * np.pi) - np.pi
                    segments["ls_dtheta_xy"] = (segments["ls_dtheta_xy"] + np.pi) % (2 * np.pi) - np.pi

                    # find the circle (radius, x_center, y_center) formed from the first three hits
                    BAD_CHI2 = 1e6
                    circle_d = 2 * (segments["ls_x_0"] * (segments["ls_y_1"] - segments["ls_y_2"]) +
                                    segments["ls_x_1"] * (segments["ls_y_2"] - segments["ls_y_0"]) +
                                    segments["ls_x_2"] * (segments["ls_y_0"] - segments["ls_y_1"]))
                    circle_x = np.divide(segments["ls_r_0"]**2 * (segments["ls_y_1"] - segments["ls_y_2"]) +
                                         segments["ls_r_1"]**2 * (segments["ls_y_2"] - segments["ls_y_0"]) +
                                         segments["ls_r_2"]**2 * (segments["ls_y_0"] - segments["ls_y_1"]),
                                         circle_d)
                    circle_y = np.divide(segments["ls_r_0"]**2 * (segments["ls_x_2"] - segments["ls_x_1"]) +
                                         segments["ls_r_1"]**2 * (segments["ls_x_0"] - segments["ls_x_2"]) +
                                         segments["ls_r_2"]**2 * (segments["ls_x_1"] - segments["ls_x_0"]),
                                         circle_d)
                    circle_r = np.sqrt((segments["ls_x_0"] - circle_x)**2 + (segments["ls_y_0"] - circle_y)**2)
                    circle_ok = circle_d != 0
                    if np.any(~circle_ok):
                        logger.warning(f"Found {np.sum(~circle_ok)} invalid circles with circle_d = 0")
                    circle_diff = np.sqrt((segments["ls_x_3"] - circle_x)**2 + (segments["ls_y_3"] - circle_y)**2) - circle_r

                    # calculate the distance from (x_3, y_3) to the circle
                    segments["ls_chi2_012"] = np.where(circle_ok, circle_diff**2, BAD_CHI2)

                    # rename some things
                    rename = {
                        "doublet_system": "ls_system",
                        "doublet_doublelayer_lower": "ls_doublelayer_lower",
                        "doublet_doublelayer_upper": "ls_doublelayer_upper",
                        "doublet_module_lower": "ls_module_lower",
                        "doublet_module_upper": "ls_module_upper",
                        "doublet_sensor_lower": "ls_sensor_lower",
                        "doublet_sensor_upper": "ls_sensor_upper",
                        "doublet_dr_lower": "ls_dr_lower",
                        "doublet_dr_upper": "ls_dr_upper",
                        "doublet_dz_lower": "ls_dz_lower",
                        "doublet_dz_upper": "ls_dz_upper",
                        "doublet_ok_lower": "ls_md_ok_lower",
                        "doublet_ok_upper": "ls_md_ok_upper",
                    }
                    segments = segments.rename(columns=rename)

                    # assign module and sensor from lower doublet (arbitrary choice)
                    segments["ls_module"] = segments["ls_module_lower"]
                    segments["ls_sensor"] = segments["ls_sensor_lower"]

                    # and drop other cols
                    dropcols = ["i_mcp_lower", "i_mcp_upper"]
                    dropcols.extend([col for col in segments.columns if col.startswith("simhit_")])
                    dropcols.extend([col for col in segments.columns if col.startswith("doublet_")])
                    dropcols.extend([col for col in segments.columns if col.startswith("mcp_") and col.endswith("_lower")])
                    dropcols.extend([col for col in segments.columns if col.startswith("mcp_") and col.endswith("_upper")])
                    segments.drop(columns=dropcols, errors="ignore", inplace=True)

                    # record some numbers
                    cutflow = {"all": len(segments)}

                    # record some cut results
                    dl = segments["ls_doublelayer"]
                    segments["ls_ok_dtheta_rz"] = np.abs(segments["ls_dtheta_rz"]) < self.LS_DTHETA_RZ_CUT[dl]
                    segments["ls_ok_dtheta_xy"] = np.abs(segments["ls_dtheta_xy"]) < self.LS_DTHETA_XY_CUT[dl]
                    segments["ls_ok_dz"] = np.abs(segments["ls_dz"]) < self.LS_DZ_CUT[dl]
                    segments["ls_ok_dr"] = np.abs(segments["ls_dr"]) < self.LS_DR_CUT[dl]
                    segments["ls_ok_dphi"] = np.abs(segments["ls_dphi"]) < np.pi / 2.0
                    segments["ls_ok_chi2_xy"] = np.abs(segments["ls_chi2_012"]) < self.LS_CHI2_XY_CUT[dl]
                    segments["ls_ok_drdz"] = segments["ls_ok_dz"] & segments["ls_ok_dr"] & segments["ls_ok_dphi"]
                    segments["ls_ok_drdzdthetarz"] = segments["ls_ok_dz"] & segments["ls_ok_dr"] & segments["ls_ok_dphi"] & segments["ls_ok_dtheta_rz"]
                    segments["ls_ok"] = (
                        segments["ls_ok_dz"] &
                        segments["ls_ok_dr"] &
                        segments["ls_ok_dphi"] &
                        segments["ls_ok_dtheta_rz"] &
                        segments["ls_ok_chi2_xy"] &
                        np.ones(len(segments), dtype=bool)
                    )

                    # remove as desired
                    if self.cut_line_segments:
                        for cut in [col for col in segments.columns if col.startswith("ls_ok")]:
                            cutflow[cut] = np.sum(segments[cut])
                        segments = segments[segments["ls_ok"]]

                    # progress bar and stats from this group
                    if i_subgroup > 0 and i_subgroup % 100 == 0:
                        logger.info(f"Processed subgroup {i_subgroup} / {n_subgroup} which has {len(segments)} line segments ...")

                    # save them
                    group_cutflows.append(cutflow)
                    all_cutflows.append(cutflow)
                    all_linesegments.append(segments)

                # group cutflow
                cutflow = pd.DataFrame(group_cutflows)
                for col in cutflow.columns:
                    logger.info(f"Line segments cutflow (group), {col}: {cutflow[col].sum()}")


        # merge them
        logger.info(f"Merging {len(all_linesegments)} groups of line segments ...")
        self.df = pd.concat(all_linesegments, ignore_index=True)

        # sort them
        sortby = [
            "file",
            "i_event",
            "i_mcp",
            "ls_doublelayer",
            "ls_module_lower",
            "ls_module_upper",
            "ls_sensor_lower",
            "ls_sensor_upper",
        ]
        self.df = self.df.sort_values(by=sortby).reset_index(drop=True)

        # announce memory
        memory = self.df.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Memory usage of line segments: {memory:.1f} MB")

        # cutflow
        cutflow = pd.DataFrame(all_cutflows)
        for col in cutflow.columns:
            logger.info(f"Line segments cutflow, {col}: {cutflow[col].sum()}")


