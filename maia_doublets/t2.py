import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)

from maia_doublets.constants import T2_DZ_CUT, T2_DR_CUT
from maia_doublets.constants import T2_DTHETA_RZ_CUT, T2_CHI2_XY_CUT
from maia_doublets.constants import BYTE_TO_MB, NO_MCP
from maia_doublets.constants import N_T2_PHI_SLICES
from maia_doublets.constants import DETECTOR_MAX_PHI, DETECTOR_MAX_ETA
from maia_doublets.constants import N_T4_PHI_SLICES, N_T4_ETA_SLICES

class T2Maker:

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
        self.T2_DZ_CUT = T2_DZ_CUT[key]
        self.T2_DR_CUT = T2_DR_CUT[key]
        self.T2_DTHETA_RZ_CUT = T2_DTHETA_RZ_CUT[key]
        self.T2_CHI2_XY_CUT = T2_CHI2_XY_CUT[key]

        self.merge_keys = [
            "file",
            "i_event",
            "doublet_phi_slice",
            "doublet_eta_slice",
        ]

        self.doublets = doublets.copy()
        self.filter_doublets()
        self.sort_doublets()
        self.make_t2s()
        # self.make_linesegments()


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


    def make_t2s(self):

        # split mds into global doublelayer once, up front
        logger.info(f"Splitting MDs ...")
        mds = {
            gdl: group.reset_index(drop=True)
            for gdl, group
            in self.doublets.groupby("doublet_gdoublelayer", sort=True)
        }

        # make T2s from neighboring doublelayers
        gdoublelayer_pairs = [
            (2, 3),
            (4, 5),
            (6, 7),
        ]
        all_t2s, all_cutflows = [], []
        for (lower, upper) in gdoublelayer_pairs:
            logger.info(f"Making T2s from gdl={lower} and gdl={upper} ...")
            if lower not in mds or upper not in mds:
                logger.warning(f"Missing gdl={lower} or gdl={upper}, skipping ...")
                continue
            lower_df = mds[lower]
            upper_df = mds[upper]
            t2s, cutflow = self.make_t2s_from_lower_upper(lower_df, upper_df)
            for col in cutflow:
                logger.info(f"T2s cutflow, gdl={lower}-{upper}, {col}: {cutflow[col]}")
            all_t2s.append(t2s)
            all_cutflows.append(cutflow)

        # merge dataframes
        logger.info(f"Merging {len(all_t2s)} groups of T2s ...")
        self.df = pd.concat(all_t2s, ignore_index=True) if len(all_t2s) > 0 else pd.DataFrame()

        # merge cutflow
        cutflow = pd.DataFrame(all_cutflows)
        for col in cutflow.columns:
            logger.info(f"T2s cutflow, {col}: {cutflow[col].sum()}")

        # announce memory
        memory = self.df.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Memory usage of T2s: {memory:.1f} MB")


    def make_t2s_from_lower_upper(self, lower: pd.DataFrame, upper: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
        """
        Given an lower df and upper df, make T2s from their combinations
            For example, lower df could be MDs which span OT layers 0-1 (global doublelayer 4),
            and upper df could be MDs which span OT layers 2-3 (global doublelayer 5).
        """

        n_t2_phi_slices = N_T2_PHI_SLICES[upper["doublet_system"]]

        # get combinations of lower and upper
        # within neighboring phi/eta slices
        cands = []
        for phi_shift in (-1, 0, 1):
            for eta_shift in (-1, 0, 1):
                shifted = upper.assign(
                    doublet_phi_slice=(upper["doublet_phi_slice"] + phi_shift) % n_t2_phi_slices,
                    doublet_eta_slice=(upper["doublet_eta_slice"] + eta_shift),
                )
                logger.info(f"Merging phi_shift={phi_shift}, eta_shift={eta_shift} ...")
                t2s = lower.merge(shifted, on=self.merge_keys, suffixes=("_lower", "_upper"))
                t2s = self.rename_t2_columns_after_merge(t2s)

                # cut some t2s early to save computations
                if self.cut_line_segments:
                    t2s = self.add_basic_t2_features(t2s)
                    sy = t2s["ls_system"]
                    dl = t2s["ls_doublelayer"]
                    t2s["ls_ok_dz"] = np.abs(t2s["ls_dz"]) < self.T2_DZ_CUT[sy, dl]
                    t2s["ls_ok_dr"] = np.abs(t2s["ls_dr"]) < self.T2_DR_CUT[sy, dl]
                    t2s = t2s[t2s["ls_ok_dz"] & t2s["ls_ok_dr"]]

                cands.append(t2s)

        # combine candidates into one dataframe
        for t2s in cands:
            memory = int(t2s.memory_usage(deep=True).sum() * BYTE_TO_MB)
            logger.info(f"Slice with phi_shift={phi_shift}, eta_shift={eta_shift} has {len(t2s)} entries and {memory} MB ...")
        logger.info(f"Concatenating T2s ...")
        t2s = pd.concat(cands, ignore_index=True)

        # calculate T2 features and cuts
        logger.info(f"Calculating T2 features for {len(t2s)} T2s ...")
        t2s, cutflow = self.consolidate_t2_features(t2s.copy())
        return t2s, cutflow


    def rename_t2_columns_after_merge(self, t2s: pd.DataFrame) -> pd.DataFrame:
        rename = {
            "doublet_doublelayer_lower": "ls_doublelayer_lower",
            "doublet_doublelayer_upper": "ls_doublelayer_upper",
            "doublet_glayer_lower": "ls_glayer_lower",
            "doublet_glayer_upper": "ls_glayer_upper",
            "doublet_module_lower": "ls_module_lower",
            "doublet_module_upper": "ls_module_upper",
            "doublet_sensor_lower": "ls_sensor_lower",
            "doublet_sensor_upper": "ls_sensor_upper",
            "doublet_system_lower": "ls_system_lower",
            "doublet_system_upper": "ls_system_upper",
            "doublet_dr_lower": "ls_dr_lower",
            "doublet_dr_upper": "ls_dr_upper",
            "doublet_dz_lower": "ls_dz_lower",
            "doublet_dz_upper": "ls_dz_upper",
            "doublet_ok_lower": "ls_md_ok_lower",
            "doublet_ok_upper": "ls_md_ok_upper",
        }
        return t2s.rename(columns=rename, copy=False)


    def add_basic_t2_features(self, t2s: pd.DataFrame) -> pd.DataFrame:

        # the system and doublelayer
        t2s["ls_system"] = t2s["ls_system_lower"].astype(np.uint8)
        t2s["ls_doublelayer"] = t2s["ls_doublelayer_lower"].astype(np.uint8)

        # rz projection
        slope_rz = np.divide(t2s["doublet_z_upper"] - t2s["doublet_z_lower"],
                             t2s["doublet_r_upper"] - t2s["doublet_r_lower"])
        t2s["ls_dz"] = t2s["doublet_z_lower"] - t2s["doublet_r_lower"] * slope_rz
        t2s["ls_theta_rz"] = np.arctan(slope_rz)

        # xy projection
        slope_xy = np.divide(t2s["doublet_y_upper"] - t2s["doublet_y_lower"],
                             t2s["doublet_x_upper"] - t2s["doublet_x_lower"])
        intercept_xy = t2s["doublet_y_lower"] - t2s["doublet_x_lower"] * slope_xy
        t2s["ls_dr"] = np.abs(intercept_xy) / np.sqrt(1 + slope_xy**2)

        return t2s


    def consolidate_t2_features(self, segments: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:

        segments = self.add_basic_t2_features(segments)

        # assign truth info
        mcp_ok = segments["i_mcp_lower"] == segments["i_mcp_upper"]
        segments["i_mcp"] = segments["i_mcp_lower"].where(mcp_ok, NO_MCP)
        if self.signal:
            segments["ls_first_exit"] = segments["doublet_first_exit_lower"] & segments["doublet_first_exit_upper"]
            segments["ls_from_fiducial_mcp"] = segments["doublet_from_fiducial_mcp_lower"] & segments["doublet_from_fiducial_mcp_upper"]
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
        segments["ls_ddr"] = segments["ls_dr_upper"] - segments["ls_dr_lower"]
        segments["ls_ddz"] = segments["ls_dz_upper"] - segments["ls_dz_lower"]
        segments["ls_deta"] = segments["doublet_eta_upper"] - segments["doublet_eta_lower"]
        segments["ls_dphi"] = segments["doublet_phi_upper"] - segments["doublet_phi_lower"]
        segments["ls_dphi"] = (segments["ls_dphi"] + np.pi) % (2 * np.pi) - np.pi
        segments["ls_dqoverpt"] = segments["doublet_qoverpt_upper"] - segments["doublet_qoverpt_lower"]

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

        # assign features from lower doublet (arbitrary choice)
        segments["ls_module"] = segments["ls_module_lower"]
        segments["ls_sensor"] = segments["ls_sensor_lower"]
        segments["ls_glayer"] = segments["ls_glayer_lower"]
        segments["ls_gdoublelayer"] = segments["ls_glayer"] // 2

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
        sy = segments["ls_system"]
        dl = segments["ls_doublelayer"]
        segments["ls_ok_dtheta_rz"] = np.abs(segments["ls_dtheta_rz"]) < self.T2_DTHETA_RZ_CUT[sy, dl]
        segments["ls_ok_dz"] = np.abs(segments["ls_dz"]) < self.T2_DZ_CUT[sy, dl]
        segments["ls_ok_dr"] = np.abs(segments["ls_dr"]) < self.T2_DR_CUT[sy, dl]
        segments["ls_ok_dphi"] = np.abs(segments["ls_dphi"]) < np.pi / 2.0
        segments["ls_ok_chi2_xy"] = np.abs(segments["ls_chi2_012"]) < self.T2_CHI2_XY_CUT[sy, dl]
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
        for cut in [col for col in segments.columns if col.startswith("ls_ok")]:
            cutflow[cut] = np.sum(segments[cut])
        if self.cut_line_segments:
            segments = segments[segments["ls_ok"]]

        # fin
        return segments, cutflow


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
                group_system = df["doublet_system"].iloc[0]

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
                        phi_ok = phi_ok % N_T2_PHI_SLICES[group_system]
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

                    # the system and doublelayer
                    segments["ls_system"] = segments["doublet_system"]
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
                        sy = segments["ls_system"]
                        dl = segments["ls_doublelayer"]
                        segments["ls_ok_dz"] = np.abs(segments["ls_dz"]) < self.T2_DZ_CUT[sy, dl]
                        segments["ls_ok_dr"] = np.abs(segments["ls_dr"]) < self.T2_DR_CUT[sy, dl]
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
                        "doublet_doublelayer_lower": "ls_doublelayer_lower",
                        "doublet_doublelayer_upper": "ls_doublelayer_upper",
                        "doublet_glayer_lower": "ls_glayer_lower",
                        "doublet_glayer_upper": "ls_glayer_upper",
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

                    # assign features from lower doublet (arbitrary choice)
                    segments["ls_module"] = segments["ls_module_lower"]
                    segments["ls_sensor"] = segments["ls_sensor_lower"]
                    segments["ls_glayer"] = segments["ls_glayer_lower"]
                    segments["ls_gdoublelayer"] = segments["ls_glayer"] // 2

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
                    sy = segments["ls_system"]
                    dl = segments["ls_doublelayer"]
                    segments["ls_ok_dtheta_rz"] = np.abs(segments["ls_dtheta_rz"]) < self.T2_DTHETA_RZ_CUT[sy, dl]
                    segments["ls_ok_dz"] = np.abs(segments["ls_dz"]) < self.T2_DZ_CUT[sy, dl]
                    segments["ls_ok_dr"] = np.abs(segments["ls_dr"]) < self.T2_DR_CUT[sy, dl]
                    segments["ls_ok_dphi"] = np.abs(segments["ls_dphi"]) < np.pi / 2.0
                    segments["ls_ok_chi2_xy"] = np.abs(segments["ls_chi2_012"]) < self.T2_CHI2_XY_CUT[sy, dl]
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
                        logger.info(f"Processed subgroup {i_subgroup} / {n_subgroup} which has {len(segments)} T2s ...")

                    # save them
                    group_cutflows.append(cutflow)
                    all_cutflows.append(cutflow)
                    all_linesegments.append(segments)

                # group cutflow
                cutflow = pd.DataFrame(group_cutflows)
                for col in cutflow.columns:
                    logger.info(f"T2s cutflow (group), {col}: {cutflow[col].sum()}")


        # merge them
        logger.info(f"Merging {len(all_linesegments)} groups of T2s ...")
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
        logger.info(f"Memory usage of T2s: {memory:.1f} MB")

        # cutflow
        cutflow = pd.DataFrame(all_cutflows)
        for col in cutflow.columns:
            logger.info(f"T2s cutflow, {col}: {cutflow[col].sum()}")
        if not self.cut_line_segments:
            col = "ls_ok"
            logger.info(f"T2s cutflow, {col}: {self.df[col].sum()}")


