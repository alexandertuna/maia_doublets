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

from maia_doublets.constants import BYTE_TO_MB, NO_MCP
from maia_doublets.constants import N_T4_PHI_SLICES
from maia_doublets.constants import N_LAYERS_IN_T4
from maia_doublets.constants import N_T8_PHI_SLICES, N_T8_ETA_SLICES, DETECTOR_MAX_ETA

class T4Maker:

    def __init__(
        self,
        signal: bool,
        cut_t4s: bool,
        calibs: dict,
        t2s: pd.DataFrame,
    ):
        self.df = None
        self.signal = signal
        self.cut_t4s = cut_t4s
        self.t2s = t2s.copy()
        memory = self.t2s.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Making T4s. T2 dataframe size: {memory:.2f} MB")

        self.T4_DZ_CUT = calibs.get("t4_dz", np.zeros((10, 10)))
        self.T4_DR_CUT = calibs.get("t4_dr", np.zeros((10, 10)))
        self.T4_DTHETA_RZ_CUT = calibs.get("t4_dtheta_rz", np.zeros((10, 10)))
        self.T4_CHI2_XY_CUT = calibs.get("t4_chi2_xy_047", np.zeros((10, 10)))
        self.T4_CHI2_SZ_CUT = calibs.get("t4_chi2_sz", np.zeros((10, 10)))

        # how to merge lower and upper T2s into T4s
        self.merge_keys = [
            "file",
            "i_event",
            "ls_phi_slice",
            "ls_eta_slice",
        ]

        self.filter_t2s()
        self.sort_t2s()
        self.make_t4s()


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
        """
        Todo: add description of this T4-making algorithm
        """

        # split t2s into global doublelayer once, up front
        logger.info(f"Splitting T2s ...")
        t2s = {
            gdl: group.reset_index(drop=True)
            for gdl, group
            in self.t2s.groupby("ls_gdoublelayer", sort=True)
        }

        # make T4s from neighboring doublelayers
        gdoublelayer_pairs = [
            (0, 2), # IT
            (4, 6), # OT
        ]
        all_t4s, all_cutflows = [], []
        for (lower, upper) in gdoublelayer_pairs:
            logger.info(f"Making T4s from gdl={lower} and gdl={upper} ...")
            lower_df = t2s[lower]
            upper_df = t2s[upper]
            t4s, cutflow = self.make_t4s_from_lower_upper(lower_df, upper_df)
            for col in cutflow:
                logger.info(f"T4s cutflow, gdl={lower}-{upper}, {col}: {cutflow[col]}")
            all_t4s.append(t4s)
            all_cutflows.append(cutflow)

        # merge dataframes
        logger.info(f"Merging {len(all_t4s)} groups of T4s ...")
        self.df = pd.concat(all_t4s, ignore_index=True)

        # merge cutflow
        cutflow = pd.DataFrame(all_cutflows)
        for col in cutflow.columns:
            logger.info(f"T4s cutflow, {col}: {cutflow[col].sum()}")

        # announce memory
        memory = self.df.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Memory usage of T4s: {memory:.1f} MB")



    def make_t4s_from_lower_upper(self, lower: pd.DataFrame, upper: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
        """
        Given an lower df and upper df, make T4s from their combinations
            For example, lower df could be T2s which span OT layers 0-3 (global doublelayer 4),
            and upper df could be T2s which span OT layers 4-7 (global doublelayer 6).
        """
        # get combinations of lower and upper
        # within neighboring phi/eta slices
        logger.info(f"Total number of lower T2s: {len(lower)}")
        logger.info(f"Total number of upper T2s: {len(upper)}")
        cands = []
        for phi_shift in (-1, 0, 1):
            for eta_shift in (-1, 0, 1):
                shifted = upper.assign(
                    ls_phi_slice=(upper["ls_phi_slice"] + phi_shift) % N_T4_PHI_SLICES,
                    ls_eta_slice=(upper["ls_eta_slice"] + eta_shift),
                )
                cands.append(lower.merge(shifted, on=self.merge_keys, suffixes=("_lower", "_upper")))

        # combine candidates into one dataframe
        t4s = pd.concat(cands, ignore_index=True)

        # calculate T4 features and cuts
        t4s, cutflow = self.consolidate_t4_features(t4s)
        return t4s, cutflow


    def consolidate_t4_features(self, t4s: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:

        # the doublelayer
        t4s["t4_doublelayer"] = t4s["ls_doublelayer_lower"]
        t4s["t4_gdoublelayer"] = t4s["ls_gdoublelayer_lower"]

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
            t4s["t4_from_fiducial_mcp"] = t4s["ls_from_fiducial_mcp_lower"] & t4s["ls_from_fiducial_mcp_upper"]
            t4s["t4_detectable"] = mcp_ok & t4s["ls_detectable_lower"] & t4s["ls_detectable_upper"]
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

        # collect new columns
        new = {}

        # pass-through the simhit positions
        for coord in ["x", "y", "r", "z"]:
            new[f"t4_{coord}_0"] = t4s[f"ls_{coord}_0_lower"]
            new[f"t4_{coord}_1"] = t4s[f"ls_{coord}_1_lower"]
            new[f"t4_{coord}_2"] = t4s[f"ls_{coord}_2_lower"]
            new[f"t4_{coord}_3"] = t4s[f"ls_{coord}_3_lower"]
            new[f"t4_{coord}_4"] = t4s[f"ls_{coord}_0_upper"]
            new[f"t4_{coord}_5"] = t4s[f"ls_{coord}_1_upper"]
            new[f"t4_{coord}_6"] = t4s[f"ls_{coord}_2_upper"]
            new[f"t4_{coord}_7"] = t4s[f"ls_{coord}_3_upper"]

        # more features
        new["t4_deta"] = t4s["ls_eta_upper"] - t4s["ls_eta_lower"]
        new["t4_dphi"] = t4s["ls_phi_upper"] - t4s["ls_phi_lower"]
        new["t4_dphi"] = (new["t4_dphi"] + np.pi) % (2 * np.pi) - np.pi
        new["t4_eta"] = 0.5 * (t4s["ls_eta_upper"] + t4s["ls_eta_lower"])
        new["t4_x"] = 0.5 * (t4s["ls_x_upper"] + t4s["ls_x_lower"])
        new["t4_y"] = 0.5 * (t4s["ls_y_upper"] + t4s["ls_y_lower"])
        new["t4_z"] = 0.5 * (t4s["ls_z_upper"] + t4s["ls_z_lower"])
        new["t4_r"] = 0.5 * (t4s["ls_r_upper"] + t4s["ls_r_lower"])
        new["t4_phi"] = np.arctan2(new["t4_y"], new["t4_x"])
        new["t4_phi_slice"] = np.floor((new["t4_phi"] + np.pi) / (2 * np.pi) * N_T8_PHI_SLICES).astype(np.int16)
        new["t4_eta_slice"] = np.floor((new["t4_eta"] + DETECTOR_MAX_ETA) / (2 * DETECTOR_MAX_ETA) * N_T8_ETA_SLICES).astype(np.int16)

        # angle differences (handle wraparound)
        new["t4_dtheta_rz"] = t4s["ls_theta_rz_upper"] - t4s["ls_theta_rz_lower"]
        new["t4_dtheta_rz"] = (new["t4_dtheta_rz"] + np.pi) % (2 * np.pi) - np.pi

        # find the circle (radius, x_center, y_center) formed from three hits of interest
        BAD_CHI2 = 1e6
        i0, i1, i2 = 0, 4, 7
        ixs = [1, 2, 3, 5, 6]
        circle_d = 2 * (new[f"t4_x_{i0}"] * (new[f"t4_y_{i1}"] - new[f"t4_y_{i2}"]) +
                        new[f"t4_x_{i1}"] * (new[f"t4_y_{i2}"] - new[f"t4_y_{i0}"]) +
                        new[f"t4_x_{i2}"] * (new[f"t4_y_{i0}"] - new[f"t4_y_{i1}"]))
        circle_x = np.divide(new[f"t4_r_{i0}"]**2 * (new[f"t4_y_{i1}"] - new[f"t4_y_{i2}"]) +
                             new[f"t4_r_{i1}"]**2 * (new[f"t4_y_{i2}"] - new[f"t4_y_{i0}"]) +
                             new[f"t4_r_{i2}"]**2 * (new[f"t4_y_{i0}"] - new[f"t4_y_{i1}"]),
                             circle_d)
        circle_y = np.divide(new[f"t4_r_{i0}"]**2 * (new[f"t4_x_{i2}"] - new[f"t4_x_{i1}"]) +
                             new[f"t4_r_{i1}"]**2 * (new[f"t4_x_{i0}"] - new[f"t4_x_{i2}"]) +
                             new[f"t4_r_{i2}"]**2 * (new[f"t4_x_{i1}"] - new[f"t4_x_{i0}"]),
                             circle_d)
        circle_r = np.sqrt((new[f"t4_x_{i0}"] - circle_x)**2 + (new[f"t4_y_{i0}"] - circle_y)**2)
        circle_ok = circle_d != 0
        if np.any(~circle_ok):
            logger.warning(f"Found {np.sum(~circle_ok)} invalid circles with circle_d = 0")

        # calculate the average diff
        diff2s = []
        for ix in ixs:
            circle_diff = np.sqrt((new[f"t4_x_{ix}"] - circle_x)**2 + (new[f"t4_y_{ix}"] - circle_y)**2) - circle_r
            diff2s.append(np.where(circle_ok, circle_diff**2, BAD_CHI2))
        new[f"t4_chi2_xy_{i0}{i1}{i2}"] = np.sum(diff2s, axis=0)

        # calculate chi2 for sz fit, where s is the arc length along the circle
        phis = [ np.arctan2(new[f"t4_y_{it}"] - circle_y, new[f"t4_x_{it}"] - circle_x) for it in range(N_LAYERS_IN_T4) ]
        dphis = [ (phi - phis[0] + np.pi) % (2 * np.pi) - np.pi for phi in phis ]
        pathlengths = [ circle_r * dphi for dphi in dphis ]
        for it in range(N_LAYERS_IN_T4):
            new[f"t4_s_{it}"] = pathlengths[it]

        # -------------------------- <Claude derivation> --------------------------
        # stack per-hit columns
        s_all = np.stack(pathlengths, axis=1)
        z_all = np.stack([ new[f"t4_z_{it}"] for it in range(N_LAYERS_IN_T4) ], axis=1)

        # row-wise least-squares line: z = z0_ref + tanlambda * s
        s_mean = s_all.mean(axis=1, keepdims=True)
        z_mean = z_all.mean(axis=1, keepdims=True)
        ds, dz = s_all - s_mean, z_all - z_mean
        s_ss = (ds * ds).sum(axis=1)                                 # (N,)
        s_sz = (ds * dz).sum(axis=1)                                 # (N,)
        tanlambda = s_sz / s_ss                                      # (N,)
        z0_ref    = z_mean.ravel() - tanlambda * s_mean.ravel()      # (N,)

        # residuals + quality metric, one per track
        resid = z_all - (z0_ref[:, None] + tanlambda[:, None] * s_all)
        resid2 = (resid ** 2).sum(axis=1)
        new[f"t4_chi2_sz"] = np.where(circle_ok, resid2, BAD_CHI2)
        # -------------------------- </Claude derivation> --------------------------

        # rename some things
        rename = {
            "ls_system": "t4_system",
            "ls_system_lower": "t4_system_lower",
            "ls_system_upper": "t4_system_upper",
            "ls_doublelayer_lower": "t4_doublelayer_lower",
            "ls_doublelayer_upper": "t4_doublelayer_upper",
            "ls_gdoublelayer_lower": "t4_gdoublelayer_lower",
            "ls_gdoublelayer_upper": "t4_gdoublelayer_upper",
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
        gdl_l = t4s["t4_gdoublelayer_lower"]
        gdl_u = t4s["t4_gdoublelayer_upper"]
        t4s["t4_ok_dz"] = np.abs(t4s["t4_dz"]) < self.T4_DZ_CUT[gdl_l, gdl_u]
        t4s["t4_ok_dr"] = np.abs(t4s["t4_dr"]) < self.T4_DR_CUT[gdl_l, gdl_u]
        t4s["t4_ok_dphi"] = np.abs(new["t4_dphi"]) < np.pi / 2.0
        t4s["t4_ok_dthetarz"] = np.abs(new["t4_dtheta_rz"]) < self.T4_DTHETA_RZ_CUT[gdl_l, gdl_u]
        t4s["t4_ok_chi2_xy"] = np.abs(new[f"t4_chi2_xy_{i0}{i1}{i2}"]) < self.T4_CHI2_XY_CUT[gdl_l, gdl_u]
        t4s["t4_ok_chi2_sz"] = np.abs(new[f"t4_chi2_sz"]) < self.T4_CHI2_SZ_CUT[gdl_l, gdl_u]
        t4s["t4_ok"] = (
            t4s["t4_ok_dphi"] &
            t4s["t4_ok_dz"] &
            t4s["t4_ok_dr"] &
            t4s["t4_ok_dthetarz"] &
            t4s["t4_ok_chi2_xy"] &
            # t4s["t4_ok_chi2_sz"] &
            np.ones(len(t4s), dtype=bool)
        )
        if self.signal:
            t4s["t4_ok_mcp"] = t4s["i_mcp"] != NO_MCP
            t4s["t4_ok_first_exit"] = t4s["t4_first_exit"] == True
            t4s["t4_ok_from_fiducial_mcp"] = t4s["t4_from_fiducial_mcp"] == True
            t4s["t4_ok_yanxi"] = (
                t4s["t4_ok_dphi"] &
                t4s["t4_ok_dz"] &
                t4s["t4_ok_dr"] &
                t4s["t4_ok_dthetarz"] &
                t4s["t4_ok_chi2_xy"] &
                t4s["t4_ok_first_exit"] &
                t4s["t4_ok_from_fiducial_mcp"] &
                t4s["t4_ok_mcp"] &
                np.ones(len(t4s), dtype=bool)
            )

        # merge new columns into the df
        t4s = pd.concat([t4s, pd.DataFrame(new, index=t4s.index)], axis=1)

        # remove as desired
        for cut in [col for col in t4s.columns if col.startswith("t4_ok")]:
            cutflow[cut] = np.sum(t4s[cut])
        if self.cut_t4s:
            t4s = t4s[t4s["t4_ok"]]

        # fin
        return t4s, cutflow


