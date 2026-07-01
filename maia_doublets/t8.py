"""
This module defines the T8Maker class, which creates a T8 from two T4s.
All T4s in layers 0-3 are considered for combination with T4s in layers 4-7.
The T4s are combined if they satisfy goodness criteria.
To avoid filling the memory with all possible combinations, combination is done with a groupby approach.

"""

import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)

from maia_doublets.constants import BYTE_TO_MB, NO_MCP
from maia_doublets.constants import N_T8_PHI_SLICES
from maia_doublets.constants import N_LAYERS_IN_T8

class T8Maker:

    def __init__(
        self,
        signal: bool,
        cut_t8s: bool,
        calibs: dict,
        t4s: pd.DataFrame,
    ):
        self.df = None
        self.signal = signal
        self.cut_t8s = cut_t8s
        self.t4s = t4s.copy()
        memory = self.t4s.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Making T8s. T4 dataframe size: {memory:.2f} MB")

        self.T8_DZ_CUT = calibs.get("t8_dz", np.zeros((10, 10)))
        self.T8_DR_CUT = calibs.get("t8_dr", np.zeros((10, 10)))

        # how to merge lower and upper T4s into T8s
        self.merge_keys = [
            "file",
            "i_event",
            "t4_phi_slice",
            "t4_eta_slice",
        ]

        self.filter_t4s()
        self.sort_t4s()
        self.make_t8s()


    def filter_t4s(self):
        # only consider "good" t4s
        logger.info("Filtering T4s for T8s ...")
        self.t4s = self.t4s[ self.t4s["t4_ok"] ]
        memory = self.t4s.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Memory usage after filtering T4s: {memory:.1f} MB")


    def sort_t4s(self):
        # sort T4s intuitively
        logger.info("Sorting T4s ...")
        cols = [
            "file",
            "i_event",
            "t4_system_lower",
            "t4_system_upper",
            "t4_doublelayer",
        ]
        print("len(self.t4s) = ", len(self.t4s))
        print("self.t4s.cols = ", self.t4s.columns)
        self.t4s = self.t4s.sort_values(by=cols).reset_index(drop=True)


    def make_t8s(self) -> None:
        """
        Todo: add description of this T8-making algorithm
        """

        # split t4s into global doublelayer once, up front
        logger.info(f"Splitting T4s ...")
        t4s = {
            gdl: group.reset_index(drop=True)
            for gdl, group
            in self.t4s.groupby("t4_gdoublelayer", sort=True)
        }

        # make T8s from neighboring doublelayers
        gdoublelayer_pairs = [
            (0, 4), # IT, OT
        ]
        all_t8s, all_cutflows = [], []
        for (lower, upper) in gdoublelayer_pairs:
            logger.info(f"Making T8s from gdl={lower} and gdl={upper} ...")
            lower_df = t4s[lower]
            upper_df = t4s[upper]
            t8s, cutflow = self.make_t8s_from_lower_upper(lower_df, upper_df)
            for col in cutflow:
                logger.info(f"T8s cutflow, gdl={lower}-{upper}, {col}: {cutflow[col]}")
            all_t8s.append(t8s)
            all_cutflows.append(cutflow)

        # merge dataframes
        logger.info(f"Merging {len(all_t8s)} groups of T8s ...")
        self.df = pd.concat(all_t8s, ignore_index=True)

        # merge cutflow
        cutflow = pd.DataFrame(all_cutflows)
        for col in cutflow.columns:
            logger.info(f"T8s cutflow, {col}: {cutflow[col].sum()}")

        # announce memory
        memory = self.df.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Memory usage of T8s: {memory:.1f} MB")


    def make_t8s_from_lower_upper(self, lower: pd.DataFrame, upper: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
        """
        Given an lower df and upper df, make T8s from their combinations
            For example, lower df could be T4s which span OT layers 0-3 (global doublelayer 4),
            and upper df could be T4s which span OT layers 4-7 (global doublelayer 6).
        """
        # get combinations of lower and upper
        # within neighboring phi/eta slices
        logger.info(f"Total number of lower T4s: {len(lower)}")
        logger.info(f"Total number of upper T4s: {len(upper)}")
        cands = []
        for phi_shift in (-1, 0, 1):
            for eta_shift in (-1, 0, 1):
                shifted = upper.assign(
                    t4_phi_slice=(upper["t4_phi_slice"] + phi_shift) % N_T8_PHI_SLICES,
                    t4_eta_slice=(upper["t4_eta_slice"] + eta_shift),
                )
                cands.append(lower.merge(shifted, on=self.merge_keys, suffixes=("_lower", "_upper")))

        # combine candidates into one dataframe
        t4s = pd.concat(cands, ignore_index=True)

        # calculate T8 features and cuts
        t8s, cutflow = self.consolidate_t8_features(t4s)
        return t8s, cutflow


    def consolidate_t8_features(self, t8s: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:

        # the doublelayer
        t8s["t8_doublelayer"] = t8s["t4_doublelayer_lower"]
        t8s["t8_gdoublelayer"] = t8s["t4_gdoublelayer_lower"]

        # rz projection
        slope_rz = np.divide(t8s["t4_z_upper"] - t8s["t4_z_lower"],
                             t8s["t4_r_upper"] - t8s["t4_r_lower"])
        t8s["t8_dz"] = t8s["t4_z_lower"] - t8s["t4_r_lower"] * slope_rz

        # xy projection
        slope_xy = np.divide(t8s["t4_y_upper"] - t8s["t4_y_lower"],
                             t8s["t4_x_upper"] - t8s["t4_x_lower"])
        intercept_xy = t8s["t4_y_lower"] - t8s["t4_x_lower"] * slope_xy
        t8s["t8_dr"] = np.abs(intercept_xy) / np.sqrt(1 + slope_xy**2)

        # assign truth info
        mcp_ok = t8s["i_mcp_lower"] == t8s["i_mcp_upper"]
        t8s["i_mcp"] = t8s["i_mcp_lower"].where(mcp_ok, NO_MCP)
        if self.signal:
            t8s["t8_first_exit"] = t8s["t4_first_exit_lower"] & t8s["t4_first_exit_upper"]
            t8s["t8_from_fiducial_mcp"] = t8s["t4_from_fiducial_mcp_lower"] & t8s["t4_from_fiducial_mcp_upper"]
            t8s["t8_detectable"] = mcp_ok & t8s["t4_detectable_lower"] & t8s["t4_detectable_upper"]
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
                t8s[attr] = t8s[f"{attr}_lower"].where(mcp_ok, 0)

        # pass-through the simhit positions
        for coord in ["x", "y", "r", "z"]:
            t8s[f"t8_{coord}_0"] = t8s[f"t4_{coord}_0_lower"]
            t8s[f"t8_{coord}_1"] = t8s[f"t4_{coord}_1_lower"]
            t8s[f"t8_{coord}_2"] = t8s[f"t4_{coord}_2_lower"]
            t8s[f"t8_{coord}_3"] = t8s[f"t4_{coord}_3_lower"]
            t8s[f"t8_{coord}_4"] = t8s[f"t4_{coord}_4_lower"]
            t8s[f"t8_{coord}_5"] = t8s[f"t4_{coord}_5_lower"]
            t8s[f"t8_{coord}_6"] = t8s[f"t4_{coord}_6_lower"]
            t8s[f"t8_{coord}_7"] = t8s[f"t4_{coord}_7_lower"]
            t8s[f"t8_{coord}_8"] = t8s[f"t4_{coord}_0_upper"]
            t8s[f"t8_{coord}_9"] = t8s[f"t4_{coord}_1_upper"]
            t8s[f"t8_{coord}_10"] = t8s[f"t4_{coord}_2_upper"]
            t8s[f"t8_{coord}_11"] = t8s[f"t4_{coord}_3_upper"]
            t8s[f"t8_{coord}_12"] = t8s[f"t4_{coord}_4_upper"]
            t8s[f"t8_{coord}_13"] = t8s[f"t4_{coord}_5_upper"]
            t8s[f"t8_{coord}_14"] = t8s[f"t4_{coord}_6_upper"]
            t8s[f"t8_{coord}_15"] = t8s[f"t4_{coord}_7_upper"]

        # more features
        t8s["t8_deta"] = t8s["t4_eta_upper"] - t8s["t4_eta_lower"]
        t8s["t8_dphi"] = t8s["t4_phi_upper"] - t8s["t4_phi_lower"]
        t8s["t8_dphi"] = (t8s["t8_dphi"] + np.pi) % (2 * np.pi) - np.pi
        t8s["t8_eta"] = 0.5 * (t8s["t4_eta_upper"] + t8s["t4_eta_lower"])
        t8s["t8_x"] = 0.5 * (t8s["t4_x_upper"] + t8s["t4_x_lower"])
        t8s["t8_y"] = 0.5 * (t8s["t4_y_upper"] + t8s["t4_y_lower"])
        t8s["t8_phi"] = np.arctan2(t8s["t8_y"], t8s["t8_x"])

        # angle differences (handle wraparound)
        # t8s["t8_dtheta_rz"] = t8s["t4_theta_rz_upper"] - t8s["t4_theta_rz_lower"]
        # t8s["t8_dtheta_rz"] = (t8s["t8_dtheta_rz"] + np.pi) % (2 * np.pi) - np.pi

        # collect new columns
        newcols = {}

        # find the circle (radius, x_center, y_center) formed from three hits of interest
        BAD_CHI2 = 1e6
        i0, i1, i2 = 0, 4, 7
        ixs = [1, 2, 3, 5, 6]
        circle_d = 2 * (t8s[f"t8_x_{i0}"] * (t8s[f"t8_y_{i1}"] - t8s[f"t8_y_{i2}"]) +
                        t8s[f"t8_x_{i1}"] * (t8s[f"t8_y_{i2}"] - t8s[f"t8_y_{i0}"]) +
                        t8s[f"t8_x_{i2}"] * (t8s[f"t8_y_{i0}"] - t8s[f"t8_y_{i1}"]))
        circle_x = np.divide(t8s[f"t8_r_{i0}"]**2 * (t8s[f"t8_y_{i1}"] - t8s[f"t8_y_{i2}"]) +
                                t8s[f"t8_r_{i1}"]**2 * (t8s[f"t8_y_{i2}"] - t8s[f"t8_y_{i0}"]) +
                                t8s[f"t8_r_{i2}"]**2 * (t8s[f"t8_y_{i0}"] - t8s[f"t8_y_{i1}"]),
                                circle_d)
        circle_y = np.divide(t8s[f"t8_r_{i0}"]**2 * (t8s[f"t8_x_{i2}"] - t8s[f"t8_x_{i1}"]) +
                                t8s[f"t8_r_{i1}"]**2 * (t8s[f"t8_x_{i0}"] - t8s[f"t8_x_{i2}"]) +
                                t8s[f"t8_r_{i2}"]**2 * (t8s[f"t8_x_{i1}"] - t8s[f"t8_x_{i0}"]),
                                circle_d)
        circle_r = np.sqrt((t8s[f"t8_x_{i0}"] - circle_x)**2 + (t8s[f"t8_y_{i0}"] - circle_y)**2)
        circle_ok = circle_d != 0
        if np.any(~circle_ok):
            logger.warning(f"Found {np.sum(~circle_ok)} invalid circles with circle_d = 0")

        # calculate the average diff
        diff2s = []
        for ix in ixs:
            circle_diff = np.sqrt((t8s[f"t8_x_{ix}"] - circle_x)**2 + (t8s[f"t8_y_{ix}"] - circle_y)**2) - circle_r
            diff2s.append(np.where(circle_ok, circle_diff**2, BAD_CHI2))
        newcols[f"t8_chi2_xy_{i0}{i1}{i2}"] = np.sum(diff2s, axis=0)

        # calculate chi2 for sz fit, where s is the arc length along the circle
        phis = [ np.arctan2(t8s[f"t8_y_{it}"] - circle_y, t8s[f"t8_x_{it}"] - circle_x) for it in range(N_LAYERS_IN_T8) ]
        dphis = [ (phi - phis[0] + np.pi) % (2 * np.pi) - np.pi for phi in phis ]
        pathlengths = [ circle_r * dphi for dphi in dphis ]
        for it in range(N_LAYERS_IN_T8):
            newcols[f"t8_s_{it}"] = pathlengths[it]

        # -------------------------- <Claude derivation> --------------------------
        # stack per-hit columns
        s_all = np.stack(pathlengths, axis=1)
        z_all = np.stack([ t8s[f"t8_z_{it}"] for it in range(N_LAYERS_IN_T8) ], axis=1)

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
        newcols[f"t8_chi2_sz"] = np.where(circle_ok, resid2, BAD_CHI2)
        # -------------------------- </Claude derivation> --------------------------

        # rename some things
        rename = {
            "t4_system": "t8_system",
            "t4_system_lower": "t8_system_lower",
            "t4_system_upper": "t8_system_upper",
            "t4_doublelayer_lower": "t8_doublelayer_lower",
            "t4_doublelayer_upper": "t8_doublelayer_upper",
            "t4_gdoublelayer_lower": "t8_gdoublelayer_lower",
            "t4_gdoublelayer_upper": "t8_gdoublelayer_upper",
            "t4_module_lower": "t8_module_lower",
            "t4_module_upper": "t8_module_upper",
            "t4_sensor_lower": "t8_sensor_lower",
            "t4_sensor_upper": "t8_sensor_upper",
            "t4_dr_lower": "t8_dr_lower",
            "t4_dr_upper": "t8_dr_upper",
            "t4_dz_lower": "t8_dz_lower",
            "t4_dz_upper": "t8_dz_upper",
            "t4_ok_lower": "t8_t4_ok_lower",
            "t4_ok_upper": "t8_t4_ok_upper",
        }
        t8s = t8s.rename(columns=rename)

        # drop other cols
        dropcols = ["i_mcp_lower", "i_mcp_upper"]
        dropcols.extend([col for col in t8s.columns if col.startswith("simhit_")])
        dropcols.extend([col for col in t8s.columns if col.startswith("doublet_")])
        dropcols.extend([col for col in t8s.columns if col.startswith("ls_")])
        dropcols.extend([col for col in t8s.columns if col.startswith("t2_")])
        dropcols.extend([col for col in t8s.columns if col.startswith("t4_")])
        dropcols.extend([col for col in t8s.columns if col.startswith("mcp_") and col.endswith("_lower")])
        dropcols.extend([col for col in t8s.columns if col.startswith("mcp_") and col.endswith("_upper")])
        t8s.drop(columns=dropcols, errors="ignore", inplace=True)

        # record some numbers
        cutflow = {"all": len(t8s)}

        # record some cut results
        gdl_l = t8s["t8_gdoublelayer_lower"]
        gdl_u = t8s["t8_gdoublelayer_upper"]
        t8s["t8_ok_dz"] = np.abs(t8s["t8_dz"]) < self.T8_DZ_CUT[gdl_l, gdl_u]
        t8s["t8_ok_dr"] = np.abs(t8s["t8_dr"]) < self.T8_DR_CUT[gdl_l, gdl_u]
        t8s["t8_ok"] = (
            t8s["t8_ok_dz"] &
            t8s["t8_ok_dr"] &
            np.ones(len(t8s), dtype=bool)
        )
        if self.signal:
            t8s["t8_ok_mcp"] = t8s["i_mcp"] != NO_MCP
            t8s["t8_ok_first_exit"] = t8s["t8_first_exit"] == True
            t8s["t8_ok_from_fiducial_mcp"] = t8s["t8_from_fiducial_mcp"] == True
            t8s["t8_ok_yanxi"] = (
                t8s["t8_ok_dz"] &
                t8s["t8_ok_dr"] &
                t8s["t8_ok_first_exit"] &
                t8s["t8_ok_from_fiducial_mcp"] &
                t8s["t8_ok_mcp"] &
                np.ones(len(t8s), dtype=bool)
            )

        # merge new columns into the df
        t8s = pd.concat([t8s, pd.DataFrame(newcols, index=t8s.index)], axis=1)

        # remove as desired
        for cut in [col for col in t8s.columns if col.startswith("t8_ok")]:
            cutflow[cut] = np.sum(t8s[cut])
        if self.cut_t8s:
            t8s = t8s[t8s["t8_ok"]]

        # fin
        return t8s, cutflow


