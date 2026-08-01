import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)

from maia_doublets.constants import BYTE_TO_MB, NO_MCP
from maia_doublets.constants import N_T2_PHI_SLICES
from maia_doublets.constants import DETECTOR_MAX_PHI, DETECTOR_MAX_ETA
from maia_doublets.constants import N_T4_PHI_SLICES, N_T4_ETA_SLICES
from maia_doublets.constants import BAD_CHI2, N_LAYERS_IN_T2

class T2Maker:

    #
    # To make doublets, we do 2 groupbys:
    #  Layers 01, 23, ... grouped by doublet_doublelayer_mod_2
    #  Layers 12, 34, ... grouped by doublet_doublelayer_plus_1_mod_2
    #

    def __init__(
            self,
            signal: bool,
            cut_t2s: bool,
            calibs: dict,
            mds: pd.DataFrame,
        ):
        self.df = None
        self.cutflow = None
        self.signal = signal
        self.cut_t2s = cut_t2s
        self.lower_suffix = "lower"
        self.upper_suffix = "upper"
        memory = mds.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Making t2s with mds memory {memory:.1f} MB ...")

        self.T2_DZ_CUT = calibs.get("t2_dz", np.zeros((10, 10)))
        self.T2_DR_CUT = calibs.get("t2_dr", np.zeros((10, 10)))
        self.T2_DTHETA_RZ_CUT = calibs.get("t2_dtheta_rz", np.zeros((10, 10)))
        self.T2_CHI2_XY_CUT = calibs.get("t2_chi2_xy", np.zeros((10, 10)))
        self.T2_CHI2_SZ_CUT = calibs.get("t2_chi2_sz", np.zeros((10, 10)))

        self.merge_keys = [
            "file",
            "i_event",
            "doublet_phi_slice",
            "doublet_eta_slice",
        ]

        self.mds = mds.copy()
        self.filter_mds()
        self.sort_mds()
        self.make_t2s()


    def filter_mds(self):
        # only consider "good" mds
        logger.info("Filtering MDs for T2s ...")
        self.mds = self.mds[ self.mds["doublet_ok"] ]
        memory = self.mds.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"Memory usage after filtering MDs: {memory:.1f} MB")


    def sort_mds(self):
        # sort mds by file, event, system, doublelayer, sensor, module
        logger.info("Sorting MDs for T2s ...")
        cols = [
            "file",
            "i_event",
            "doublet_system",
            "doublet_doublelayer",
            "doublet_sensor",
            "doublet_module",
        ]
        self.mds = self.mds.sort_values(by=cols).reset_index(drop=True)


    def make_t2s(self):

        # split mds into global doublelayer once, up front
        logger.info(f"Splitting MDs ...")
        mds = {
            gdl: group.reset_index(drop=True)
            for gdl, group
            in self.mds.groupby("doublet_gdoublelayer", sort=True)
        }

        # make T2s from neighboring doublelayers
        gdoublelayer_pairs = [
            (0, 1), # IT
            (2, 3), # IT
            (4, 5), # OT
            (6, 7), # OT
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
        self.cutflow = pd.DataFrame(all_cutflows)
        for col in self.cutflow.columns:
            logger.info(f"T2s cutflow, {col}: {self.cutflow[col].sum()}")

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
                if self.cut_t2s:
                    t2s = self.add_basic_t2_features(t2s)
                    sy = t2s["t2_system"]
                    dl = t2s["t2_doublelayer"]
                    t2s["t2_ok_dz"] = np.abs(t2s["t2_dz"]) < self.T2_DZ_CUT[sy, dl]
                    t2s["t2_ok_dr"] = np.abs(t2s["t2_dr"]) < self.T2_DR_CUT[sy, dl]
                    t2s = t2s[t2s["t2_ok_dz"] & t2s["t2_ok_dr"]]

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
            "doublet_doublelayer_lower": "t2_doublelayer_lower",
            "doublet_doublelayer_upper": "t2_doublelayer_upper",
            "doublet_glayer_lower": "t2_glayer_lower",
            "doublet_glayer_upper": "t2_glayer_upper",
            "doublet_module_lower": "t2_module_lower",
            "doublet_module_upper": "t2_module_upper",
            "doublet_sensor_lower": "t2_sensor_lower",
            "doublet_sensor_upper": "t2_sensor_upper",
            "doublet_system_lower": "t2_system_lower",
            "doublet_system_upper": "t2_system_upper",
            "doublet_dr_lower": "t2_dr_lower",
            "doublet_dr_upper": "t2_dr_upper",
            "doublet_dz_lower": "t2_dz_lower",
            "doublet_dz_upper": "t2_dz_upper",
            "doublet_ok_lower": "t2_md_ok_lower",
            "doublet_ok_upper": "t2_md_ok_upper",
        }
        return t2s.rename(columns=rename, copy=False)


    def add_basic_t2_features(self, t2s: pd.DataFrame) -> pd.DataFrame:

        # the system and doublelayer
        t2s["t2_system"] = t2s["t2_system_lower"].astype(np.uint8)
        t2s["t2_doublelayer"] = t2s["t2_doublelayer_lower"].astype(np.uint8)

        # rz projection
        slope_rz = np.divide(t2s["doublet_z_upper"] - t2s["doublet_z_lower"],
                             t2s["doublet_r_upper"] - t2s["doublet_r_lower"])
        t2s["t2_dz"] = t2s["doublet_z_lower"] - t2s["doublet_r_lower"] * slope_rz
        t2s["t2_theta_rz"] = np.arctan(slope_rz)

        # xy projection
        slope_xy = np.divide(t2s["doublet_y_upper"] - t2s["doublet_y_lower"],
                             t2s["doublet_x_upper"] - t2s["doublet_x_lower"])
        intercept_xy = t2s["doublet_y_lower"] - t2s["doublet_x_lower"] * slope_xy
        t2s["t2_dr"] = np.abs(intercept_xy) / np.sqrt(1 + slope_xy**2)

        return t2s


    def consolidate_t2_features(self, t2s: pd.DataFrame) -> tuple[pd.DataFrame, dict]:

        t2s = self.add_basic_t2_features(t2s)

        # assign truth info
        mcp_ok = t2s["i_mcp_lower"] == t2s["i_mcp_upper"]
        t2s["i_mcp"] = t2s["i_mcp_lower"].where(mcp_ok, NO_MCP)
        if self.signal:
            t2s["t2_first_exit"] = t2s["doublet_first_exit_lower"] & t2s["doublet_first_exit_upper"]
            t2s["t2_from_fiducial_mcp"] = t2s["doublet_from_fiducial_mcp_lower"] & t2s["doublet_from_fiducial_mcp_upper"]
            t2s["t2_detectable"] = mcp_ok & t2s["doublet_detectable_lower"] & t2s["doublet_detectable_upper"]
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
                t2s[attr] = t2s[f"{attr}_lower"].where(mcp_ok, 0)

        # collect new columns
        new = {}

        # features: position
        new["t2_r"] = (t2s["doublet_r_lower"] + t2s["doublet_r_upper"]) / 2
        new["t2_z"] = (t2s["doublet_z_lower"] + t2s["doublet_z_upper"]) / 2
        new["t2_x"] = (t2s["doublet_x_lower"] + t2s["doublet_x_upper"]) / 2
        new["t2_y"] = (t2s["doublet_y_lower"] + t2s["doublet_y_upper"]) / 2
        new["t2_phi"] = np.arctan2(new["t2_y"], new["t2_x"])
        new["t2_theta"] = np.arctan2(new["t2_r"], new["t2_z"])
        new["t2_eta"] = -np.log(np.tan(new["t2_theta"] / 2))
        new["t2_phi_slice"] = np.floor((new["t2_phi"] + DETECTOR_MAX_PHI) / (2 * DETECTOR_MAX_PHI) * N_T4_PHI_SLICES).astype(np.int16)
        new["t2_eta_slice"] = np.floor((new["t2_eta"] + DETECTOR_MAX_ETA) / (2 * DETECTOR_MAX_ETA) * N_T4_ETA_SLICES).astype(np.int16)

        # assign more features
        new["t2_ddr"] = t2s["t2_dr_upper"] - t2s["t2_dr_lower"]
        new["t2_ddz"] = t2s["t2_dz_upper"] - t2s["t2_dz_lower"]
        new["t2_deta"] = t2s["doublet_eta_upper"] - t2s["doublet_eta_lower"]
        new["t2_dphi"] = t2s["doublet_phi_upper"] - t2s["doublet_phi_lower"]
        new["t2_dphi"] = (new["t2_dphi"] + np.pi) % (2 * np.pi) - np.pi
        new["t2_dqoverpt"] = t2s["doublet_qoverpt_upper"] - t2s["doublet_qoverpt_lower"]

        # pass-through the simhit positions as float64 for now
        for coord in ["x", "y", "r", "z"]:
            new[f"t2_{coord}_0"] = t2s[f"doublet_{coord}_0_lower"].astype(np.float64)
            new[f"t2_{coord}_1"] = t2s[f"doublet_{coord}_1_lower"].astype(np.float64)
            new[f"t2_{coord}_2"] = t2s[f"doublet_{coord}_0_upper"].astype(np.float64)
            new[f"t2_{coord}_3"] = t2s[f"doublet_{coord}_1_upper"].astype(np.float64)

        # angle differences (handle wraparound)
        new["t2_dtheta_rz"] = t2s["doublet_theta_rz_upper"] - t2s["doublet_theta_rz_lower"]
        new["t2_dtheta_xy"] = t2s["doublet_theta_xy_upper"] - t2s["doublet_theta_xy_lower"]
        new["t2_dtheta_rz"] = (new["t2_dtheta_rz"] + np.pi) % (2 * np.pi) - np.pi
        new["t2_dtheta_xy"] = (new["t2_dtheta_xy"] + np.pi) % (2 * np.pi) - np.pi

        # find the circle (radius, x_center, y_center) formed from the first three hits
        circle_d = 2 * (new["t2_x_0"] * (new["t2_y_1"] - new["t2_y_2"]) +
                        new["t2_x_1"] * (new["t2_y_2"] - new["t2_y_0"]) +
                        new["t2_x_2"] * (new["t2_y_0"] - new["t2_y_1"]))
        circle_x = np.divide(new["t2_r_0"]**2 * (new["t2_y_1"] - new["t2_y_2"]) +
                             new["t2_r_1"]**2 * (new["t2_y_2"] - new["t2_y_0"]) +
                             new["t2_r_2"]**2 * (new["t2_y_0"] - new["t2_y_1"]),
                             circle_d)
        circle_y = np.divide(new["t2_r_0"]**2 * (new["t2_x_2"] - new["t2_x_1"]) +
                             new["t2_r_1"]**2 * (new["t2_x_0"] - new["t2_x_2"]) +
                             new["t2_r_2"]**2 * (new["t2_x_1"] - new["t2_x_0"]),
                             circle_d)
        circle_r = np.sqrt((new["t2_x_0"] - circle_x)**2 + (new["t2_y_0"] - circle_y)**2)
        circle_ok = circle_d != 0
        if np.any(~circle_ok):
            logger.warning(f"Found {np.sum(~circle_ok)} invalid circles with circle_d = 0")
        circle_diff = np.sqrt((new["t2_x_3"] - circle_x)**2 + (new["t2_y_3"] - circle_y)**2) - circle_r

        # calculate the distance from (x_3, y_3) to the circle
        new["t2_chi2_xy"] = np.where(circle_ok, circle_diff**2, BAD_CHI2)

        # calculate chi2 for sz fit, where s is the arc length along the circle
        phis = [ np.arctan2(new[f"t2_y_{it}"] - circle_y, new[f"t2_x_{it}"] - circle_x) for it in range(N_LAYERS_IN_T2) ]
        dphis = [ (phi - phis[0] + np.pi) % (2 * np.pi) - np.pi for phi in phis ]
        pathlengths = [ circle_r * dphi for dphi in dphis ]
        for it in range(N_LAYERS_IN_T2):
            new[f"t2_s_{it}"] = pathlengths[it]

        # -------------------------- <Claude derivation> --------------------------
        # stack per-hit columns
        s_all = np.stack(pathlengths, axis=1)
        z_all = np.stack([ new[f"t2_z_{it}"] for it in range(N_LAYERS_IN_T2) ], axis=1)

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
        new[f"t2_chi2_sz"] = np.where(circle_ok, resid2, BAD_CHI2)
        # -------------------------- </Claude derivation> --------------------------

        # downscope the simhit positions to float32
        for coord in ["x", "y", "r", "z"]:
            new[f"t2_{coord}_0"] = new[f"t2_{coord}_0"].astype(np.float32)
            new[f"t2_{coord}_1"] = new[f"t2_{coord}_1"].astype(np.float32)
            new[f"t2_{coord}_2"] = new[f"t2_{coord}_2"].astype(np.float32)
            new[f"t2_{coord}_3"] = new[f"t2_{coord}_3"].astype(np.float32)

        # assign features from lower doublet (arbitrary choice)
        t2s["t2_module"] = t2s["t2_module_lower"]
        t2s["t2_sensor"] = t2s["t2_sensor_lower"]
        t2s["t2_glayer"] = t2s["t2_glayer_lower"]
        t2s["t2_gdoublelayer"] = t2s["t2_glayer"] // 2

        # and drop other cols
        dropcols = ["i_mcp_lower", "i_mcp_upper"]
        dropcols.extend([col for col in t2s.columns if col.startswith("simhit_")])
        dropcols.extend([col for col in t2s.columns if col.startswith("doublet_")])
        dropcols.extend([col for col in t2s.columns if col.startswith("mcp_") and col.endswith("_lower")])
        dropcols.extend([col for col in t2s.columns if col.startswith("mcp_") and col.endswith("_upper")])
        t2s.drop(columns=dropcols, errors="ignore", inplace=True)

        # record some numbers
        cutflow = {"all": len(t2s)}

        # record some cut results
        sy = t2s["t2_system"]
        dl = t2s["t2_doublelayer"]
        t2s["t2_ok_dz"] = np.abs(t2s["t2_dz"]) < self.T2_DZ_CUT[sy, dl] # NB: these two cuts are special
        t2s["t2_ok_dr"] = np.abs(t2s["t2_dr"]) < self.T2_DR_CUT[sy, dl] # because theyre already defined above
        new["t2_ok_dtheta_rz"] = np.abs(new["t2_dtheta_rz"]) < self.T2_DTHETA_RZ_CUT[sy, dl]
        new["t2_ok_dphi"] = np.abs(new["t2_dphi"]) < np.pi / 2.0
        new["t2_ok_chi2_xy"] = np.abs(new["t2_chi2_xy"]) < self.T2_CHI2_XY_CUT[sy, dl]
        new["t2_ok_chi2_sz"] = np.abs(new["t2_chi2_sz"]) < self.T2_CHI2_SZ_CUT[sy, dl]
        new["t2_ok_drdz"] = t2s["t2_ok_dz"] & t2s["t2_ok_dr"] & new["t2_ok_dphi"]
        new["t2_ok_drdzdthetarz"] = t2s["t2_ok_dz"] & t2s["t2_ok_dr"] & new["t2_ok_dphi"] & new["t2_ok_dtheta_rz"]
        new["t2_ok"] = (
            t2s["t2_ok_dz"] &
            t2s["t2_ok_dr"] &
            new["t2_ok_dphi"] &
            new["t2_ok_chi2_xy"] &
            new["t2_ok_chi2_sz"] &
            np.ones(len(t2s), dtype=bool)
        )

        # merge new columns into the df
        t2s = pd.concat([t2s, pd.DataFrame(new, index=t2s.index)], axis=1)

        # remove as desired
        for cut in [col for col in t2s.columns if col.startswith("t2_ok")]:
            cutflow[cut] = np.sum(t2s[cut])
        if self.cut_t2s:
            t2s = t2s[t2s["t2_ok"]]

        # fin
        return t2s, cutflow


