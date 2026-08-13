import numpy as np
import pandas as pd
import time
import logging
logger = logging.getLogger(__name__)

from maia_doublets.constants import MAGNETIC_FIELD, SPEED_OF_LIGHT
from maia_doublets.constants import BYTE_TO_MB, NO_MCP
from maia_doublets.constants import N_T2_PHI_SLICES, N_T2_ETA_SLICES, DETECTOR_MAX_ETA, DETECTOR_MAX_PHI
from maia_doublets.constants import BAD_CHI2

class MDMaker:

    def __init__(
        self,
        signal: bool,
        cut_mds: bool,
        fast_merge: bool,
        calibs: dict,
        hits: pd.DataFrame,
    ):
        self.signal = signal
        self.cut_mds = cut_mds
        self.MD_DZ_CUT = calibs.get("md_dz", np.zeros((10, 10)))
        self.MD_DR_CUT = calibs.get("md_dr", np.zeros((10, 10)))
        self.md_cols = [
            "file",
            "i_event", # the event
            "hit_system", # the system (IT, OT)
            "hit_layer_div_2", # the double layer
            "hit_module", # the phi-module
            "hit_sensor", # the z-sensor
        ]
        self.fast_merge = fast_merge
        self.df, self.cutflow = self.make_mds(hits)


    def make_mds(self, df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        logger.info("Making mds ...")

        groupby_cols = [
            "file",
        ]
        if not self.signal:
            groupby_cols += [
                "i_event", # the event
                "hit_system", # the system (IT, OT)
                "hit_layer_div_2", # the double layer
                "hit_module", # the phi-module
                "hit_sensor", # the z-sensor
            ]

        # group loop
        logger.info("Grouping hits ...")
        groups = df.groupby(groupby_cols)
        all_mds, all_cutflows = [], []

        for i_group, (cols, group) in enumerate(groups):

            mds, cutflow = self.make_mds_from_group(group)

            all_mds.append(mds)
            all_cutflows.append(cutflow)

            if (self.signal and i_group % 100 == 0) or (not self.signal and i_group % 4 == 0):
                length = len(mds)
                size = mds.memory_usage(deep=True).sum() * BYTE_TO_MB
                logger.info(f"Processed group {i_group}/{len(groups)}, md size = {size:.1f} MB, n(mds) = {length} ...")

        # concatenate mds and cutflows
        logger.info(f"Concatenating mds ...")
        mds = pd.concat(all_mds, ignore_index=True)
        cutflow = pd.DataFrame(all_cutflows)
        for col in cutflow.columns:
            logger.info(f"Mds cutflow, {col}: {cutflow[col].sum()}")
        if len(mds) == 0:
            raise ValueError("No mds found in the DataFrame")

        # announcements
        logger.info(f"Total mds: {len(mds)}")
        logger.info(f"Total mds size: {mds.memory_usage(deep=True).sum() * BYTE_TO_MB:.1f} MB")
        counts = mds.groupby(["md_system",
                                   "md_doublelayer"]).size()
        for (system, doublelayer), total in counts.items():
            logger.info(f"n(mds) for system {system}, doublelayer {doublelayer}: {total}")

        return mds, cutflow


    def make_mds_from_group(self, group: pd.DataFrame) -> tuple[pd.DataFrame, dict]:

        lower_mask = group["hit_layer_mod_2"] == 0
        upper_mask = group["hit_layer_mod_2"] == 1
        lower = group[lower_mask]
        upper = group[upper_mask]

        # inner join to find mds
        if self.cut_mds and self.fast_merge:
            mds, n_full = self.merge_binned(lower, upper)
        else:
            mds = pd.merge(
                lower,
                upper,
                on=self.md_cols,
                how="inner",
                suffixes=("_lower", "_upper"),
            )
            n_full = len(mds)

        # md feature: xy, dr at point of closest approach to origin
        slope_xy = np.divide(mds["hit_y_upper"] - mds["hit_y_lower"],
                             mds["hit_x_upper"] - mds["hit_x_lower"])
        intercept_xy = mds["hit_y_lower"] - slope_xy * mds["hit_x_lower"]
        mds["md_dr"] = np.abs(intercept_xy) / np.sqrt(1 + slope_xy**2)

        # md feature: rz
        slope_rz = np.divide(mds["hit_z_upper"] - mds["hit_z_lower"],
                             mds["hit_r_upper"] - mds["hit_r_lower"])
        mds["md_dz"] = mds["hit_z_lower"] - mds["hit_r_lower"] * slope_rz
        mds["md_theta_rz"] = np.arctan(slope_rz)

        # deal with xy slope is NaN
        feature = "md_dr"
        isnull = mds[feature].isnull()
        if isnull.any():
            if np.any(isnull != (mds["hit_x_upper"] == mds["hit_x_lower"])):
                raise ValueError(f"Found NaN in {feature}, but not all are due to infinite slope")
            mds[feature] = mds[feature].fillna(mds["hit_x_lower"])

        # announce any nans
        for feature in ["md_dr", "md_dz"]:
            if mds[feature].isnull().any():
                n_nan = mds[feature].isnull().sum()
                msg = f"Found {n_nan} unexpected NaN in {feature}"
                logger.error(msg)
                raise ValueError(msg)

        # record some numbers
        cutflow = {"all": len(mds)}
        mask = {}

        # record some cut results
        sy = mds["hit_system"]
        dl = mds["hit_layer_div_2"]
        mds["md_ok_dr"] = np.abs(mds["md_dr"]) < self.MD_DR_CUT[sy, dl]
        mds["md_ok_dz"] = np.abs(mds["md_dz"]) < self.MD_DZ_CUT[sy, dl]
        mds["md_ok"] = mds["md_ok_dr"] & mds["md_ok_dz"]

        # remove as desired
        if self.cut_mds:
            cutflow["md_ok_dr"] = np.sum(mds["md_ok_dr"])
            cutflow["md_ok_dz"] = np.sum(mds["md_ok_dz"])
            cutflow["md_ok"] = np.sum(mds["md_ok"])
            mds = mds[mds["md_ok"]]

        # rename some columns
        rename = {
            "hit_system": "md_system",
            "hit_layer_div_2": "md_doublelayer",
            "hit_sensor": "md_sensor",
            "hit_module": "md_module",
        }
        mds = mds.rename(columns=rename)
        mds["md_glayer"] = mds["hit_glayer_lower"]
        mds["md_gdoublelayer"] = mds["hit_glayer_lower"] // 2

        # md feature, xy dphi
        phi_local = np.arctan2(mds["hit_y_upper"] - mds["hit_y_lower"],
                                mds["hit_x_upper"] - mds["hit_x_lower"])
        phi_global = np.arctan2((mds["hit_y_lower"] + mds["hit_y_upper"]) / 2.0,
                                (mds["hit_x_lower"] + mds["hit_x_upper"]) / 2.0)
        mds["md_dphi"] = phi_local - phi_global
        mds["md_dphi"] = (mds["md_dphi"] + np.pi) % (2 * np.pi) - np.pi
        mds["md_theta_xy"] = phi_local

        # md features: position
        mds["md_r"] = (mds["hit_r_lower"] + mds["hit_r_upper"]) / 2
        mds["md_z"] = (mds["hit_z_lower"] + mds["hit_z_upper"]) / 2
        mds["md_x"] = (mds["hit_x_lower"] + mds["hit_x_upper"]) / 2
        mds["md_y"] = (mds["hit_y_lower"] + mds["hit_y_upper"]) / 2
        mds["md_phi"] = np.arctan2(mds["md_y"], mds["md_x"])
        mds["md_theta"] = np.arctan2(mds["md_r"], mds["md_z"])
        mds["md_eta"] = -np.log(np.tan(mds["md_theta"] / 2))

        # divide the eta/phi space into slices, to be used in T2 seeding
        n_phi_slices = N_T2_PHI_SLICES[mds["md_system"]]
        n_eta_slices = N_T2_ETA_SLICES[mds["md_system"]]
        mds["md_phi_slice"] = np.floor((mds["md_phi"] + DETECTOR_MAX_PHI) / (2 * DETECTOR_MAX_PHI) * n_phi_slices).astype(np.int16)
        mds["md_eta_slice"] = np.floor((mds["md_eta"] + DETECTOR_MAX_ETA) / (2 * DETECTOR_MAX_ETA) * n_eta_slices).astype(np.int16)

        # guess charge from dphi:
        # positively charged particles have negative dphi, and vice versa
        mds["md_q"] = (-1*np.sign(mds["md_dphi"])).astype(np.int8)

        # pass-through the hit positions
        for coord in ["x", "y", "r", "z"]:
            mds[f"md_{coord}_0"] = mds[f"hit_{coord}_lower"]
            mds[f"md_{coord}_1"] = mds[f"hit_{coord}_upper"]

        # md feature: radius of circle composed of the two hits and the origin. R = abc/4K
        # then get pt from R
        circle_a = mds["hit_r_lower"]
        circle_b = mds["hit_r_upper"]
        circle_c = np.sqrt((mds["hit_x_upper"] - mds["hit_x_lower"])**2 +
                            (mds["hit_y_upper"] - mds["hit_y_lower"])**2)
        circle_K = 0.5 * np.abs(mds["hit_x_lower"] * mds["hit_y_upper"] -
                                mds["hit_x_upper"] * mds["hit_y_lower"])
        mds["md_circle_radius"] = np.divide(circle_a * circle_b * circle_c, 4.0 * circle_K)
        mds["md_pt"] = SPEED_OF_LIGHT * MAGNETIC_FIELD * mds["md_circle_radius"] * 1e-6
        mds["md_qoverpt"] = mds["md_q"] / mds["md_pt"]

        # md feature: truth info
        mcp_ok = mds["i_mcp_lower"] == mds["i_mcp_upper"]
        mds["i_mcp"] = mds["i_mcp_lower"].where(mcp_ok, NO_MCP)
        if self.signal:
            mds["md_first_exit"] = mds["hit_first_exit_lower"] & mds["hit_first_exit_upper"]
            mds["md_from_fiducial_mcp"] = mds["hit_from_fiducial_mcp_lower"] & mds["hit_from_fiducial_mcp_upper"]
            mds["md_detectable"] = mcp_ok & mds["hit_detectable_lower"] & mds["hit_detectable_upper"]
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
                mds[attr] = mds[f"{attr}_lower"].where(mcp_ok, 0)

        # drop columns which arent used downstream
        dropcols = ["i_mcp_lower", "i_mcp_upper"]
        dropcols.extend([col for col in mds.columns if col.startswith("hit_")])
        dropcols.extend([col for col in mds.columns if col.startswith("mcp_") and col.endswith("_lower")])
        dropcols.extend([col for col in mds.columns if col.startswith("mcp_") and col.endswith("_upper")])
        mds.drop(columns=dropcols, inplace=True)

        return mds, cutflow


    def merge_binned(self, lower: pd.DataFrame, upper: pd.DataFrame) -> tuple[pd.DataFrame, int]:
        """
        Experimental binned merge to speed up the md making, courtesy of Claude

        Binned-merge equivalent of pd.merge(lower, upper, on=md_cols,
        how="inner", suffixes=("_lower","_upper")), restricted to upper hits whose
        z bin is within +/-1 of the lower hit's predicted z bin. Returns
        (merge_equivalent_frame, n_full_crossproduct).

        Assumes a single cell (system / doublelayer / module / sensor constant).

        Bin width: within this cell the dz cut |md_dz| < DZ is exactly
        |z_lo*r_up - z_up*r_lo| < DZ*(r_up-r_lo), which confines z_up to an interval
        around the radial projection of the lower hit. The bin width is the widest
        such half-interval the cut allows in the cell, so the predicted bin plus its
        two neighbours are guaranteed to contain every dz survivor.
        """
        n_lo = len(lower)
        n_up = len(upper)
        n_full = n_lo * n_up

        system = int(lower["hit_system"].iloc[0])
        doublelayer = int(lower["hit_layer_div_2"].iloc[0])
        dz_cut = float(self.MD_DZ_CUT[system, doublelayer])

        z_lo = lower["hit_z"].to_numpy(np.float64)
        r_lo = lower["hit_r"].to_numpy(np.float64)
        z_up = upper["hit_z"].to_numpy(np.float64)
        r_up = upper["hit_r"].to_numpy(np.float64)

        r_up_min = r_up.min()
        r_up_max = r_up.max()

        # guard: closed form assumes the upper layer sits strictly outside the
        # lower one. If radii overlap, fall back to the exact full merge.
        if r_up_min <= r_lo.max():
            logger.warning(f"MD binned merge: weird data, falling back to full merge")
            mds = pd.merge(lower, upper, on=self.md_cols, how="inner",
                                suffixes=("_lower", "_upper"))
            return mds, n_full

        # allowed z_up interval per lower hit (union over r_up in [min, max])
        fmin_a = ((z_lo - dz_cut) * r_up_min + dz_cut * r_lo) / r_lo
        fmin_b = ((z_lo - dz_cut) * r_up_max + dz_cut * r_lo) / r_lo
        fmax_a = ((z_lo + dz_cut) * r_up_min - dz_cut * r_lo) / r_lo
        fmax_b = ((z_lo + dz_cut) * r_up_max - dz_cut * r_lo) / r_lo
        win_lo = np.minimum(fmin_a, fmin_b)
        win_hi = np.maximum(fmax_a, fmax_b)

        center = 0.5 * (win_lo + win_hi)            # = radial projection of the lower hit
        half = 0.5 * (win_hi - win_lo)              # half-window the dz cut permits
        bin_width = float(half.max())               # widest half-window in this cell

        # degenerate (e.g. dz_cut == 0): nothing can pass, fall back is safe
        if not (bin_width > 0.0):
            logger.warning(f"MD binned merge: degenerate bin width, falling back to full merge")
            mds = pd.merge(lower, upper, on=self.md_cols, how="inner",
                                suffixes=("_lower", "_upper"))
            return mds, n_full

        # discretise z (shared grid origin) and predict each lower hit's bin
        z0 = z_up.min()
        upper_bin = np.floor((z_up - z0) / bin_width).astype(np.int64)
        lower_bin = np.floor((center - z0) / bin_width).astype(np.int64)

        # binned inner merge over {bin-1, bin, bin+1}: replicate the lower side
        # across the three neighbour keys, merge on the bin key. Each upper hit
        # lives in exactly one bin, so it matches at most one lower copy -> no dupes.
        lower_keys = pd.DataFrame({
            "_zbin": np.concatenate([lower_bin - 1, lower_bin, lower_bin + 1]),
            "_lpos": np.tile(np.arange(n_lo), 3),
        })
        upper_keys = pd.DataFrame({
            "_zbin": upper_bin,
            "_upos": np.arange(n_up),
        })
        matched = lower_keys.merge(upper_keys, on="_zbin", how="inner")

        if len(matched) == 0:
            logger.warning(f"MD binned merge: no matches, falling back to full merge")
            mds = pd.merge(lower.iloc[:0], upper.iloc[:0], on=self.md_cols,
                                how="inner", suffixes=("_lower", "_upper"))
            return mds, n_full

        lower_pos = matched["_lpos"].to_numpy()
        upper_pos = matched["_upos"].to_numpy()

        # assemble the merge-equivalent frame (keys once, everything else suffixed)
        L = lower.iloc[lower_pos].reset_index(drop=True)
        U = upper.iloc[upper_pos].reset_index(drop=True)
        keys = L[self.md_cols]
        L_other = L.drop(columns=self.md_cols).add_suffix("_lower")
        U_other = U.drop(columns=self.md_cols).add_suffix("_upper")
        mds = pd.concat([keys, L_other, U_other], axis=1)

        return mds, n_full
