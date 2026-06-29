import numpy as np
import pandas as pd
import time
import logging
logger = logging.getLogger(__name__)

from maia_doublets.constants import MD_DZ_CUT, MD_DR_CUT
from maia_doublets.constants import MAGNETIC_FIELD, SPEED_OF_LIGHT
from maia_doublets.constants import BYTE_TO_MB, MEV_TO_GEV, NO_MCP
from maia_doublets.constants import N_T2_PHI_SLICES, N_T2_ETA_SLICES, DETECTOR_MAX_ETA, DETECTOR_MAX_PHI

class MDMaker:

    def __init__(
        self,
        geometry_version: str,
        sim: bool,
        smear: str,
        signal: bool,
        cut_doublets: bool,
        fast_merge: bool,
        simhits: pd.DataFrame,
    ):
        self.signal = signal
        self.cut_doublets = cut_doublets
        key = (geometry_version, "sim") if sim else (geometry_version, "digi", smear)
        self.MD_DZ_CUT = MD_DZ_CUT[key]
        self.MD_DR_CUT = MD_DR_CUT[key]
        self.doublet_cols = [
            "file",
            "i_event", # the event
            "simhit_system", # the system (IT, OT)
            "simhit_layer_div_2", # the double layer
            "simhit_module", # the phi-module
            "simhit_sensor", # the z-sensor
        ]
        self.fast_merge = fast_merge
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
                "simhit_sensor", # the z-sensor
            ]

        # group loop
        logger.info("Grouping simhits ...")
        groups = df.groupby(groupby_cols)
        all_doublets, all_cutflows = [], []

        for i_group, (cols, group) in enumerate(groups):

            doublets, cutflow = self.make_doublets_from_group(group)

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


    def make_doublets_from_group(self, group: pd.DataFrame) -> tuple[pd.DataFrame, dict]:

        lower_mask = group["simhit_layer_mod_2"] == 0
        upper_mask = group["simhit_layer_mod_2"] == 1
        lower = group[lower_mask]
        upper = group[upper_mask]

        # inner join to find doublets
        if self.cut_doublets and self.fast_merge:
            doublets, n_full = self.merge_binned(lower, upper)
        else:
            doublets = pd.merge(
                lower,
                upper,
                on=self.doublet_cols,
                how="inner",
                suffixes=("_lower", "_upper"),
            )
            n_full = len(doublets)

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
        sy = doublets["simhit_system"]
        dl = doublets["simhit_layer_div_2"]
        mask["dr"] = np.abs(doublets["doublet_dr"]) < self.MD_DR_CUT[sy, dl]
        mask["dz"] = np.abs(doublets["doublet_dz"]) < self.MD_DZ_CUT[sy, dl]
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
        doublets["doublet_glayer"] = doublets["simhit_glayer_lower"]
        doublets["doublet_gdoublelayer"] = doublets["simhit_glayer_lower"] // 2

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

        # divide the eta/phi space into slices, to be used in T2 seeding
        n_phi_slices = N_T2_PHI_SLICES[doublets["doublet_system"]]
        n_eta_slices = N_T2_ETA_SLICES[doublets["doublet_system"]]
        doublets["doublet_phi_slice"] = np.floor((doublets["doublet_phi"] + DETECTOR_MAX_PHI) / (2 * DETECTOR_MAX_PHI) * n_phi_slices).astype(np.int16)
        doublets["doublet_eta_slice"] = np.floor((doublets["doublet_eta"] + DETECTOR_MAX_ETA) / (2 * DETECTOR_MAX_ETA) * n_eta_slices).astype(np.int16)

        # guess charge from dphi:
        # positively charged particles have negative dphi, and vice versa
        doublets["doublet_q"] = (-1*np.sign(doublets["doublet_dphi"])).astype(np.int8)

        # pass-through the simhit positions
        for coord in ["x", "y", "r", "z"]:
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
            doublets["doublet_from_fiducial_mcp"] = doublets["simhit_from_fiducial_mcp_lower"] & doublets["simhit_from_fiducial_mcp_upper"]
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


    def merge_binned(self, lower: pd.DataFrame, upper: pd.DataFrame) -> tuple[pd.DataFrame, int]:
        """
        Experimental binned merge to speed up the doublet making, courtesy of Claude

        Binned-merge equivalent of pd.merge(lower, upper, on=doublet_cols,
        how="inner", suffixes=("_lower","_upper")), restricted to upper hits whose
        z bin is within +/-1 of the lower hit's predicted z bin. Returns
        (merge_equivalent_frame, n_full_crossproduct).

        Assumes a single cell (system / doublelayer / module / sensor constant).

        Bin width: within this cell the dz cut |doublet_dz| < DZ is exactly
        |z_lo*r_up - z_up*r_lo| < DZ*(r_up-r_lo), which confines z_up to an interval
        around the radial projection of the lower hit. The bin width is the widest
        such half-interval the cut allows in the cell, so the predicted bin plus its
        two neighbours are guaranteed to contain every dz survivor.
        """
        n_lo = len(lower)
        n_up = len(upper)
        n_full = n_lo * n_up

        system = int(lower["simhit_system"].iloc[0])
        doublelayer = int(lower["simhit_layer_div_2"].iloc[0])
        dz_cut = float(self.MD_DZ_CUT[system, doublelayer])

        z_lo = lower["simhit_z"].to_numpy(np.float64)
        r_lo = lower["simhit_r"].to_numpy(np.float64)
        z_up = upper["simhit_z"].to_numpy(np.float64)
        r_up = upper["simhit_r"].to_numpy(np.float64)

        r_up_min = r_up.min()
        r_up_max = r_up.max()

        # guard: closed form assumes the upper layer sits strictly outside the
        # lower one. If radii overlap, fall back to the exact full merge.
        if r_up_min <= r_lo.max():
            self.logger.warning(f"MD binned merge: weird data, falling back to full merge")
            doublets = pd.merge(lower, upper, on=self.doublet_cols, how="inner",
                                suffixes=("_lower", "_upper"))
            return doublets, n_full

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
            self.logger.warning(f"MD binned merge: degenerate bin width, falling back to full merge")
            doublets = pd.merge(lower, upper, on=self.doublet_cols, how="inner",
                                suffixes=("_lower", "_upper"))
            return doublets, n_full

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
            self.logger.warning(f"MD binned merge: no matches, falling back to full merge")
            doublets = pd.merge(lower.iloc[:0], upper.iloc[:0], on=self.doublet_cols,
                                how="inner", suffixes=("_lower", "_upper"))
            return doublets, n_full

        lower_pos = matched["_lpos"].to_numpy()
        upper_pos = matched["_upos"].to_numpy()

        # assemble the merge-equivalent frame (keys once, everything else suffixed)
        L = lower.iloc[lower_pos].reset_index(drop=True)
        U = upper.iloc[upper_pos].reset_index(drop=True)
        keys = L[self.doublet_cols]
        L_other = L.drop(columns=self.doublet_cols).add_suffix("_lower")
        U_other = U.drop(columns=self.doublet_cols).add_suffix("_upper")
        doublets = pd.concat([keys, L_other, U_other], axis=1)

        return doublets, n_full
