import argparse
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import logging
logger = logging.getLogger(__name__)

N_ITER = 1
MAX_RADIUS_L0 = 127.5
MIN_RADIUS_L1 = 128.5

def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    ops = options()
    group_df = pd.read_pickle(ops.group)
    results_df = pd.read_pickle(ops.results)

    weird_hit_mask = (
        (group_df["simhit_x"] > MAX_RADIUS_L0) &
        (group_df["simhit_x"] < MIN_RADIUS_L1)
    )
    group_df = group_df[~weird_hit_mask]

    for binned in [
        False,
        True,
    ]:
        for it in range(N_ITER):
            with Timer() as timer:
                doublets, cutflow = make_doublets_from_group(group_df, binned)
            if it == 0:
                print_comparison(results_df, doublets)
            logger.info(f"T(binned={binned}) = {timer.duration:.2f}")


def options():
    parser = argparse.ArgumentParser(description="Compare algorithms for making MDs from a group dataframe with two layers")
    parser.add_argument("--group", type=str, required=True, help="Group dataframe (pkl) with two layers to be converted into mds")
    parser.add_argument("--results", type=str, required=True, help="Results dataframe (pkl) of the expected output mds")
    return parser.parse_args()


def postprocess(df: pd.DataFrame):
    columns = [
        "doublet_eta",
        "doublet_phi",
    ]
    df = df.reset_index(drop=True)[sorted(df.columns)]
    df = df.sort_values(by=columns).reset_index(drop=True)
    return df


def print_comparison(df1, df2):
    df1 = postprocess(df1)
    df2 = postprocess(df2)
    try:
        assert_frame_equal(df1, df2, check_dtype=False, check_column_type=False)
        logger.info("Dataframes are equal according to assert_frame_equal")
    except AssertionError as e:
        logger.error(f"Dataframes are not equal: {e}")


def make_doublets_from_group(
    group: pd.DataFrame,
    use_binned: bool,
) -> tuple[pd.DataFrame, dict]:

    doublet_cols = [
        "file",
        "i_event", # the event
        "simhit_system", # the system (IT, OT)
        "simhit_layer_div_2", # the double layer
        "simhit_module", # the phi-module
        "simhit_sensor", # the z-sensor
    ]
    lower_mask = group["simhit_layer_mod_2"] == 0
    upper_mask = group["simhit_layer_mod_2"] == 1
    lower = group[lower_mask]
    upper = group[upper_mask]

    if use_binned:
        doublets, n_full = merge_binned(lower, upper)
    else:
        doublets = pd.merge(
            lower,
            upper,
            on=doublet_cols,
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
    mask["dr"] = np.abs(doublets["doublet_dr"]) < MD_DR_CUT[sy, dl]
    mask["dz"] = np.abs(doublets["doublet_dz"]) < MD_DZ_CUT[sy, dl]
    mask["and"] = mask["dr"] & mask["dz"]
    doublets["doublet_ok"] = mask["and"].astype(bool)

    # remove as desired
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

    # drop columns which arent used downstream
    dropcols = ["i_mcp_lower", "i_mcp_upper"]
    dropcols.extend([col for col in doublets.columns if col.startswith("simhit_")])
    dropcols.extend([col for col in doublets.columns if col.startswith("mcp_") and col.endswith("_lower")])
    dropcols.extend([col for col in doublets.columns if col.startswith("mcp_") and col.endswith("_upper")])
    doublets.drop(columns=dropcols, inplace=True)

    return doublets, cutflow


def merge_binned(lower: pd.DataFrame, upper: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """
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
    doublet_cols = [
        "file",
        "i_event", # the event
        "simhit_system", # the system (IT, OT)
        "simhit_layer_div_2", # the double layer
        "simhit_module", # the phi-module
        "simhit_sensor", # the z-sensor
    ]

    n_lo = len(lower)
    n_up = len(upper)
    n_full = n_lo * n_up

    system = int(lower["simhit_system"].iloc[0])
    doublelayer = int(lower["simhit_layer_div_2"].iloc[0])
    dz_cut = float(MD_DZ_CUT[system, doublelayer])

    z_lo = lower["simhit_z"].to_numpy(np.float64)
    r_lo = lower["simhit_r"].to_numpy(np.float64)
    z_up = upper["simhit_z"].to_numpy(np.float64)
    r_up = upper["simhit_r"].to_numpy(np.float64)

    r_up_min = r_up.min()
    r_up_max = r_up.max()

    # guard: closed form assumes the upper layer sits strictly outside the
    # lower one. If radii overlap, fall back to the exact full merge.
    if r_up_min <= r_lo.max():
        print("f"*30)
        doublets = pd.merge(lower, upper, on=doublet_cols, how="inner",
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
    logger.info(f"Found bin_width={bin_width:.5f}")

    # degenerate (e.g. dz_cut == 0): nothing can pass, fall back is safe
    if not (bin_width > 0.0):
        doublets = pd.merge(lower, upper, on=doublet_cols, how="inner",
                            suffixes=("_lower", "_upper"))
        return doublets, n_full

    # discretise z (shared grid origin) and predict each lower hit's bin
    z0 = z_up.min()
    upper_bin = np.floor((z_up - z0) / bin_width).astype(np.int64)
    lower_bin = np.floor((center - z0) / bin_width).astype(np.int64)
    logger.info(f"Found min lower_bin={lower_bin.min()}")
    logger.info(f"Found max lower_bin={lower_bin.max()}")
    logger.info(f"Found min upper_bin={upper_bin.min()}")
    logger.info(f"Found max upper_bin={upper_bin.max()}")

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
        doublets = pd.merge(lower.iloc[:0], upper.iloc[:0], on=doublet_cols,
                            how="inner", suffixes=("_lower", "_upper"))
        return doublets, n_full

    lower_pos = matched["_lpos"].to_numpy()
    upper_pos = matched["_upos"].to_numpy()

    # assemble the merge-equivalent frame (keys once, everything else suffixed)
    L = lower.iloc[lower_pos].reset_index(drop=True)
    U = upper.iloc[upper_pos].reset_index(drop=True)
    keys = L[doublet_cols]
    L_other = L.drop(columns=doublet_cols).add_suffix("_lower")
    U_other = U.drop(columns=doublet_cols).add_suffix("_upper")
    doublets = pd.concat([keys, L_other, U_other], axis=1)

    return doublets, n_full


class Timer:
    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, *args):
        self.end = time.perf_counter()
        self.duration = self.end - self.start


SPEED_OF_LIGHT = 299.792458  # mm/ns
MAGNETIC_FIELD = 5.0 # Tesla
NO_MCP = np.uint32(0xffffffff)
MD_DZ_CUT = np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0.62, 1.02, 5.49, 7.37],
        [0, 0, 0, 0],
        [22, 29, 112, 137],
        [0, 0, 0, 0],
    ])
MD_DR_CUT = np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [7.2, 12.1, 102.0, 118.5],
        [0, 0, 0, 0],
        [260, 313, 718, 806],
        [0, 0, 0, 0],
    ])
DETECTOR_MAX_PHI = np.pi
DETECTOR_MAX_ETA = 2.5 # Must include all background hits
MAX_T2_DPHI = np.array([
    2 * DETECTOR_MAX_PHI,
    2 * DETECTOR_MAX_PHI,
    2 * DETECTOR_MAX_PHI,
    0.10 * 2,
    2 * DETECTOR_MAX_PHI,
    0.10 * 2,
    2 * DETECTOR_MAX_PHI,
])
MAX_T2_DETA = np.array([
    2 * DETECTOR_MAX_ETA,
    2 * DETECTOR_MAX_ETA,
    2 * DETECTOR_MAX_ETA,
    0.01 * 2,
    2 * DETECTOR_MAX_ETA,
    0.01 * 2,
    2 * DETECTOR_MAX_ETA,
])
N_T2_PHI_SLICES = (2 * DETECTOR_MAX_PHI / MAX_T2_DPHI).astype(int)
N_T2_ETA_SLICES = (2 * DETECTOR_MAX_ETA / MAX_T2_DETA).astype(int)


if __name__ == "__main__":
    main()
