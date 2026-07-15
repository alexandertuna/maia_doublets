"""
A script to plot the calibration data for maia doublets.
For example, find the keys of:
    calibs/v01_digi_10um.json
    calibs/v05_digi_10um.json
And overlay the plots for comparison.

The data looks like:
> cat calibs/v01_digi_10um.json  | head -n 20
{
    "doublet_dz": {
        "3": {
            "0": 3.0,
            "1": 4.0,
            "2": 13.0,
            "3": 15.0
        },
        "5": {
            "0": 31.0,
            "1": 38.0,
            "2": 120.0,
            "3": 145.0
        }
    },
    "doublet_dr": {
        "3": {
            "0": 9.0,
            "1": 14.0,
            "2": 108.0,
"""

import argparse
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams

FNAMES = {
    # ("v01", "sim"): "calibs/v01_sim.json",
    ("v01", "10um"): "calibs/v01_digi_10um.json",
    ("v05", "10um"): "calibs/v05_digi_10um.json",
}
KEYS = [
    "doublet_dz",
    "doublet_dr",
    "ls_dz",
    "ls_dr",
    "ls_chi2_xy",
    "ls_chi2_sz",
    "t4_dz",
    "t4_dr",
    "t4_chi2_xy",
    "t4_chi2_sz",
    "t8_dz",
    "t8_dr",
    "t8_chi2_xy",
    "t8_chi2_sz",
]


def main():
    with PdfPages("calibs.pdf") as pdf:
        plot_calibs(pdf)
        write_calibs(pdf)


def plot_calibs(pdf):
    all_global_layers = []
    x_row, y_row, row_sep = 0.55, 0.4, 0.05
    row_len = 400

    # get the calib data
    data = {}
    for fname in FNAMES.values():
        with open(fname, "r") as f:
            data[fname] = json.load(f)

    # make 1 plot for each key (quantity)
    for key in KEYS:
        lo, hi = marker_lo_hi(key)
        fig, ax = plt.subplots()
        all_values = []

        # overlay the plots for each detector and smearing
        for i_det, ((detector, smear), fname) in enumerate(FNAMES.items()):
            calib = data[fname][key]
            global_layers, calib_values = interpret_calib_data(calib, lo, hi, row_len)
            ax.scatter(global_layers, calib_values, marker="s", s=1, color=color(detector, smear))
            ax.text(x_row, y_row - i_det * row_sep, name(detector, smear), transform=ax.transAxes, color=color(detector, smear))
            all_global_layers.extend(global_layers)
            all_values.extend(calib_values)

        # set the bin labels
        unique_global_layers = np.unique(np.array(all_global_layers).astype(int))
        if not np.array_equal(unique_global_layers, np.arange(len(unique_global_layers))):
            raise ValueError(f"Unexpected global layers: {unique_global_layers}")
        ax.set_xticks(range(len(unique_global_layers)))
        ax.set_xticklabels([bin_name(gl) for gl in unique_global_layers])
        ax.set_xlabel("Double-layer (global)")
        ax.set_ylabel(f"3-sigma interval {unit(key)}")
        ax.set_title(nickname(key))
        if semilogy(key):
            ax.semilogy()
        else:
            top = max(all_values) * 1.1
            ax.set_ylim(bottom=0, top=top)
        pdf.savefig()
        plt.close()


def write_calibs(pdf):
    lo = hi = 0.0
    row_len = 1

    # get the calib data
    data = {}
    for fname in FNAMES.values():
        with open(fname, "r") as f:
            data[fname] = json.load(f)

    # write one page of calibs for each detector
    for ((detector, smear), fname) in FNAMES.items():

        text = f"Intervals for {name(detector, smear)} ({detector})\n"
        text += f"{'Feature':15s} {'First layer':15s} {'Value':12s}\n"

        # make one line for each feature
        for key in KEYS:

            calib = data[fname][key]
            global_doublelayers, calib_values = interpret_calib_data(calib, lo, hi, row_len)
            key = key.replace("doublet_", "md_").replace("ls_", "t2_")
            for (gdl, value) in zip(global_doublelayers, calib_values):
                gdl = int(gdl)
                text += f"{key:15s} {bin_name(gdl):15s} {value:12.6f}\n"

        # plot the text on a page
        fig, ax = plt.subplots(figsize=(8, 8))
        args = {"ha":"left", "va":"top", "fontfamily":"monospace", "fontsize":10}
        ax.text(0.0, 1.0, text, **args)
        ax.axis("off")
        pdf.savefig()
        plt.close()


def name(detector: str, smear: str) -> str:
    if detector == "v01":
        if smear == "sim":
            return "Super-layers (sim)"
        else:
            return r"Super-layers, $\sigma$=10um"
    elif detector == "v05":
        return r"Equal-spaced, $\sigma$=10um"
    else:
        raise ValueError(f"Unknown detector: {detector}")


def semilogy(feature: str) -> bool:
    return (
        feature.startswith("doublet_") or
        feature.startswith("md_") or
        feature.startswith("ls_") or
        feature.startswith("t2_")
    ) and (
        "dz" in feature or
        "dr" in feature
    )

def nickname(feature: str) -> str:
    if feature == "doublet_dz":
        return r"Doublet $\Delta z$"
    elif feature == "doublet_dr":
        return r"Doublet $\Delta r$"
    elif feature == "ls_dz":
        return r"T2 $\Delta z$"
    elif feature == "ls_dr":
        return r"T2 $\Delta r$"
    elif feature == "ls_chi2_xy":
        return r"T2 $\chi^2_{xy}$"
    elif feature == "ls_chi2_sz":
        return r"T2 $\chi^2_{sz}$"
    elif feature == "t4_dz":
        return r"T4 $\Delta z$"
    elif feature == "t4_dr":
        return r"T4 $\Delta r$"
    elif feature == "t4_chi2_xy":
        return r"T4 $\chi^2_{xy}$"
    elif feature == "t4_chi2_sz":
        return r"T4 $\chi^2_{sz}$"
    elif feature == "t8_dz":
        return r"T8 $\Delta z$"
    elif feature == "t8_dr":
        return r"T8 $\Delta r$"
    elif feature == "t8_chi2_xy":
        return r"T8 $\chi^2_{xy}$"
    elif feature == "t8_chi2_sz":
        return r"T8 $\chi^2_{sz}$"
    else:
        raise ValueError(f"Unknown feature: {feature}")


def unit(feature: str) -> str:
    if feature in [
        "doublet_dz",
        "doublet_dr",
        "ls_dz",
        "ls_dr",
        "t4_dz",
        "t4_dr",
        "t8_dz",
        "t8_dr",
    ]:
        return "[mm]"
    elif feature in [
        "ls_chi2_xy",
        "ls_chi2_sz",
        "t4_chi2_xy",
        "t4_chi2_sz",
        "t8_chi2_xy",
        "t8_chi2_sz",
    ]:
        return "[mm^2]"
    else:
        raise ValueError(f"Unknown feature: {feature}")


def marker_lo_hi(feature: str) -> tuple[float, float]:
    if feature.startswith("doublet_") or feature.startswith("md_"):
        return 0.5 - 0.1, 0.5 - 0.1
    elif feature.startswith("ls_"):
        return 0.5 - 0.1, 1.5 - 0.1
    elif feature.startswith("t4_"):
        return 0.5 - 0.1, 3.5 - 0.1
    elif feature.startswith("t8_"):
        return 0.5 - 0.1, 7.5 - 0.1
    else:
        raise ValueError(f"Unknown feature: {feature}")


def bin_name(global_layer: int) -> str:
    if global_layer == 0:
        return "IT01"
    elif global_layer == 1:
        return "IT23"
    elif global_layer == 2:
        return "IT45"
    elif global_layer == 3:
        return "IT67"
    elif global_layer == 4:
        return "OT01"
    elif global_layer == 5:
        return "OT23"
    elif global_layer == 6:
        return "OT45"
    elif global_layer == 7:
        return "OT67"
    else:
        raise ValueError(f"Unknown global layer: {global_layer}")


def color(detector: str, smear: str) -> str:
    if detector == "v01":
        if smear == "sim":
            return "black"
        else:
            return "blue"
    elif detector == "v05":
        return "green"
    else:
        raise ValueError(f"Unknown detector: {detector}")


def convert_system_layer_to_global(system, layer):
    offsets = {
        3: 0,
        5: 4,
    }
    return offsets[system] + layer


def interpret_calib_data(calib: dict, lo: float, hi: float, row_len: int) -> tuple[list, list]:
    global_doublelayers = []
    calib_values = []
    for key, value in calib.items():
        if (isinstance(value, float) or
            isinstance(value, int)):
            value = float(value)
            global_doublelayer = int(key)
            for gdl in np.linspace(global_doublelayer - lo,
                                  global_doublelayer + hi,
                                  row_len):
                global_doublelayers.append(gdl)
                calib_values.append(value)

        elif isinstance(value, dict):
            system = int(key)
            for doublelayer_str, calib_value in value.items():
                doublelayer = int(doublelayer_str)
                global_doublelayer = convert_system_layer_to_global(system, doublelayer)
                for gdl in np.linspace(global_doublelayer - lo,
                                       global_doublelayer + hi,
                                       row_len):
                    global_doublelayers.append(gdl)
                    calib_values.append(calib_value)
        else:
            raise ValueError(f"Unknown type for value: {type(value)}")

    return global_doublelayers, calib_values


rcParams.update({
    "font.size": 16,
    "figure.figsize": (8, 8),
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "axes.grid": True,
    "axes.grid.which": "both",
    # "axes.axisbelow": True,
    "grid.linewidth": 0.5,
    "grid.alpha": 0.1,
    "grid.color": "gray",
    "figure.subplot.left": 0.15,
    "figure.subplot.bottom": 0.09,
    "figure.subplot.right": 0.97,
    "figure.subplot.top": 0.95,
})

if __name__ == "__main__":
    main()
