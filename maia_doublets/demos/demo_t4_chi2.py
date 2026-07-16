import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams

FNAME_V01 = "v01_signal_digi_10um/t4s.pkl"
FNAME_V05 = "v05_signal_digi_10um/t4s.pkl"
N_LAYERS = 8
LAYERS = list(range(N_LAYERS))
DET = {
    "IT": 0,
    "OT": 4,
}
XY_COLS = [f"t4_chi2_xy_{it}" for it in LAYERS] + ["t4_chi2_xy"]
SZ_COLS = [f"t4_chi2_sz_{it}" for it in LAYERS] + ["t4_chi2_sz"]
COLOR = {
    "v01": "blue",
    "v05": "red",
}
GAP = 2.0 # mm
RADIUS = {
    ("v01", "IT"): [127, 127+GAP, 167, 167+GAP, 510, 510+GAP, 550, 550+GAP],
    ("v01", "OT"): [819, 819+GAP, 899, 899+GAP, 1366, 1366+GAP, 1446, 1446+GAP],
    ("v05", "IT"): [127, 127+GAP, 268, 268+GAP, 409, 409+GAP, 550, 550+GAP],
    ("v05", "OT"): [819, 819+GAP, 1028, 1028+GAP, 1237, 1237+GAP, 1446, 1446+GAP],
} # mm
EPSILON = 1e-7

def main():
    df = get_dataframes()
    with PdfPages("demo_t4_chi2.pdf") as pdf:
        # plot_pt(df, pdf)
        plot_chi2_cols(df, pdf)


def get_dataframes():
    df = {}
    df["v01"] = pd.read_pickle(FNAME_V01)
    df["v05"] = pd.read_pickle(FNAME_V05)
    for version in df:
        mask = (
            df[version]["t4_ok_mcp"] & 
            df[version]["t4_ok_first_exit"] & 
            df[version]["t4_ok_from_fiducial_mcp"]
        )
        df[version] = df[version][mask]
    return df


def plot_pt(df, pdf):
    args = {
        "fontsize": 14,
        "ha": "right",
        "va": "top",
    }
    # IT, OT
    for name, gdl in DET.items():
        fig, ax = plt.subplots()
        for i_version, version in enumerate(df):
            mask = df[version]["t4_gdoublelayer"] == gdl
            series = df[version][mask]["t4_pt"]
            color = COLOR[version]
            median = series.median()
            perc997 = np.abs(series).quantile(0.997)
            label = f"{version}, med={median:.3f}, 99.7%={perc997:.3f}"
            ax.hist(series, bins=100, alpha=0.5, color=color, histtype="stepfilled", edgecolor="black")
            ax.text(0.98, 0.98 - 0.05 * i_version, label, color=color, transform=ax.transAxes, **args)
        ax.set_xlabel("pT [GeV]")
        ax.set_ylabel("T4s")
        ax.set_title(f"{name} T4s")
        pdf.savefig()
        plt.close()


def plot_chi2_cols(df, pdf):

    bins = np.logspace(-6, 2, 201)
    histtype = "stepfilled"
    edgecolor = "black"

    for cols in [
        XY_COLS,
        SZ_COLS,
    ]:

        # IT, OT
        for name, gdl in DET.items():

            medians, percs = {}, {}

            # column
            for col in cols:

                fig, ax = plt.subplots()

                # v01, v05
                for i_version, version in enumerate(df):
                    color = COLOR[version]
                    mask = df[version]["t4_gdoublelayer"] == gdl
                    series = df[version][mask][col]
                    median = series.median()
                    perc997 = np.abs(series).quantile(0.997)
                    medians[version, col] = median
                    percs[version, col] = perc997
                    label = f"{version}, med={median:.2e}, 99.7%={perc997:.2e}"
                    ax.hist(series, bins=bins, alpha=0.5, color=color, histtype=histtype, edgecolor=edgecolor)
                    ax.text(0.98, 0.98 - 0.05 * i_version, label, color=color, transform=ax.transAxes, fontsize=14, ha="right", va="top")

                print(f"Plotting {name}, GDL {gdl}, {col} ... ")
                ax.set_title(f"{name}, GDL {gdl}, {col}")
                ax.set_xlabel(col)
                ax.set_ylabel("Entries")
                ax.semilogx()
                ax.semilogy()
                y_lo, y_hi = ax.get_ylim()
                ax.set_ylim(y_lo, y_hi * 2)
                pdf.savefig()
                plt.close()

            # table of 99.7% quantiles
            fig, ax = plt.subplots()
            ax.axis("off")
            x, y = 0.1, 0.9
            cell_h, cell_w = 0.06, 0.20
            header = ["", "v01 99.7%", "v05 99.7%", "v01 R [mm]", "v05 R [mm]"]
            n_rows = len(cols) + 1
            n_cols = len(header)
            for i_row in range(n_rows):
                if i_row == 0:
                    row = header
                else:
                    col = cols[i_row - 1]
                    row = [col] + [f"{percs[version, col]:0.5f}" for version in df]
                    try:
                        index = int(col[-1])
                        row.append(f"{RADIUS['v01', name][index]:.0f}")
                        row.append(f"{RADIUS['v05', name][index]:.0f}")
                    except:
                        row.extend(["", ""])
                for i_col in range(n_cols):
                    text = row[i_col]
                    ax.text(x + i_col * cell_w, y - i_row * cell_h, text, fontsize=14, ha="right", va="top")
                print(", ".join(row))
            pdf.savefig()
            plt.close()


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
