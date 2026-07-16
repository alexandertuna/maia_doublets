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
XY_COLS = ["t4_chi2_xy"] + [f"t4_chi2_xy_{it}" for it in LAYERS]
SZ_COLS = ["t4_chi2_sz"] + [f"t4_chi2_sz_{it}" for it in LAYERS]
COLOR = {
    "v01": "blue",
    "v05": "red",
}
EPSILON = 1e-7

def main():
    df = get_dataframes()
    with PdfPages("demo_t4_chi2.pdf") as pdf:
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


def plot_chi2_cols(df, pdf):

    bins = np.logspace(-6, 2, 201)

    for cols in [
        XY_COLS,
        # SZ_COLS,
    ]:

        # IT, OT
        for name, gdl in DET.items():

            # column
            for col in cols:

                fig, ax = plt.subplots()

                # v01, v05
                for i_version, version in enumerate(df):
                    color = COLOR[version]
                    mask = df[version]["t4_gdoublelayer"] == gdl
                    series = df[version][mask][col]
                    label = f"{version}, median={series.median():.2e}"
                    ax.hist(series, bins=bins, alpha=0.5, color=color, histtype="stepfilled")
                    ax.text(0.95, 0.95 - 0.05 * i_version, label, transform=ax.transAxes, ha="right", va="top")

                print(f"Plotting {name}, GDL {gdl}, {col} ... ")
                ax.set_title(f"{name}, GDL {gdl}, {col}")
                ax.set_xlabel(col)
                ax.set_ylabel("Entries")
                ax.semilogx()
                ax.semilogy()
                pdf.savefig()
                plt.close()


def plot_chi2_sz(df, pdf):
    pass


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
