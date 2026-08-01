import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
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
    "grid.linewidth": 0.5,
    "grid.alpha": 0.1,
    "grid.color": "gray",
    "figure.subplot.left": 0.15,
    "figure.subplot.bottom": 0.09,
    "figure.subplot.right": 0.97,
    "figure.subplot.top": 0.95,
})

import logging
logger = logging.getLogger(__name__)

SMEARS = ["00um", "05um", "10um", "20um"]
COLORS = {
    "00um": "white",
    "05um": "red",
    "10um": "blue",
    "20um": "green",
}

def main():
    logging.basicConfig(level=logging.INFO)
    with PdfPages("demo_chi2s.pdf") as pdf:
        plot_chi2_xy(pdf)


def plot_chi2_xy(pdf: PdfPages) -> None:
    data = {smear: get_data(smear) for smear in SMEARS}
    bins = np.logspace(-8, 0, 100)
    fig, ax = plt.subplots()
    for smear, series in data.items():
        ax.hist(
            series["t2_chi2_xy"],
            bins=bins,
            histtype="stepfilled",
            label=f"Smear={smear}",
            color=COLORS[smear],
            edgecolor="black",
            alpha=0.3,
        )
    ax.semilogx()
    ax.semilogy()
    ax.set_xlabel("Diff^2 between circle(xy, 012) and 3 [mm2]")
    ax.set_ylabel("T2s")
    ax.legend()
    pdf.savefig()
    plt.close()



def get_data(smear: str) -> pd.Series:
    df = pd.read_pickle(f"v01_signal_digi_{smear}/t2s.pkl")
    mask = df["i_mcp"] != 0xffff_ffff
    series = df[mask][["t2_chi2_xy"]]
    return series


if __name__ == "__main__":
    main()
