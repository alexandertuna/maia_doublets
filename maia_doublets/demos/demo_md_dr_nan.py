from glob import glob
import pandas as pd
import logging
logger = logging.getLogger(__name__)

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams

PKL_DIR = "/ceph/users/atuna/work/maia/maia_doublets/run/debug"

def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    df = get_df()
    with PdfPages("demo_md_dr_nan.pdf") as pdf:
        plot_all(df, pdf)
        plot_event_displays(df, pdf)


def get_df():
    logger.info(f"Reading pkl files from {PKL_DIR} ... ")
    pkl_files = sorted(glob(f"{PKL_DIR}/*.pkl"))

    logger.info(f"Concat {len(pkl_files)} pkl files ... ")
    df = pd.concat([pd.read_pickle(f) for f in pkl_files], ignore_index=True)
    return df


def plot_all(df, pdf):
    size = 1
    logger.info(f"Plotting all ... ")
    fig, ax = plt.subplots()
    ax.scatter(df['simhit_x_lower'], df['simhit_y_lower'], s=size, color="blue", label="lower")
    ax.scatter(df['simhit_x_upper'], df['simhit_y_upper'], s=size, color="red", label="upper")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    pdf.savefig()
    plt.close() 
   

def plot_event_displays(df, pdf):
    colors = {
        "lower": "blue",
        "upper": "red",
    }
    n_row = len(df)
    r_edge = max(
        abs(df['simhit_x_lower'].max()),
        abs(df['simhit_x_lower'].min()),
        abs(df['simhit_x_upper'].max()),
        abs(df['simhit_x_upper'].min()),
        abs(df['simhit_y_lower'].max()),
        abs(df['simhit_y_lower'].min()),
        abs(df['simhit_y_upper'].max()),
        abs(df['simhit_y_upper'].min()),
    )

    # One page per row
    for i_row, row in df.iterrows():

        if i_row % 10 == 0:
            logger.info(f"Plotting event display for row {i_row+1}/{n_row} ... ")
        if i_row > 50:
            break

        fig, ax = plt.subplots()
        x = {
            "lower": row['simhit_x_lower'],
            "upper": row['simhit_x_upper'],
        }
        y = {
            "lower": row['simhit_y_lower'],
            "upper": row['simhit_y_upper'],
        }
        for key in x:
            ax.scatter(x[key], y[key], color=colors[key])
        ax.set_xlim(-r_edge, r_edge)
        ax.set_ylim(-r_edge, r_edge)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
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
