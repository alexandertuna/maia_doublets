import argparse
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
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
    # "axes.axisbelow": True,
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


def main():
    logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")
    args = options()
    check_options(args)
    plotter = OccupancyPlotter(get_input_files(args.input), args.output, args.n)
    plotter.plot()


def options():
    parser = argparse.ArgumentParser(description="Demo for 3D scatter plot")
    parser.add_argument("--input", type=str, help="Comma-separated list of input file paths")
    parser.add_argument("--output", type=str, help="Output pdf file path")
    parser.add_argument("-n", type=int, help="Sample n random points to speed things up")
    return parser.parse_args()


def check_options(args):
    if not args.input:
        raise ValueError("Input file(s) must be specified")
    if not args.output:
        raise ValueError("Output file must be specified")


def get_input_files(input_str):
    fnames = []
    for pattern in input_str.split(","):
        fnames.extend(glob(pattern))
    return fnames


GAP = 2.0
RADIUS = {
    ("v01", "IT"): [127, 127+GAP, 167, 167+GAP, 510, 510+GAP, 550, 550+GAP],
    ("v01", "OT"): [819, 819+GAP, 899, 899+GAP, 1366, 1366+GAP, 1446, 1446+GAP],
    ("v05", "IT"): [127, 127+GAP, 268, 268+GAP, 409, 409+GAP, 550, 550+GAP],
    ("v05", "OT"): [819, 819+GAP, 1028, 1028+GAP, 1237, 1237+GAP, 1446, 1446+GAP],
} # mm


class OccupancyPlotter:

    def __init__(self,
                 input_files_hits: list[str],
                 output_file: str,
                 n_random: int,
                 ):
        self.input_files_hits = input_files_hits
        self.output_file = output_file
        self.keys = ["simhit_x", "simhit_y", "simhit_z"]
        self.hits_cols = {
            "simhit_x": "x",
            "simhit_y": "y",
            "simhit_z": "z",
        }
        self.cargs = {
            # "cmap": "gist_rainbow",
            # "cmap": "nipy_spectral",
            # "cmap": "rainbow",
            "cmap": "jet",
            "cmin": 0.5,
            "norm": colors.LogNorm(vmin=5e3, vmax=1.5e5),
        }
        bins_r = [
            0,
            127-2*GAP, 127+4*GAP, 167-2*GAP, 167+4*GAP, 510-2*GAP, 510+4*GAP, 550-2*GAP, 550+4*GAP,
            819-2*GAP, 819+4*GAP, 899-2*GAP, 899+4*GAP, 1366-2*GAP, 1366+4*GAP, 1446-2*GAP, 1446+4*GAP,
            1550,
        ]
        bins_z = np.concatenate(([-1550], np.linspace(-1264.2, 1264.2, 250), [1550]))
        self.bins_rz = [bins_z, bins_r]
        self.bins_xy = [np.linspace(-1600, 1600, 400), np.linspace(-1600, 1600, 400)]
        self.get_input_data()


    def get_input_data(self):
        for fi in self.input_files_hits:
            logger.info(f"Reading input file: {fi}")
        self.n_events = len(self.input_files_hits)
        self.hits_df = pd.concat([pd.read_pickle(fi)[self.keys] for fi in self.input_files_hits], ignore_index=True)
        self.hits_df = self.hits_df.rename(columns=self.hits_cols)
        self.hits_df["w"] = 1 / self.n_events
        logger.info(f"Total hits read: {len(self.hits_df)}")
        logger.info(f"Total events read: {self.n_events}")
        logger.info(f"Total hits per event: {len(self.hits_df)/self.n_events:.2f}")


    def plot(self):
        with PdfPages(self.output_file) as pdf:
            # self.plot_xy_occupancy(pdf)
            self.plot_rz_occupancy(pdf)


    def plot_rz_occupancy(self, pdf):
        fig, ax = plt.subplots()
        _, _, _, im = ax.hist2d(self.hits_df["z"],
                                np.sqrt(self.hits_df["x"]**2 + self.hits_df["y"]**2),
                                bins=self.bins_rz,
                                weights=self.hits_df["w"],
                                **self.cargs,
                                )
        fig.colorbar(im, ax=ax, pad=0.01, label="Number of hits")
        ax.set_xlabel("z [mm]")
        ax.set_ylabel("r [mm]")
        ax.set_title(f"Average hits per event: {len(self.hits_df)/self.n_events:.1e}")
        pdf.savefig()
        plt.close()


    def plot_xy_occupancy(self, pdf):
        fig, ax = plt.subplots()
        _, _, _, im = ax.hist2d(self.hits_df["x"],
                                self.hits_df["y"],
                                bins=self.bins_xy,
                                weights=self.hits_df["w"],
                                **self.cargs,
                                )
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        fig.colorbar(im, ax=ax, pad=0.01, label="Number of hits")
        pdf.savefig()
        plt.close()


if __name__ == "__main__":
    main()
