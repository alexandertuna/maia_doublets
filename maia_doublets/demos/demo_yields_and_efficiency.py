import argparse
from glob import glob
import numpy as np
import pandas as pd
# import datashader as ds
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages
import logging
logger = logging.getLogger(__name__)

MUON = 13
ZERO_POINT_ZERO_ONE_MM = 0.01
ONE_POINT_FIVE_GEV = 1.5
BARREL_TRACKER_MAX_ETA = 0.65

def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    args = options()

    paths = {
        "background_mcps": args.background_mcps,
        "background_hits": args.background_hits,
        "background_mds": args.background_mds,
        "background_t2s": args.background_t2s,
        "background_t4s": args.background_t4s,
        "background_t8s": args.background_t8s,
        "signal_mcps": args.signal_mcps,
        "signal_hits": args.signal_hits,
        "signal_mds": args.signal_mds,
        "signal_t2s": args.signal_t2s,
        "signal_t4s": args.signal_t4s,
        "signal_t8s": args.signal_t8s,
    }

    data = DataGrabber(paths, args.pdf)
    data.plot()


class DataGrabber:

    def __init__(self,
                 data_paths: dict,
                 pdf_path: str,
                 ):
        self.paths = {}
        self.nfiles_per_key = {}
        for key, objs in data_paths.items():
            self.paths[key] = self.get_input_filenames(objs)
            self.nfiles_per_key[key] = len(self.paths[key])
        self.pdf_path = pdf_path
        self.unique_keys = ["file", "i_event", "i_mcp"]
        self.load_data()
        self.calculate_efficiency()
        self.pad_r = 10
        self.color_r = "#1f77b4"
        self.color_l = "lightgrey"
        self.hatch = None


    def get_input_filenames(self, input_str):
        fnames = []
        if isinstance(input_str, list):
            for pattern in input_str:
                for patt in pattern.split(","):
                    fnames.extend(glob(patt))
        else:
            for patt in input_str.split(","):
                fnames.extend(glob(patt))
        return fnames


    def load_data(self):
        self.df = {}
        for key, file_paths in self.paths.items():
            logger.info(f"Loading {key} from {file_paths} ...")
            for fi in file_paths:
                if len(pd.read_pickle(fi)) > 0:
                    break
            else:
                logger.warning(f"No non-empty files found for {key} in {file_paths}. Skipping.")
                self.df[key] = pd.DataFrame(columns=self.unique_keys)
                continue
            keys = self.unique_keys + (["mcp_pdg", "mcp_q", "mcp_pt", "mcp_eta", "mcp_vertex_r", "mcp_vertex_z"] if key == "signal_mcps" else [])
            self.df[key] = pd.concat([pd.read_pickle(fi)[keys] for fi in file_paths], ignore_index=True)
        self.announce_sizes()


    def announce_sizes(self):
        for key in self.df.keys():
            logger.info(f"{key}: {len(self.df[key])} entries")
        for key, nfiles in self.nfiles_per_key.items():
            logger.info(f"{key}: {nfiles} files")


    def denominator_mask(self):
        mask = (
            (np.abs(self.df["signal_mcps"]["mcp_pdg"]) == MUON) &
            (self.df["signal_mcps"]["mcp_q"] != 0) &
            (self.df["signal_mcps"]["mcp_pt"] > ONE_POINT_FIVE_GEV) &
            (self.df["signal_mcps"]["mcp_vertex_r"] < ZERO_POINT_ZERO_ONE_MM) &
            (np.abs(self.df["signal_mcps"]["mcp_vertex_z"]) < ZERO_POINT_ZERO_ONE_MM) &
            (np.abs(self.df["signal_mcps"]["mcp_eta"]) < BARREL_TRACKER_MAX_ETA)
        )
        return mask


    def calculate_efficiency(self):
        mask = self.denominator_mask()
        if mask.sum() == 0:
            msg = "No signal mcps found that satisfy denominator criteria"
            logger.error(msg)
            raise ValueError(msg)
        denom = self.df["signal_mcps"][mask]
        logger.info(f"Denominator mask: {mask.sum()} entries")
        self.efficiency = {}
        for key in self.df.keys():
            if key.startswith("background"):
                continue
            df = self.df[key].drop_duplicates(subset=self.unique_keys)
            # for each row in df, check if it exists in mcps
            merged = pd.merge(df, denom, on=self.unique_keys, how="inner")
            eff = len(merged) / len(denom)
            logger.info(f"Efficiency for {key}: {len(merged)} / {len(denom)} = {eff}")
            self.efficiency[key] = eff


    def plot(self):
        with PdfPages(self.pdf_path) as pdf:
            self.plot_yields_and_efficiency(pdf, do_signal=True)
            self.plot_yields_and_efficiency(pdf, do_signal=False)


    def plot_yields_and_efficiency(self, pdf, do_signal):
        y_vals_l = [
            len(self.df["background_hits"]) / self.nfiles_per_key["background_hits"],
            len(self.df["background_mds"]) / self.nfiles_per_key["background_mds"],
            len(self.df["background_t2s"]) / self.nfiles_per_key["background_t2s"],
            len(self.df["background_t4s"]) / self.nfiles_per_key["background_t4s"],
            len(self.df["background_t8s"]) / self.nfiles_per_key["background_t8s"],
        ]
        if do_signal:
            y_vals_r = [
                self.efficiency["signal_hits"],
                self.efficiency["signal_mds"],
                self.efficiency["signal_t2s"],
                self.efficiency["signal_t4s"],
                self.efficiency["signal_t8s"]
            ]
        else:
            y_vals_r = [0] * len(y_vals_l)
        x_vals = ["Hits", "MDs", "T2s", "T4s", "T8s"]
        bins = np.arange(len(x_vals) + 1) - 0.5

        # check
        if len(y_vals_l) != len(y_vals_r) or len(y_vals_l) != len(x_vals):
            msg = f"Length mismatch: y_vals_l={len(y_vals_l)}, y_vals_r={len(y_vals_r)}, x_vals={len(x_vals)}"
            logger.error(msg)
            raise ValueError(msg)

        # plot background yields
        fig, ax_l = plt.subplots(figsize=(8, 8))
        ax_l.hist(
            x_vals,
            bins=bins,
            weights=y_vals_l,
            histtype="stepfilled",
            hatch=self.hatch,
            color=self.color_l,
            edgecolor="black",
            linewidth=1.0,
        )
        ax_l.set_xlabel("")
        ax_l.set_ylabel("Average BIB yield per event")
        ax_l.semilogy()
        ymin, ymax = ax_l.get_ylim()
        ax_l.set_ylim(0.08, ymax * 1.35)

        # force y-axis ticks at 1e-1 to 1e7, with minor ticks at 2, 3, 4, 5, 6, 7, 8, 9 times each power of ten
        ax_l.set_yticks([10 ** i for i in range(-1, 8)], minor=False)
        ax_l.set_yticks([j * 10 ** i for i in range(-1, 8) for j in range(2, 10)], minor=True)

        # plot signal efficiency
        if do_signal:
            ax_r = ax_l.twinx()
            ax_r.plot(x_vals, y_vals_r, marker="o", color=self.color_r)
            ax_r.set_ylabel(r"Muon efficiency for $\geq 1$ reconstructed object", color=self.color_r, labelpad=self.pad_r)
            ax_r.set_ylim(0.5, 1.02)
            ax_r.tick_params(axis="y", colors=self.color_r)

            # add bonus text
            # textfill = {"bbox": dict(boxstyle="square,pad=0.2", ec="black", fc=self.color_l)}
            ax_l.text(0.02, 1.015, "Background", color="black", transform=ax_l.transAxes)
            ax_r.text(0.87, 1.010, "Signal", color=self.color_r, transform=ax_r.transAxes)

            # add line to right y-axis for efficiency = 1
            ax_r.axhline(y=1.0, xmin=0.5, color=self.color_r, linestyle="--", linewidth=1.0)
        else:
            ax_l.grid(True)
            ax_l.set_title("Inner and outer trackers, barrel")

        # save
        pdf.savefig(fig)
        plt.close(fig)



def options():
    _num = "1"
    _default = {
        "--background-mcps": f"v01_background100_digi_10um/mcps_{_num}.pkl",
        "--background-hits": f"v01_background100_digi_10um/simhits_{_num}.pkl",
        "--background-mds": f"v01_background100_digi_10um/mds_{_num}.pkl",
        "--background-t2s": f"v01_background100_digi_10um/t2s_{_num}.pkl",
        "--background-t4s": f"v01_background100_digi_10um/t4s_{_num}.pkl",
        "--background-t8s": f"v01_background100_digi_10um/t8s_{_num}.pkl",
        "--signal-mcps": f"v01_signal_digi_10um/mcps.pkl",
        "--signal-hits": f"v01_signal_digi_10um/hits.pkl",
        "--signal-mds": f"v01_signal_digi_10um/mds.pkl",
        "--signal-t2s": f"v01_signal_digi_10um/t2s.pkl",
        "--signal-t4s": f"v01_signal_digi_10um/t4s.pkl",
        "--signal-t8s": f"v01_signal_digi_10um/t8s.pkl",
    }
    parser = argparse.ArgumentParser(description="Demo for 2D scatter plot")
    parser.add_argument("--background-mcps", default=_default["--background-mcps"], help="Comma-separated list of input mcps file paths")
    parser.add_argument("--background-hits", default=_default["--background-hits"], help="Comma-separated list of input hits file paths")
    parser.add_argument("--background-mds", default=_default["--background-mds"], help="Comma-separated list of input mds file paths")
    parser.add_argument("--background-t2s", default=_default["--background-t2s"], help="Comma-separated list of input t2s file paths")
    parser.add_argument("--background-t4s", default=_default["--background-t4s"], help="Comma-separated list of input t4s file paths")
    parser.add_argument("--background-t8s", default=_default["--background-t8s"], help="Comma-separated list of input t8s file paths")
    parser.add_argument("--signal-mcps", default=_default["--signal-mcps"], help="Comma-separated list of input signal mcps file paths")
    parser.add_argument("--signal-hits", default=_default["--signal-hits"], help="Comma-separated list of input signal file paths")
    parser.add_argument("--signal-mds", default=_default["--signal-mds"], help="Comma-separated list of input signal mds file paths")
    parser.add_argument("--signal-t2s", default=_default["--signal-t2s"], help="Comma-separated list of input signal t2s file paths")
    parser.add_argument("--signal-t4s", default=_default["--signal-t4s"], help="Comma-separated list of input signal t4s file paths")
    parser.add_argument("--signal-t8s", default=_default["--signal-t8s"], help="Comma-separated list of input signal t8s file paths")
    parser.add_argument("--pdf", default="yields_and_efficiency.pdf", help="Output pdf file path")
    return parser.parse_args()


rcParams.update({
    "font.size": 16,
    "figure.figsize": (8, 8),
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    # "axes.grid": True,
    # "axes.grid.which": "both",
    "axes.axisbelow": True,
    "grid.linewidth": 0.5,
    "grid.alpha": 0.1,
    "grid.color": "gray",
    "figure.subplot.left": 0.13,
    "figure.subplot.bottom": 0.07,
    "figure.subplot.right": 0.90,
    "figure.subplot.top": 0.95,
})

if __name__ == "__main__":
    main()
