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
        for key, objs in data_paths.items():
            self.paths[key] = self.get_input_filenames(objs)
        self.pdf_path = pdf_path
        self.load_data()


    def get_input_filenames(self, input_str):
        fnames = []
        for pattern in input_str.split(","):
            fnames.extend(glob(pattern))
        return fnames


    def load_data(self):
        save_keys = ["file", "i_event", "i_mcp"]
        self.df = {}
        for key, file_paths in self.paths.items():
            logger.info(f"Loading {key} from {file_paths} ...")
            # tmp
            if key == "background_t8s":
                logger.info(f"EXCEPTIONALLY SKIPPING {key} ...")
                self.df[key] = pd.DataFrame(columns=save_keys)
                continue
            # /tmp
            self.df[key] = pd.concat([pd.read_pickle(fi)[save_keys] for fi in file_paths], ignore_index=True)
        self.announce_sizes()


    def announce_sizes(self):
        for key in self.df.keys():
            logger.info(f"{key}: {len(self.df[key])} entries")


    def plot(self):
        with PdfPages(self.pdf_path) as pdf:
            self.plot_yields_and_efficiency(pdf)


    def plot_yields_and_efficiency(self, pdf):
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_title("BIB Yields and Muon Efficiency")
        ax.set_xlabel("Object abstraction")
        ax.set_ylabel("Average per event")
        ax.plot([len(self.df["background_hits"]),
                 len(self.df["background_mds"]),
                 len(self.df["background_t2s"]),
                 len(self.df["background_t4s"]),
                 len(self.df["background_t8s"])],
                label="Background",
                marker="o",
                alpha=0.5,
        )
        ax.semilogy()
        ax.set_ylim(0.8, None)
        pdf.savefig(fig)
        plt.close(fig)



def options():
    _num = 0
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
