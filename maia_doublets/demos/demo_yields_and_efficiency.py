"""

"""
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
ZERO_EVENTS_AT_95_CL = 3.0

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

    data = DataGrabber(paths, args.pdf, args.precision_timing)
    data.plot()


class DataGrabber:

    def __init__(self,
                 data_paths: dict,
                 pdf_path: str,
                 do_precision_timing: bool,
                 ):
        self.paths = {}
        self.nfiles_per_key = {}
        for key, objs in data_paths.items():
            self.paths[key] = self.get_input_filenames(objs)
            self.nfiles_per_key[key] = len(self.paths[key])
        self.pdf_path = pdf_path
        self.do_precision_timing = do_precision_timing
        self.unique_cols = ["file", "i_event", "i_mcp"]
        self.columns_background = ["file"]
        self.columns_signal = self.unique_cols + ["mcp_pdg",
                                                  "mcp_q",
                                                  "mcp_pt",
                                                  "mcp_eta",
                                                  "mcp_phi",
                                                  "mcp_vertex_r",
                                                  "mcp_vertex_z"]
        self.columns_mcp = self.columns_signal + ["mcp_detectable_ITB",
                                                  "mcp_detectable_OTB"]
        self.load_data()
        self.load_precision_timing_data()
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
            columns = self.columns_signal if key.startswith("signal") else self.columns_background
            if key == "signal_mcps":
                columns = self.columns_mcp
            logger.info(f"Loading {key} from {file_paths} ...")
            for fi in file_paths:
                if len(pd.read_pickle(fi)) > 0:
                    break
            else:
                logger.warning(f"No non-empty files found for {key} in {file_paths}. Skipping.")
                self.df[key] = pd.DataFrame(columns=columns)
                continue
            self.df[key] = pd.concat([pd.read_pickle(fi)[columns] for fi in file_paths], ignore_index=True)
        self.announce_sizes()


    def announce_sizes(self):
        for key in self.df.keys():
            logger.info(f"{key}: {len(self.df[key])} entries")
        for key, nfiles in self.nfiles_per_key.items():
            logger.info(f"{key}: {nfiles} files")


    def load_precision_timing_data(self):
        """
        See: https://github.com/MuonColliderSoft/MAIAConfig/blob/main/MAIAConfig/TrackerDigi/
        """
        if not self.do_precision_timing:
            return
        key = "background_hits"
        time_key = "simhit_t_corrected"
        TimeWindowMax = 0.3
        TimeWindowMin = -0.18
        logger.info(f"Loading precision timing data for {key} ...")
        time_df = pd.concat([pd.read_pickle(fi)[time_key] for fi in self.paths[key]], ignore_index=True)
        n_pass = ((time_df > TimeWindowMin) & (time_df < TimeWindowMax)).sum()
        n_files = self.nfiles_per_key[key]
        logger.info(f"Precision timing: {len(time_df)} hits in {n_files} files")
        logger.info(f"Precision timing: {n_pass} hits within time window")
        logger.info(f"Precision timing: {n_pass / n_files} hits per file within time window")
        self.n_precision_timing = n_pass / n_files


    def get_denominator_mask(self):
        mask = (
            (np.abs(self.df["signal_mcps"]["mcp_pdg"]) == MUON) &
            (self.df["signal_mcps"]["mcp_q"] != 0) &
            (self.df["signal_mcps"]["mcp_pt"] > ONE_POINT_FIVE_GEV) &
            (self.df["signal_mcps"]["mcp_vertex_r"] < ZERO_POINT_ZERO_ONE_MM) &
            (np.abs(self.df["signal_mcps"]["mcp_vertex_z"]) < ZERO_POINT_ZERO_ONE_MM) &
            (np.abs(self.df["signal_mcps"]["mcp_eta"]) < BARREL_TRACKER_MAX_ETA)
        )
        return mask


    def get_algorithmic_mask(self):
        mask = (
            self.get_denominator_mask() &
            (self.df["signal_mcps"]["mcp_detectable_ITB"] == True) &
            (self.df["signal_mcps"]["mcp_detectable_OTB"] == True)
        )
        return mask


    def calculate_efficiency(self):
        overall_mask = self.get_denominator_mask()
        algorithmic_mask = self.get_algorithmic_mask()
        if overall_mask.sum() == 0:
            msg = "No signal mcps found that satisfy denominator criteria"
            logger.error(msg)
            raise ValueError(msg)

        overall_denom = self.df["signal_mcps"][overall_mask]
        algorithmic_denom = self.df["signal_mcps"][algorithmic_mask]
        logger.info(f"Denominator mask: {overall_mask.sum()} entries")
        logger.info(f"Algorithmic mask: {algorithmic_mask.sum()} entries")

        self.overall_efficiency = {}
        self.algorithmic_efficiency = {}
        for key in self.df.keys():
            if key.startswith("background"):
                continue

            # for each row in df, check if it exists in mcps
            df = self.df[key].drop_duplicates(subset=self.unique_cols)
            overall_merged = pd.merge(df, overall_denom, on=self.unique_cols, how="inner")
            algorithmic_merged = pd.merge(df, algorithmic_denom, on=self.unique_cols, how="inner")

            # calculate efficiency
            overall_eff = len(overall_merged) / len(overall_denom)
            algorithmic_eff = len(algorithmic_merged) / len(algorithmic_denom)
            self.overall_efficiency[key] = overall_eff
            self.algorithmic_efficiency[key] = algorithmic_eff
            logger.info(f"Overall efficiency for {key}: {len(overall_merged)} / {len(overall_denom)} = {overall_eff:.6f}")
            logger.info(f"Algorithmic efficiency for {key}: {len(algorithmic_merged)} / {len(algorithmic_denom)} = {algorithmic_eff:.6f}")


    def plot(self):
        with PdfPages(self.pdf_path) as pdf:
            self.plot_yields_and_efficiency(pdf, do_signal=True)
            self.plot_yields_and_efficiency(pdf, do_signal=False)
            if self.do_precision_timing:
                self.plot_yields_and_efficiency(pdf, do_signal=True, do_precision_timing=True)
            self.plot_efficiency_vs_kinematics(pdf)


    def plot_yields_and_efficiency(self,
                                   pdf,
                                   do_signal,
                                   do_precision_timing=False,
                                   ):
        y_vals_l = [
            len(self.df["background_hits"]) / self.nfiles_per_key["background_hits"],
            len(self.df["background_mds"]) / self.nfiles_per_key["background_mds"],
            len(self.df["background_t2s"]) / self.nfiles_per_key["background_t2s"],
            len(self.df["background_t4s"]) / self.nfiles_per_key["background_t4s"],
            len(self.df["background_t8s"]) / self.nfiles_per_key["background_t8s"],
        ]
        if do_signal:
            y_vals_r_overall = [
                self.overall_efficiency["signal_hits"],
                self.overall_efficiency["signal_mds"],
                self.overall_efficiency["signal_t2s"],
                self.overall_efficiency["signal_t4s"],
                self.overall_efficiency["signal_t8s"],
            ]
            y_vals_r_algorithmic = [
                self.algorithmic_efficiency["signal_hits"],
                self.algorithmic_efficiency["signal_mds"],
                self.algorithmic_efficiency["signal_t2s"],
                self.algorithmic_efficiency["signal_t4s"],
                self.algorithmic_efficiency["signal_t8s"],
            ]
        else:
            y_vals_r_overall = [0] * len(y_vals_l)
            y_vals_r_algorithmic = [0] * len(y_vals_l)
        x_vals = ["Hits", "MDs", "T2s", "T4s", "T8s"]
        bins = np.arange(len(x_vals) + 1) - 0.5

        # check
        if (len(y_vals_l) != len(y_vals_r_overall) or
            len(y_vals_l) != len(y_vals_r_algorithmic) or
            len(y_vals_l) != len(x_vals)
            ):
            msg = f"Length mismatch: y_vals_l={len(y_vals_l)}, y_vals_r_overall={len(y_vals_r_overall)}, y_vals_r_algorithmic={len(y_vals_r_algorithmic)}, x_vals={len(x_vals)}"
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
        ymin, ymax = ax_l.get_ylim()

        # make a horizontal dotted line in case no tracks are found in the last bin
        if y_vals_l[-1] == 0:
            upper_limit = ZERO_EVENTS_AT_95_CL / self.nfiles_per_key["background_t8s"]
            logger.info(f"Upper limit for last bin: {upper_limit} (95% CL tracks per event)")
            if upper_limit > ymin:
                center = len(x_vals) - 1
                # ax_l.hlines(y=upper_limit,
                #             xmin=center - 0.15,
                #             xmax=center + 0.15,
                #             color="black")
                # ax_l.annotate("",
                #               xy=(center, upper_limit),
                #               xytext=(center, upper_limit / 2.5),
                #               arrowprops=dict(arrowstyle="<-"),
                #               )

        # force y-axis ticks at 1e-1 to 1e7, with minor ticks at 2, 3, 4, 5, 6, 7, 8, 9 times each power of ten
        ax_l.set_yticks([10 ** i for i in range(-1, 8)], minor=False)
        ax_l.set_yticks([j * 10 ** i for i in range(-1, 8) for j in range(2, 10)], minor=True)

        # add line where number of hits with precision timing would be
        if do_precision_timing:
            ax_l.axhline(y=self.n_precision_timing,
                        color="red",
                        xmin=0.050, # suffering
                        xmax=0.222, # suffering
                        linestyle="--",
                        linewidth=2.0,
                        )

        # plot signal efficiency
        if do_signal:
            ax_r = ax_l.twinx()

            # plot overall efficiency as filled circle
            ax_r.plot(x_vals, y_vals_r_overall, marker="o", color=self.color_r)

            # plot algorithmic efficiency as unfilled circle
            ax_r.plot(x_vals, y_vals_r_algorithmic, marker="o", color=self.color_r, fillstyle="none")
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


    def plot_efficiency_vs_kinematics(self, pdf: PdfPages):

        # kinematics
        kin_cols = [
            # "mcp_q",
            "mcp_pt",
            "mcp_eta",
            "mcp_phi",
        ]
        self.bins = {
            "mcp_q": np.array([-2, 0, 2]),
            "mcp_pt": np.unique(np.concatenate([np.linspace(0, 2, 21), # 0.1
                                                np.linspace(2, 3, 6), # 0.2
                                                np.linspace(3, 5, 7), # 0.333
                                                np.linspace(5, 10, 11), # 0.5
                                                ])),
            "mcp_eta": np.linspace(-0.7, 0.7, 141),
            "mcp_phi": np.linspace(-3.2, 3.2, 161),
        }
        self.xlabel = {
            "mcp_q": "Inclusive",
            "mcp_pt": r"Muon $p_T$ [GeV]",
            "mcp_eta": r"Muon $\eta$",
            "mcp_phi": r"Muon $\phi$ [rad]",
        }

        # denominators
        overall_mask = self.get_denominator_mask()
        algorithmic_mask = self.get_algorithmic_mask()
        overall_denom = self.df["signal_mcps"][overall_mask][self.unique_cols + kin_cols]
        algorithmic_denom = self.df["signal_mcps"][algorithmic_mask][self.unique_cols + kin_cols]
        if overall_denom.duplicated().any():
            raise ValueError("Denominator has duplicated rows!")

        for obj in ["signal_t4s", "signal_t8s"]:

            objects = self.df[obj][self.unique_cols].drop_duplicates()
            overall_numer = pd.merge(overall_denom, objects, on=self.unique_cols, how="inner")
            algorithmic_numer = pd.merge(algorithmic_denom, objects, on=self.unique_cols, how="inner")
            print(f"All numerator for {obj}: {len(overall_numer)} entries")
            print(f"Alg numerator for {obj}: {len(algorithmic_numer)} entries")
            print(f"All denomator for {obj}: {len(overall_denom)} entries")
            print(f"Alg denomator for {obj}: {len(algorithmic_denom)} entries")

            # if not overall_numer.equals(algorithmic_numer):
            #     raise ValueError(f"Overall and algorithmic numerators are different for {obj}!")

            for kin in kin_cols:
                n_denom_overall, edges = np.histogram(overall_denom[kin], bins=self.bins[kin])
                n_numer_overall, edges = np.histogram(overall_numer[kin], bins=self.bins[kin])
                eff_overall = np.divide(n_numer_overall, n_denom_overall, out=np.zeros_like(n_numer_overall, dtype=float), where=n_denom_overall!=0)
                n_denom_algorithmic, edges = np.histogram(algorithmic_denom[kin], bins=self.bins[kin])
                n_numer_algorithmic, edges = np.histogram(algorithmic_numer[kin], bins=self.bins[kin])
                eff_algorithmic = np.divide(n_numer_algorithmic, n_denom_algorithmic, out=np.zeros_like(n_numer_algorithmic, dtype=float), where=n_denom_algorithmic!=0)
                centers = 0.5 * (edges[1:] + edges[:-1])
                fig, ax = plt.subplots()
                color_alg = "#83BAE2"
                pargs = dict(marker="o", markersize=1, linewidth=4.0, linestyle="-")
                targs = dict(fontsize=20, va="top", ha="left", transform=ax.transAxes)
                ax.plot(centers, eff_algorithmic, color=color_alg, **pargs)
                ax.plot(centers, eff_overall, color=self.color_r, **pargs)
                ax.text(0.50, 0.30, "Algorithm efficiency", color=color_alg, **targs)
                ax.text(0.50, 0.24, "Overall efficiency", color=self.color_r, **targs)
                ax.set_xlabel(self.xlabel[kin])
                ax.set_ylabel("Efficiency")
                ax.set_title(f"{obj}")
                ax.set_ylim(0.65, 1.04)
                ax.grid(True)
                pdf.savefig()
                plt.close()



def options():
    _num = "1"
    _default = {
        "--background-mcps": f"v01_background100_digi_10um/mcps_{_num}.pkl",
        "--background-hits": f"v01_background100_digi_10um/simhits_{_num}.pkl",
        "--background-mds": f"v01_background100_digi_10um/mds_{_num}.pkl",
        "--background-t2s": f"v01_background100_digi_10um/t2s_{_num}.pkl",
        "--background-t4s": f"v01_background100_digi_10um/t4s_{_num}.pkl",
        "--background-t8s": f"v01_background100_digi_10um/t8s_{_num}.pkl",
        # "--signal-mcps": f"v01_signal_digi_10um/mcps.pkl",
        # "--signal-hits": f"v01_signal_digi_10um/hits.pkl",
        # "--signal-mds": f"v01_signal_digi_10um/mds.pkl",
        # "--signal-t2s": f"v01_signal_digi_10um/t2s.pkl",
        # "--signal-t4s": f"v01_signal_digi_10um/t4s.pkl",
        # "--signal-t8s": f"v01_signal_digi_10um/t8s.pkl",
        "--signal-mcps": f"v01_muonGun_pT_0_10_digi_10um/mcps.pkl",
        "--signal-hits": f"v01_muonGun_pT_0_10_digi_10um/hits.pkl",
        "--signal-mds": f"v01_muonGun_pT_0_10_digi_10um/mds.pkl",
        "--signal-t2s": f"v01_muonGun_pT_0_10_digi_10um/t2s.pkl",
        "--signal-t4s": f"v01_muonGun_pT_0_10_digi_10um/t4s.pkl",
        "--signal-t8s": f"v01_muonGun_pT_0_10_digi_10um/t8s.pkl",
    }
    parser = argparse.ArgumentParser(description="Demo for 2D scatter plot")
    parser.add_argument("--precision-timing", action="store_true", help="Mention precision timing for the plots")
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
