import inspect
import textwrap
import time
import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)

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

from constants import MUON
from constants import BARREL_TRACKER_MAX_ETA
from constants import BARREL_TRACKER_MAX_RADIUS
from constants import ONE_POINT_FIVE_GEV, ONE_MM, ZERO_POINT_ZERO_ONE_MM
from constants import NICKNAMES, OUTER_TRACKER_BARREL
from constants import MD_DZ_CUT, MD_DR_CUT
from constants import REQ_PASSTHROUGH, REQ_RZ, REQ_XY, REQ_RZ_XY
from constants import DOUBLET_REQS, NO_MCP
from constants import LS_REQS, LS_REQ_DR_POS, LS_REQ_DZ_POS, LS_REQ_XY_CHI2, LS_REQ_RZ_ANG, LS_REQ_ALL
from constants import LS_DZ_CUT, LS_DR_CUT, LS_DTHETA_RZ_CUT, LS_DTHETA_XY_CUT, LS_CHI2_XY_CUT
from constants import MIN_COSTHETA, MIN_SIMHIT_PT_FRACTION, MAX_TIME
from constants import N_LS_PHI_SLICES
from constants import T4_DR_CUT, T4_DZ_CUT, T4_DTHETA_RZ_CUT, T4_CHI2_XY_CUT


class Plotter:

    def __init__(
        self,
        geometry_version: str,
        sim: bool,
        smear: str,
        signal: bool,
        mcps: pd.DataFrame,
        simhits: pd.DataFrame,
        doublets: pd.DataFrame,
        linesegments: pd.DataFrame,
        t4s: pd.DataFrame,
        pdf: str,
    ):
        self.signal = signal
        self.mcps = mcps
        self.simhits = simhits
        self.doublets = doublets
        self.linesegments = linesegments
        self.t4s = t4s
        self.pdf = pdf

        # shorthands for cuts
        key = (geometry_version, "sim") if sim else (geometry_version, "digi", smear)
        self.MD_DZ_CUT = MD_DZ_CUT[key]
        self.MD_DR_CUT = MD_DR_CUT[key]
        self.LS_DZ_CUT = LS_DZ_CUT[key]
        self.LS_DR_CUT = LS_DR_CUT[key]
        self.LS_DTHETA_RZ_CUT = LS_DTHETA_RZ_CUT[key]
        self.LS_DTHETA_XY_CUT = LS_DTHETA_XY_CUT[key]
        self.LS_CHI2_XY_CUT = LS_CHI2_XY_CUT[key]
        self.T4_DR_CUT = T4_DR_CUT[key]
        self.T4_DZ_CUT = T4_DZ_CUT[key]
        self.T4_DTHETA_RZ_CUT = T4_DTHETA_RZ_CUT[key]
        self.T4_CHI2_XY_CUT = T4_CHI2_XY_CUT[key]


    def plot(self):
        logger.info(f"Writing plots to {self.pdf} ...")
        with PdfPages(self.pdf) as pdf:
            self.plot_numbers_for_comparison(pdf)
            self.write_date(pdf)
            self.write_quality_cuts(pdf)
            self.plot_time(pdf)
            # self.plot_layer_occupancy_1d(pdf)
            # self.plot_layer_occupancy_2d(pdf)
            self.plot_radius_vs_layer(pdf)
            # self.plot_doublet_occupancy(pdf)
            # self.plot_doublet_features(pdf)
            # self.plot_linesegment_features(pdf)
            # self.plot_t4_features(pdf)
            if self.signal:
                self.write_denominator_info(pdf)
                self.plot_detectable_efficiency_vs_kinematics(pdf)
                self.plot_doublet_efficiency_vs_kinematics_2(pdf)
                # self.plot_doublet_efficiency_vs_kinematics(pdf)
                # self.write_doublet_denominator_info(pdf)
                # self.plot_doublet_quality_efficiency(pdf)
                # self.plot_segment_efficiency_vs_kinematics(pdf)
                # self.plot_segment_quality_efficiency(pdf)
                # self.plot_t4_efficiency_vs_kinematics(pdf)
                # self.plot_t4_quality_efficiency(pdf)


    def plot_numbers_for_comparison(self, pdf: PdfPages):
        """
        Cutflow comparison!
        """
        if self.signal:
            self.plot_numbers_for_comparison_signal(pdf)
        else:
            self.plot_numbers_for_comparison_background(pdf)


    def plot_numbers_for_comparison_signal(self, pdf: PdfPages):

        # part 1: simhits
        mask = np.ones(len(self.simhits), dtype=bool)
        for [req, label] in [
            [self.simhits["simhit_layer"].isin([0, 1]), "All simhits in layers 0 and 1"],
            [np.abs(self.simhits["mcp_pdg"]) == MUON, "abs(pdg) == muon"],
            [self.simhits["mcp_q"] != 0, "q is not 0"],
            [self.simhits["mcp_pt"] > ONE_POINT_FIVE_GEV, "pT > 1.5 GeV"],
            [np.abs(self.simhits["mcp_eta"]) < BARREL_TRACKER_MAX_ETA, f"abs(eta) < {BARREL_TRACKER_MAX_ETA}"],
            [self.simhits["mcp_vertex_r"] < ZERO_POINT_ZERO_ONE_MM, "vertex r < 0.01 mm"],
            [np.abs(self.simhits["mcp_vertex_z"]) < ZERO_POINT_ZERO_ONE_MM, "abs(vertex z) < 0.01 mm"],
            [self.simhits["simhit_t_corrected"] < MAX_TIME, f"corrected t < {MAX_TIME} ns"],
            [self.simhits["simhit_costheta"] > MIN_COSTHETA, f"costheta > {MIN_COSTHETA}"],
            [self.simhits["simhit_p"] / self.simhits["mcp_p"] > MIN_SIMHIT_PT_FRACTION, f"simhit p / mcp p > {MIN_SIMHIT_PT_FRACTION}"],
            # [self.simhits["simhit_sensor"] == 20, "z-sensor 20"],
            # [self.simhits["simhit_module"] == 0, "phi-module 0"],
        ]:
            mask &= req
            logger.info(f"* {label:<30} :: {mask.sum():>10}")

        # part 2a: doublets by hand
        doublet_cols = [
            "file",
            "i_event", # the event
            "simhit_system", # the system (IT, OT)
            "simhit_layer_div_2", # the double layer
            "simhit_module", # the phi-module
            "simhit_sensor", # the z-sensor
        ]
        lower_mask = mask & (self.simhits["simhit_layer_mod_2"] == 0)
        upper_mask = mask & (self.simhits["simhit_layer_mod_2"] == 1)
        logger.info(f"* {'Lower hit':<30} :: {lower_mask.sum():>10}")
        logger.info(f"* {'Upper hit':<30} :: {upper_mask.sum():>10}")

        # logger.info("Getting lower and upper hits ...")
        lower = self.simhits[lower_mask].rename(columns={"i_mcp": "i_mcp_lower"})[doublet_cols + ["i_mcp_lower"]]
        upper = self.simhits[upper_mask].rename(columns={"i_mcp": "i_mcp_upper"})[doublet_cols + ["i_mcp_upper"]]
        doublets = lower.merge(upper, on=doublet_cols, how="inner")
        logger.info(f"* {'Doublets':<30} :: {len(doublets):>10}")

        same_mcp = doublets["i_mcp_lower"] == doublets["i_mcp_upper"]
        logger.info(f"* {'Doublets from same MCP':<30} :: {same_mcp.sum():>10}")

        # part 2b: doublets with dr and dz cuts
        doublelayer = self.doublets["doublet_doublelayer"]
        dl_0 = doublelayer == 0
        dl_1 = doublelayer == 1
        baseline_cuts = (
            (self.doublets["i_mcp"] >= 0) &
            (np.abs(self.doublets["mcp_pdg"]) == MUON) &
            (self.doublets["mcp_q"] != 0) &
            (self.doublets["mcp_pt"] > ONE_POINT_FIVE_GEV) &
            (np.abs(self.doublets["mcp_eta"]) < BARREL_TRACKER_MAX_ETA) &
            (self.doublets["mcp_vertex_r"] < ZERO_POINT_ZERO_ONE_MM) &
            (np.abs(self.doublets["mcp_vertex_z"]) < ZERO_POINT_ZERO_ONE_MM) &
            (self.doublets["doublet_first_exit"])
        )
        quality_cuts = baseline_cuts & self.doublets["doublet_ok"]
        doublets_0 = self.doublets[baseline_cuts & dl_0]
        doublets_1 = self.doublets[baseline_cuts & dl_1]
        logger.info(f"* {'Doublets, baseline, L01':<30} :: {len(doublets_0):>10}")
        logger.info(f"* {'Doublets, baseline, L23':<30} :: {len(doublets_1):>10}")
        doublets_0 = self.doublets[quality_cuts & dl_0]
        doublets_1 = self.doublets[quality_cuts & dl_1]
        logger.info(f"* {'Doublets, drdz cuts, L01':<30} :: {len(doublets_0):>10}")
        logger.info(f"* {'Doublets, drdz cuts, L23':<30} :: {len(doublets_1):>10}")

        # part 3: line segments
        keys = [
            "file",
            "i_event",
            "i_mcp",
        ]
        segments = doublets_0.merge(
            doublets_1,
            on=keys,
            how="inner",
            validate="many_to_many",
            suffixes=("_lower", "_upper"),
        )
        logger.info(f"* {'Line segments from doublets':<30} :: {len(segments):>10}")


    def plot_numbers_for_comparison_background(self, pdf: PdfPages):
        layers = [0, 1]
        # layers = [2, 3]
        the_doublelayer = layers[0] // 2

        # part 1: simhits
        mask = np.ones(len(self.simhits), dtype=bool)
        for [req, label] in [
            [self.simhits["simhit_layer"].isin(layers), f"All simhits in layers {layers}"],
            # [self.simhits["simhit_sensor"] == 20, "z-sensor 20"],
            # [self.simhits["simhit_module"] == 0, "phi-module 0"],
        ]:
            mask &= req
            logger.info(f"* {label:<30} :: {mask.sum():>10}")

        # part 2: doublets
        doublet_cols = [
            "file",
            "i_event", # the event
            "simhit_system", # the system (IT, OT)
            "simhit_layer_div_2", # the double layer
            "simhit_module", # the phi-module
            "simhit_sensor", # the z-sensor
        ]
        lower_mask = mask & (self.simhits["simhit_layer_mod_2"] == 0)
        upper_mask = mask & (self.simhits["simhit_layer_mod_2"] == 1)
        logger.info(f"* {'Lower hit':<30} :: {lower_mask.sum():>10}")
        logger.info(f"* {'Upper hit':<30} :: {upper_mask.sum():>10}")

        # number of doublets by hand
        # lower = self.simhits[lower_mask][doublet_cols]
        # upper = self.simhits[upper_mask][doublet_cols]
        # doublets = lower.merge(upper, on=doublet_cols, how="inner")
        # logger.info(f"* {'Doublets (by hand)':<30} :: {len(doublets):>10}")

        # number of doublets
        mask = np.ones(len(self.doublets), dtype=bool)
        for [req, label] in [
            [self.doublets["doublet_system"] == OUTER_TRACKER_BARREL, "Doublets in OTB"],
            [self.doublets["doublet_doublelayer"] == the_doublelayer, f"Doublets in layers {layers}"],
            # [self.doublets["doublet_sensor"] == 20, "z-sensor 20"],
            # [self.doublets["doublet_module"] == 0, "phi-module 0"],
            [np.abs(self.doublets["doublet_dz"]) < self.MD_DZ_CUT[the_doublelayer], f"Doublets with |dz| < {self.MD_DZ_CUT[the_doublelayer]}mm"],
            [np.abs(self.doublets["doublet_dr"]) < self.MD_DR_CUT[the_doublelayer], f"Doublets with |dr| < {self.MD_DR_CUT[the_doublelayer]}mm"],
        ]:
            mask &= req
            logger.info(f"* {label:<30} :: {mask.sum():>10}")

        # part 3: line segments
        if self.linesegments is None:
            logger.info("No line segments, skipping ...")
        else:
            mask = np.ones(len(self.linesegments), dtype=bool)
            for [req, label] in [
                [self.linesegments["ls_system"] == OUTER_TRACKER_BARREL, "LS in OTB"],
                [self.linesegments["ls_doublelayer"] == the_doublelayer, f"LS starting on layer {the_doublelayer}"],
                [np.abs(self.linesegments["ls_dz"]) < self.LS_DZ_CUT[the_doublelayer], f"LS with |dz| < {self.LS_DZ_CUT[the_doublelayer]}mm"],
                [np.abs(self.linesegments["ls_dr"]) < self.LS_DR_CUT[the_doublelayer], f"LS with |dr| < {self.LS_DR_CUT[the_doublelayer]}mm"],
                [np.abs(self.linesegments["ls_dtheta_rz"]) < self.LS_DTHETA_RZ_CUT[the_doublelayer], f"LS with |dtheta_rz| < {self.LS_DTHETA_RZ_CUT[the_doublelayer]}"],
                [np.abs(self.linesegments["ls_dtheta_xy"]) < self.LS_DTHETA_XY_CUT[the_doublelayer], f"LS with |dtheta_xy| < {self.LS_DTHETA_XY_CUT[the_doublelayer]}"],
                [np.abs(self.linesegments["ls_chi2_012"]) < self.LS_CHI2_XY_CUT[the_doublelayer], f"LS with |chi2_xy| < {self.LS_CHI2_XY_CUT[the_doublelayer]}"],
            ]:
                mask &= req
                logger.info(f"* {label:<30} :: {mask.sum():>10}")


    def write_date(self, pdf: PdfPages):
        logger.info(f"Writing the date and time")
        now = time.strftime("%Y_%m_%d_%Hh%Mm%Ss")
        text = f"Plots made at {now}"
        fig, ax = plt.subplots(figsize=(8, 8))
        args = {"ha":"left", "va":"top", "fontfamily":"monospace"}
        ax.text(0.0, 0.7, text, **args, fontsize=16)
        ax.axis("off")
        pdf.savefig()
        plt.close()


    def write_quality_cuts(self, pdf: PdfPages):
        logger.info(f"Writing quality cuts")
        text = textwrap.dedent(f"""\
            Quality cuts for MDs:
            * |dr| < {self.MD_DR_CUT[0]} mm (DL0)
            * |dr| < {self.MD_DR_CUT[1]} mm (DL1)
            * |dr| < {self.MD_DR_CUT[2]} mm (DL2)
            * |dr| < {self.MD_DR_CUT[3]} mm (DL3)
            * |dz| < {self.MD_DZ_CUT[0]} mm (DL0)
            * |dz| < {self.MD_DZ_CUT[1]} mm (DL1)
            * |dz| < {self.MD_DZ_CUT[2]} mm (DL2)
            * |dz| < {self.MD_DZ_CUT[3]} mm (DL3)

            Quality cuts for T2s:
            * |dr| < {self.LS_DR_CUT[0]} mm (DL0)
            * |dr| < {self.LS_DR_CUT[1]} mm (DL1)
            * |dr| < {self.LS_DR_CUT[2]} mm (DL2)
            * |dz| < {self.LS_DZ_CUT[0]} mm (DL0)
            * |dz| < {self.LS_DZ_CUT[1]} mm (DL1)
            * |dz| < {self.LS_DZ_CUT[2]} mm (DL2)
            * |dtheta_rz| < {self.LS_DTHETA_RZ_CUT[0]} (DL0)
            * |dtheta_rz| < {self.LS_DTHETA_RZ_CUT[1]} (DL1)
            * |dtheta_rz| < {self.LS_DTHETA_RZ_CUT[2]} (DL2)
            * |chi2_xy| < {self.LS_CHI2_XY_CUT[0]} (DL0)
            * |chi2_xy| < {self.LS_CHI2_XY_CUT[1]} (DL1)
            * |chi2_xy| < {self.LS_CHI2_XY_CUT[2]} (DL2)

            Quality cuts for T4s:
            * |dr| < {self.T4_DR_CUT[0]} mm (DL0)
            * |dz| < {self.T4_DZ_CUT[0]} mm (DL0)
            * |dtheta_rz| < {self.T4_DTHETA_RZ_CUT[0]} (DL0)
            * |chi2_xy| < {self.T4_CHI2_XY_CUT[0]} (DL0)
        """)
        fig, ax = plt.subplots(figsize=(8, 8))
        args = {"ha":"left", "va":"top", "fontfamily":"monospace"}
        ax.text(0.0, 0.9, text, **args, fontsize=12)
        ax.axis("off")
        pdf.savefig()
        plt.close()


    def plot_time(self, pdf: PdfPages):
        logger.info(f"Plotting time")
        xlabel = "Sim. hit time [ns]" + r" minus $R/c$"
        for (system, simhits) in self.simhits.groupby("simhit_system"):
            bins = np.linspace(-10, 20, 301)
            fig, ax = plt.subplots()
            ax.hist(
                simhits["simhit_t_corrected"],
                bins=bins,
                histtype="stepfilled",
                color="yellow",
                edgecolor="black",
                linewidth=1.0,
                alpha=0.9,
            )
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Sim. hits")
            ax.set_title(f"{NICKNAMES[system]}")
            ax.semilogy()
            ax.set_ylim(0.8, None)
            pdf.savefig()
            plt.close()


    def plot_layer_occupancy_1d(self, pdf: PdfPages):
        logger.info(f"Plotting layer occupancy (1d)")
        for (system, simhits) in self.simhits.groupby("simhit_system"):
            bins = np.arange(simhits["simhit_layer"].min()-0.5,
                             simhits["simhit_layer"].max()+1.0,
                             1)
            fig, ax = plt.subplots()
            ax.hist(
                simhits["simhit_layer"],
                bins=bins,
                histtype="stepfilled",
                color="yellow",
                edgecolor="black",
                linewidth=1.0,
                alpha=0.9,
            )
            # ax.set_ylim(0, 280e3)
            ax.set_xlabel("Layer")
            ax.set_ylabel("Sim. hits")
            ax.set_title(f"{NICKNAMES[system]}")
            pdf.savefig()
            plt.close()


    def plot_layer_occupancy_2d(self, pdf: PdfPages):
        logger.info(f"Plotting layer occupancy (2d)")
        for ((system, layer), simhits) in self.simhits.groupby(["simhit_system",
                                                                "simhit_layer",
                                                                ]):
            logger.info(f"Occupancy of {NICKNAMES[system]} layer {layer}: {len(simhits)} sim hits")
            if len(simhits) == 0:
                continue
            bins = [
                np.arange(-0.5, simhits["simhit_module"].max()+1.5, 1),
                np.arange(-0.5, simhits["simhit_sensor"].max()+1.5, 1),
            ]
            fig, ax = plt.subplots()
            _, _, _, im = ax.hist2d(
                simhits["simhit_module"],
                simhits["simhit_sensor"],
                bins=bins,
                cmap="gist_rainbow",
                norm=colors.LogNorm(vmin=0.9),
            )
            ax.set_xlabel("Phi module")
            ax.set_ylabel("Z sensor")
            ax.set_title(f"{NICKNAMES[system]}, layer {layer}")
            fig.colorbar(im, ax=ax, pad=0.01, label="Number of sim. hits")
            pdf.savefig()
            plt.close()


    def plot_radius_vs_layer(self, pdf: PdfPages):
        logger.info(f"Plotting radius vs layer")
        for (system, simhits) in self.simhits.groupby("simhit_system"):
            bins = [
                np.arange(simhits["simhit_layer"].min() - 0.5,
                          simhits["simhit_layer"].max() + 1.0, 1),
                np.arange(simhits["simhit_r"].min() - 100,
                          simhits["simhit_r"].max() + 100, 1),
                # np.linspace(0, 1600, 1600),
            ]
            fig, ax = plt.subplots()
            _, _, _, im = ax.hist2d(
                simhits["simhit_layer"],
                simhits["simhit_r"],
                bins=bins,
                cmap="gist_rainbow",
                cmin=0.5,
            )
            fig.colorbar(im, ax=ax, label="Number of hits", pad=0.01)
            ax.set_xlabel("Layer")
            ax.set_ylabel("Radius")
            ax.set_title(f"{NICKNAMES[system]}")
            fig.subplots_adjust(right=0.95)
            pdf.savefig()
            plt.close()


    def doublet_requirements(self, doublets: pd.DataFrame, req: str) -> tuple[str, pd.DataFrame]:
        # return description and mask
        doublelayer = doublets["doublet_doublelayer"]
        if len(doublelayer.unique()) != 1:
            raise ValueError(f"Multiple doublelayers found: {doublelayer.unique()}")
        doublelayer = doublelayer.iloc[0]
        if req == REQ_PASSTHROUGH:
            text = "No requirement"
            mask = np.ones(len(doublets), dtype=bool)
        elif req == REQ_XY:
            text = f"|dr| < {self.MD_DR_CUT[doublelayer]}mm"
            mask = np.abs(doublets["doublet_dr"]) < self.MD_DR_CUT[doublelayer]
        elif req == REQ_RZ:
            text = f"|dz| < {self.MD_DZ_CUT[doublelayer]}mm"
            mask = np.abs(doublets["doublet_dz"]) < self.MD_DZ_CUT[doublelayer]
        elif req == REQ_RZ_XY:
            text = f"|dr| < {self.MD_DR_CUT[doublelayer]}mm, |dz| < {self.MD_DZ_CUT[doublelayer]}mm"
            mask = (
                (np.abs(doublets["doublet_dz"]) < self.MD_DZ_CUT[doublelayer]) &
                (np.abs(doublets["doublet_dr"]) < self.MD_DR_CUT[doublelayer])
            )
        else:
            raise ValueError(f"Unknown requirement: {req}")
        return text, mask


    def plot_doublet_occupancy(self, pdf: PdfPages):
        for ((system, doublelayer), group) in self.doublets.groupby(["doublet_system",
                                                                     "doublet_doublelayer",
                                                                    ]):
            zmax = None
            for req in DOUBLET_REQS:

                req_text, req_mask = self.doublet_requirements(group, req)
                doublets = group[req_mask]

                logger.info(f"Occupancy of {NICKNAMES[system]} doublelayer {doublelayer}, {req}: {len(doublets)} doublets")
                if len(doublets) == 0:
                    continue

                layers = [doublelayer * 2, doublelayer * 2 + 1]

                # skip these plots to save time
                # if not self.signal:
                #     continue

                # doublets = self.doublets[mask]
                bins = [
                    np.arange(-0.5, doublets["doublet_module"].max()+1.5, 1),
                    np.arange(-0.5, doublets["doublet_sensor"].max()+1.5, 1),
                ]
                fig, ax = plt.subplots()
                h2d, _, _, im = ax.hist2d(
                    doublets["doublet_module"],
                    doublets["doublet_sensor"],
                    bins=bins,
                    cmap="gist_rainbow",
                    norm=colors.LogNorm(vmin=0.9),
                )
                if zmax is None:
                    zmax = h2d.max()
                im.set_clim(0.9, zmax)
                ax.set_xlabel("Phi module")
                ax.set_ylabel("Z sensor")
                ax.set_title(f"{NICKNAMES[system]}, layers {layers}, {req_text}")
                fig.colorbar(im, ax=ax, label="Number of doublets", pad=0.01)
                pdf.savefig()
                plt.close()


    def write_denominator_info(self, pdf: PdfPages):
        logger.info(f"Writing efficiency denominator info")
        text = f"Efficiency denominator:"
        function = inspect.getsource(self.get_denominator_mask)
        function = textwrap.dedent(function)
        fig, ax = plt.subplots(figsize=(8, 8))
        args = {"ha":"left", "va":"top", "fontfamily":"monospace"}
        ax.text(0.0, 0.9, text, **args, fontsize=16)
        ax.text(0.0, 0.8, function, **args, fontsize=10)
        ax.axis("off")
        pdf.savefig()
        plt.close()


    def plot_detectable_efficiency_vs_kinematics(self, pdf: PdfPages):

        bins = {
            "mcp_pt": np.linspace(0.0, 10.0, 201),
            "mcp_eta": np.linspace(-0.7, 0.7, 281),
            "mcp_phi": np.linspace(-3.2, 3.2, 321),
        }
        xlabel = {
            "mcp_pt": r"Muon $p_T$ [GeV]",
            "mcp_eta": r"Muon $\eta$",
            "mcp_phi": r"Muon $\phi$ [rad]",
        }

        # denominator
        dmask = self.get_denominator_mask()
        denom = self.mcps[dmask]

        # numerator: the denominator is detectable
        numers = {}
        numers["mcp_detectable_OTB_01"] = self.mcps[dmask & (self.mcps["mcp_detectable_OTB_01"] == True)]
        numers["mcp_detectable_OTB_23"] = self.mcps[dmask & (self.mcps["mcp_detectable_OTB_23"] == True)]
        numers["mcp_detectable_OTB_45"] = self.mcps[dmask & (self.mcps["mcp_detectable_OTB_45"] == True)]
        numers["mcp_detectable_OTB_67"] = self.mcps[dmask & (self.mcps["mcp_detectable_OTB_67"] == True)]
        numers["mcp_detectable_OTB"] = self.mcps[dmask & (self.mcps["mcp_detectable_OTB"] == True)]

        for kin in ["mcp_pt", "mcp_eta", "mcp_phi"]:

            for numer in numers:

                n_denom, edges = np.histogram(denom[kin], bins=bins[kin])
                n_numer, edges = np.histogram(numers[numer][kin], bins=bins[kin])
                efficiency = np.divide(n_numer, n_denom, out=np.zeros_like(n_numer, dtype=float), where=n_denom!=0)
                centers = 0.5 * (edges[1:] + edges[:-1])
                fig, ax = plt.subplots()
                ax.plot(
                    centers,
                    efficiency,
                    marker="o",
                    markersize=1,
                    linestyle="-",
                    color="dodgerblue",
                )
                ax.set_xlabel(xlabel[kin])
                ax.set_ylabel("Detector efficiency")
                ax.set_title(f"{numer}")
                ax.set_ylim(0.7, 1.03)
                pdf.savefig()
                plt.close()


    def plot_doublet_efficiency_vs_kinematics_2(self, pdf: PdfPages):

        bins = {
            "mcp_pt": np.linspace(0.0, 10.0, 201),
            "mcp_eta": np.linspace(-0.7, 0.7, 281),
            "mcp_phi": np.linspace(-3.2, 3.2, 321),
        }
        xlabel = {
            "mcp_pt": r"Muon $p_T$ [GeV]",
            "mcp_eta": r"Muon $\eta$",
            "mcp_phi": r"Muon $\phi$ [rad]",
        }

        # denominator
        dmask = self.get_denominator_mask() & (self.mcps["mcp_detectable_OTB"] == True)
        denom = self.mcps[dmask][["file", "i_event", "i_mcp", "mcp_pt", "mcp_eta", "mcp_phi"]]
        if denom.duplicated().any():
            raise ValueError("Denominator has duplicated rows!")

        # numerator
        doublet_cols = [
            "file",
            "i_event", # the event
            "doublet_system", # the system (IT, OT)
            "doublet_doublelayer", # the double layer
            "doublet_module", # the phi-module
            "doublet_sensor", # the z-sensor
        ]

        # filter doublets to only those with same parent mcp
        pass_cuts = (
            (self.doublets["i_mcp"] != NO_MCP) &
            self.doublets["doublet_ok"] &
            self.doublets["doublet_first_exit"]
        )
        doublets = self.doublets[pass_cuts][ doublet_cols + ["i_mcp"] ].drop_duplicates()

        # check if doublets's [file, i_event, i_mcp] is in denominator
        for kin in ["mcp_pt", "mcp_eta", "mcp_phi"]:
            for ((system, doublelayer), group) in doublets.groupby(["doublet_system",
                                                                    "doublet_doublelayer",
                                                                    ]):
                layers = [doublelayer * 2, doublelayer * 2 + 1]

                doublet_keys = group[["file", "i_event", "i_mcp"]].drop_duplicates()
                merged = denom.merge(doublet_keys, on=["file", "i_event", "i_mcp"], how="inner")

                n_denom, edges = np.histogram(denom[kin], bins=bins[kin])
                n_numer, edges = np.histogram(merged[kin], bins=bins[kin])
                efficiency = np.divide(n_numer, n_denom, out=np.zeros_like(n_numer, dtype=float), where=n_denom!=0)
                centers = 0.5 * (edges[1:] + edges[:-1])
                fig, ax = plt.subplots()
                ax.plot(
                    centers,
                    efficiency,
                    marker="o",
                    markersize=1,
                    linestyle="-",
                    color="dodgerblue",
                )
                ax.set_xlabel(xlabel[kin])
                ax.set_ylabel("Doublet algorithm efficiency")
                ax.set_title(f"{NICKNAMES[system]}, layers {layers}")
                ax.set_ylim(0.7, 1.03)
                pdf.savefig()
                plt.close()


    def plot_doublet_efficiency_vs_kinematics(self, pdf: PdfPages):

        bins = {
            "mcp_pt": np.linspace(0.0, 10.0, 201),
            "mcp_eta": np.linspace(-0.7, 0.7, 281),
            "mcp_phi": np.linspace(-3.2, 3.2, 321),
        }
        xlabel = {
            "mcp_pt": r"Muon $p_T$ [GeV]",
            "mcp_eta": r"Muon $\eta$",
            "mcp_phi": r"Muon $\phi$ [rad]",
        }

        # denominator
        dmask = self.get_denominator_mask()
        denom = self.mcps[dmask][["file", "i_event", "i_mcp", "mcp_pt", "mcp_eta", "mcp_phi"]]
        if denom.duplicated().any():
            raise ValueError("Denominator has duplicated rows!")

        # numerator
        doublet_cols = [
            "file",
            "i_event", # the event
            "doublet_system", # the system (IT, OT)
            "doublet_doublelayer", # the double layer
            "doublet_module", # the phi-module
            "doublet_sensor", # the z-sensor
        ]

        # filter doublets to only those with same parent mcp
        same_parent = self.doublets["i_mcp"] != NO_MCP
        doublets = self.doublets[same_parent][ doublet_cols + ["i_mcp"] ].drop_duplicates()

        # check if doublets's [file, i_event, i_mcp] is in denominator
        for kin in ["mcp_pt", "mcp_eta", "mcp_phi"]:
            for ((system, doublelayer), group) in doublets.groupby(["doublet_system",
                                                                    "doublet_doublelayer",
                                                                    ]):
                layers = [doublelayer * 2, doublelayer * 2 + 1]

                doublet_keys = group[["file", "i_event", "i_mcp"]].drop_duplicates()
                merged = denom.merge(doublet_keys, on=["file", "i_event", "i_mcp"], how="inner")

                n_denom, edges = np.histogram(denom[kin], bins=bins[kin])
                n_numer, edges = np.histogram(merged[kin], bins=bins[kin])
                efficiency = np.divide(n_numer, n_denom, out=np.zeros_like(n_numer, dtype=float), where=n_denom!=0)
                centers = 0.5 * (edges[1:] + edges[:-1])
                fig, ax = plt.subplots()
                ax.plot(
                    centers,
                    efficiency,
                    marker="o",
                    markersize=1,
                    linestyle="-",
                    color="dodgerblue",
                )
                ax.set_xlabel(xlabel[kin])
                ax.set_ylabel("Doublet finding efficiency")
                ax.set_title(f"{NICKNAMES[system]}, layers {layers}")
                ax.set_ylim(0.7, 1.03)
                pdf.savefig()
                plt.close()


    def get_denominator_mask(self):
        mask = (
            (np.abs(self.mcps["mcp_pdg"]) == MUON) &
            (self.mcps["mcp_q"] != 0) &
            (self.mcps["mcp_pt"] > ONE_POINT_FIVE_GEV) &
            (self.mcps["mcp_vertex_r"] < ZERO_POINT_ZERO_ONE_MM) &
            (np.abs(self.mcps["mcp_vertex_z"]) < ZERO_POINT_ZERO_ONE_MM) &
            (np.abs(self.mcps["mcp_eta"]) < BARREL_TRACKER_MAX_ETA)
        )
        return mask


    def plot_doublet_features(self, pdf: PdfPages):
        logger.info("Plotting doublet features ...")
        baseline = self.baseline_doublet_mask() if self.signal else np.ones(len(self.doublets), dtype=bool)

        bins = {
            "doublet_dz": np.linspace(-150, 150, 301) if self.signal else np.linspace(-49e3, 49e3, 101),
            "doublet_dr": np.linspace(0, 1000, 101) if self.signal else np.linspace(0, 1500, 101),
            "doublet_dphi": np.linspace(-1.0, 1.0, 201) if self.signal else np.linspace(-3.2, 3.2, 201),
            "doublet_pt": np.linspace(0, 10, 101),
            "doublet_qoverpt": np.linspace(-0.8, 0.8, 161),
            "doublet_phi_slice": np.linspace(-1, N_LS_PHI_SLICES+1, N_LS_PHI_SLICES+3),
            "mcp_qoverpt": np.linspace(-0.8, 0.8, 161),
            "mc_pt": np.linspace(0, 10, 101),
        }
        xlabel = {
            "doublet_dz": r"dz in rz-plane [mm]",
            "doublet_dr": r"dr in xy-plane [mm]",
            "doublet_dphi": r"dphi in xy-plane [rad]",
            "doublet_pt": r"pT [GeV]",
            "doublet_qoverpt": r"Doublet q/pT [1/GeV]",
            "doublet_phi_slice": r"Phi slice",
            "mcp_qoverpt": r"MC q/pT [1/GeV]",
            "mc_pt": r"MC pT [GeV]",
        }
        formatting = {
            "doublet_dz": ".1f",
            "doublet_dr": ".0f",
            "doublet_dphi": ".3f",
            "doublet_pt": ".1f",
            "doublet_qoverpt": ".3f",
            "doublet_phi_slice": ".0f",
            "mcp_qoverpt": ".3f",
        }

        # 1d histograms
        for feature in [
            "doublet_dz",
            "doublet_dr",
            "doublet_dphi",
            "doublet_pt",
            "doublet_phi_slice",
        ]:

            for semilogy in [
                False,
                # True,
            ]:

                for ((system, doublelayer), group) in self.doublets[baseline].groupby(["doublet_system",
                                                                                       "doublet_doublelayer",
                                                                                       ]):

                        # logger.info(f"Plotting signal doublet feature {feature}, system {system}, doublelayer {doublelayer} ...")
                        layers = [doublelayer * 2, doublelayer * 2 + 1]
                        if len(group) == 0:
                            continue

                        fig, ax = plt.subplots()
                        ax.hist(
                            group[feature],
                            bins=bins[feature],
                            histtype="stepfilled",
                            color="crimson",
                            edgecolor="black",
                            linewidth=1.0,
                            alpha=0.9,
                        )
                        if semilogy:
                            ax.semilogy()
                        num = len(group)
                        mean = np.mean(group[feature])
                        rms = np.sqrt(np.mean((group[feature] - mean) ** 2))
                        p997 = np.percentile(np.abs(group[feature]), 99.7)
                        fmt = formatting[feature]
                        ax.set_ylim(0.8 if semilogy else 0, None)
                        ax.set_xlabel(xlabel[feature])
                        ax.set_ylabel("Doublets")
                        ax.set_title(f"{NICKNAMES[system]} layers {layers}. N={num}, Mean={mean:{fmt}}, RMS={rms:{fmt}}")
                        ax.text(0.05, 0.95, f"99.7% in {p997:{fmt}}", transform=ax.transAxes)
                        logger.info(f"{NICKNAMES[system]} doublelayer {doublelayer} {feature}: 99.7% in {p997:{fmt}}")
                        pdf.savefig()
                        plt.close()


        # 2d histograms
        for feature_x, feature_y in [
            # ("doublet_dphi", "doublet_dr"),
            # ("doublet_dphi", "mcp_qoverpt"),
            # ("doublet_qoverpt", "mcp_qoverpt"),
        ]:

            if not self.signal and any(["mcp" in feat for feat in [feature_x, feature_y]]):
                continue

            for ((system, doublelayer), group) in self.doublets[baseline].groupby(["doublet_system",
                                                                                   "doublet_doublelayer",
                                                                                    ]):

                logger.info(f"Plotting signal doublet features {feature_x} vs {feature_y}, system {system}, doublelayer {doublelayer} ...")
                layers = [doublelayer * 2, doublelayer * 2 + 1]
                if len(group) == 0:
                    continue

                fig, ax = plt.subplots()
                h2d, _, _, im = ax.hist2d(
                    group[feature_x],
                    group[feature_y],
                    bins=[bins[feature_x], bins[feature_y]],
                    cmap="gist_rainbow",
                    cmin=0.5,
                )
                if np.nansum(h2d) == 0:
                    raise ValueError(f"No entries in 2d histogram. Fix the binning!")
                fig.colorbar(im, ax=ax, label="Doublets", pad=0.01)
                num = len(group)
                ax.set_xlabel(xlabel[feature_x])
                ax.set_ylabel(xlabel[feature_y])
                ax.set_title(f"{NICKNAMES[system]} layers {layers}. N={num}")
                pdf.savefig()
                plt.close()


    def baseline_doublet_mask(self) -> pd.Series:
        return (
            (self.doublets["i_mcp"] != NO_MCP) &
            (self.doublets["doublet_first_exit"]) &
            (np.abs(self.doublets["mcp_pdg"]) == MUON) &
            (self.doublets["mcp_q"] != 0) &
            (self.doublets["mcp_pt"] > ONE_POINT_FIVE_GEV) &
            (np.abs(self.doublets["mcp_eta"]) < BARREL_TRACKER_MAX_ETA) &
            (self.doublets["mcp_vertex_r"] < ZERO_POINT_ZERO_ONE_MM) &
            (np.abs(self.doublets["mcp_vertex_z"]) < ZERO_POINT_ZERO_ONE_MM) &
            np.ones(len(self.doublets), dtype=bool)
        )


    def plot_doublet_quality_efficiency(self, pdf: PdfPages):

        bins = {
            "mcp_pt": np.linspace(0.0, 10.0, 201),
            "mcp_eta": np.linspace(-0.7, 0.7, 281),
            "mcp_phi": np.linspace(-3.2, 3.2, 321),
        }
        xlabel = {
            "mcp_pt": r"Muon $p_T$ [GeV]",
            "mcp_eta": r"Muon $\eta$",
            "mcp_phi": r"Muon $\phi$ [rad]",
        }

        # only consider truth-match doublets
        baseline = self.baseline_doublet_mask()
        logger.info(f"Doublet efficiency: total doublets: {len(self.doublets)}")
        logger.info(f"Doublet efficiency: total doublets in baseline: {baseline.sum()}")

        # todo: add comment
        for i_kin, kin in enumerate([
            "mcp_pt",
            "mcp_eta",
            "mcp_phi"
        ]):

            for ((system, doublelayer), group) in self.doublets[baseline].groupby(["doublet_system",
                                                                                   "doublet_doublelayer",
                                                                                   ]):

                logger.info(f"Plotting doublet quality efficiency vs {kin}, system {system}, doublelayer {doublelayer} ...")
                layers = [doublelayer * 2, doublelayer * 2 + 1]

                for req in DOUBLET_REQS:
                    req_text, req_mask = self.doublet_requirements(group, req)
                    denom = group
                    numer = group[req_mask]
                    if i_kin == 0:
                        logger.info(f"Denom for system {system} layers {layers} {req}: {len(denom)} doublets")
                        logger.info(f"Numer for system {system} layers {layers} {req}: {len(numer)} doublets")

                    n_denom, edges = np.histogram(denom[kin], bins=bins[kin])
                    n_numer, edges = np.histogram(numer[kin], bins=bins[kin])
                    efficiency = np.divide(n_numer, n_denom, out=np.zeros_like(n_numer, dtype=float), where=n_denom!=0)
                    centers = 0.5 * (edges[1:] + edges[:-1])

                    fig, ax = plt.subplots()
                    ax.plot(
                        centers,
                        efficiency,
                        marker="o",
                        markersize=1,
                        linestyle="-",
                        color="dodgerblue",
                    )
                    ax.set_xlabel(xlabel[kin])
                    ax.set_ylabel("Doublet quality efficiency")
                    ax.set_title(f"{NICKNAMES[system]} layers {layers}: {req_text}")
                    ax.set_ylim(0.965, 1.004)
                    pdf.savefig()
                    plt.close()


    def write_doublet_denominator_info(self, pdf: PdfPages):
        text = f"Double quality efficiency denominator:"
        fig, ax = plt.subplots(figsize=(8, 8))
        args = {"ha":"left", "va":"top", "fontfamily":"monospace"}
        ax.text(-0.1, 0.9, text, **args, fontsize=16)
        ax.text(-0.1, 0.4, f"MIN_COSTHETA = {MIN_COSTHETA}")
        ax.text(-0.1, 0.3, f"MAX_TIME = {MAX_TIME} ns")
        ax.text(-0.1, 0.2, f"MIN_SIMHIT_PT_FRACTION = {MIN_SIMHIT_PT_FRACTION}")
        ax.axis("off")
        pdf.savefig()
        plt.close()


    def baseline_linesegment_mask(self) -> pd.Series:
        return (
            (self.linesegments["i_mcp"] != NO_MCP) &
            (self.linesegments["ls_first_exit"]) &
            (np.abs(self.linesegments["mcp_pdg"]) == MUON) &
            (self.linesegments["mcp_q"] != 0) &
            (self.linesegments["mcp_pt"] > ONE_POINT_FIVE_GEV) &
            (np.abs(self.linesegments["mcp_eta"]) < BARREL_TRACKER_MAX_ETA) &
            (self.linesegments["mcp_vertex_r"] < ZERO_POINT_ZERO_ONE_MM) &
            (np.abs(self.linesegments["mcp_vertex_z"]) < ZERO_POINT_ZERO_ONE_MM) &
            (self.linesegments["ls_md_ok_lower"]) &
            (self.linesegments["ls_md_ok_upper"]) &
            np.ones(len(self.linesegments), dtype=bool)
        )


    def plot_linesegment_features(self, pdf: PdfPages):
        logger.info("Plotting linesegment features ...")
        baseline = self.baseline_linesegment_mask() if self.signal else np.ones(len(self.linesegments), dtype=bool)

        bins = {
            "ls_deta": np.linspace(-3.2, 3.2, 641) if not self.signal else np.linspace(-0.012, 0.012, 241),
            "ls_dphi": np.linspace(-3.2, 3.2, 321) if not self.signal else np.linspace(-0.12, 0.12, 241),
            "ls_dr": np.linspace(0, 1500, 501) if not self.signal else np.linspace(0, 1000, 401),
            "ls_dz": np.linspace(-30000, 30000, 201) if not self.signal else np.linspace(-200, 200, 201),
            "ls_ddr": np.linspace(-300, 300, 601),
            "ls_ddz": np.linspace(-60, 60, 601),
            "ls_dqoverpt": np.linspace(-0.2, 0.2, 201),
            "ls_dtheta_rz": np.linspace(-0.024, 0.024, 241),
            "ls_dtheta_xy": np.linspace(-0.12, 0.12, 241),
            # "ls_chi2_012": np.linspace(0, 2.0, 201),
            "ls_chi2_012": np.linspace(0, 0.01, 201),
        }
        xlabel = {
            "ls_deta": r"upper doublet eta - lower doublet eta",
            "ls_dphi": r"upper doublet phi - lower doublet phi [rad]",
            "ls_dr": "line segment dr [mm]",
            "ls_dz": "line segment dz [mm]",
            "ls_ddr": "upper doublet dr - lower doublet dr",
            "ls_ddz": "upper doublet dz - lower doublet dz",
            "ls_dqoverpt": "upper doublet q/pt - lower doublet q/pt",
            "ls_dtheta_rz": "upper doublet theta_rz - lower doublet theta_rz",
            "ls_dtheta_xy": "upper doublet theta_xy - lower doublet theta_xy",
            "ls_chi2_012": "Diff^2 between circle(xy, 012) and 3 [mm^2]",
        }
        formatting = {
            "ls_deta": ".5f",
            "ls_dphi": ".5f",
            "ls_dr": ".1f",
            "ls_dz": ".1f",
            "ls_ddr": ".3f",
            "ls_ddz": ".3f",
            "ls_dqoverpt": ".3f",
            "ls_dtheta_rz": ".5f",
            "ls_dtheta_xy": ".4f",
            "ls_chi2_012": ".5f",
        }
        color = "cornflowerblue" if self.signal else "crimson"

        # 1d histograms
        for feature in [
            "ls_deta",
            "ls_dphi",
            "ls_ddr",
            "ls_ddz",
            "ls_dqoverpt",
            "ls_dr",
            "ls_dz",
            "ls_dtheta_rz",
            "ls_dtheta_xy",
            "ls_chi2_012",
        ]:

            for semilogy in [
                False,
                # True,
            ]:

                for ((system, doublelayer), group) in self.linesegments[baseline].groupby(["ls_system",
                                                                                           "ls_doublelayer",
                                                                                           ]):

                        # logger.info(f"Plotting signal linesegment feature {feature}, system {system}, doublelayer {doublelayer} ...")

                        fig, ax = plt.subplots()
                        ax.hist(
                            group[feature],
                            bins=bins[feature],
                            histtype="stepfilled",
                            color=color,
                            edgecolor="black",
                            linewidth=1.0,
                            alpha=0.9,
                        )
                        if semilogy or feature in ["ls_chi2_012"]:
                            ax.semilogy()
                        num = len(group)
                        mean = np.mean(group[feature])
                        rms = np.sqrt(np.mean((group[feature] - mean) ** 2))
                        p997 = np.percentile(np.abs(group[feature]), 99.7)
                        fmt = formatting[feature]
                        ax.set_ylim(0.8 if ax.get_yscale() == "log" else 0, None)
                        ax.set_xlabel(xlabel[feature])
                        ax.set_ylabel("Line Segments")
                        ax.set_title(f"{NICKNAMES[system]}. DL={doublelayer}. N={num}, Mean={mean:{fmt}}, RMS={rms:{fmt}}")
                        ax.text(0.30, 0.92, f"99.7% in {p997:{fmt}}", transform=ax.transAxes, fontsize=16)
                        logger.info(f"{NICKNAMES[system]} doublelayer {doublelayer} {feature}: 99.7% in {p997:{fmt}}")
                        pdf.savefig()
                        plt.close()

        # 2d histograms
        for feature_x, feature_y in [
            # ("ls_dphi", "ls_dr"),
        ]:

            for ((system, doublelayer), group) in self.linesegments[baseline].groupby(["ls_system",
                                                                                       "ls_doublelayer",
                                                                                       ]):

                logger.info(f"Plotting signal linesegment features {feature_x} vs {feature_y}, system {system}, doublelayer {doublelayer} ...")
                if len(group) == 0:
                    logger.info(f"No linesegments in {NICKNAMES[system]} doublelayer {doublelayer} passing baseline, skipping feature plot")
                    continue

                fig, ax = plt.subplots()
                _, _, _, im = ax.hist2d(
                    group[feature_x],
                    group[feature_y],
                    bins=[bins[feature_x], bins[feature_y]],
                    cmap="gist_rainbow",
                    cmin=0.5,
                )
                fig.colorbar(im, ax=ax, label="Line Segments", pad=0.01)
                num = len(group)
                ax.set_xlabel(xlabel[feature_x])
                ax.set_ylabel(xlabel[feature_y])
                ax.set_title(f"{NICKNAMES[system]} DL={doublelayer}. N={num}")
                pdf.savefig()
                plt.close()


    def plot_segment_efficiency_vs_kinematics(self, pdf: PdfPages):

        bins = {
            "mcp_pt": np.linspace(0.0, 10.0, 201),
            "mcp_eta": np.linspace(-0.7, 0.7, 281),
            "mcp_phi": np.linspace(-3.2, 3.2, 321),
        }
        xlabel = {
            "mcp_pt": r"Muon $p_T$ [GeV]",
            "mcp_eta": r"Muon $\eta$",
            "mcp_phi": r"Muon $\phi$ [rad]",
        }

        # denominator
        dmask = self.get_denominator_mask()
        denom = self.mcps[dmask][["file", "i_event", "i_mcp", "mcp_pt", "mcp_eta", "mcp_phi"]]
        if denom.duplicated().any():
            raise ValueError("Denominator has duplicated rows!")

        # numerator
        segment_cols = [
            "file", # the file
            "i_event", # the event
            "i_mcp", # the parent mc particle
            "ls_system", # the system (IT, OT)
            "ls_doublelayer", # the first double layer
        ]

        # filter doublets to only those with same parent mcp
        same_parent = self.linesegments["i_mcp"] != NO_MCP
        segments = self.linesegments[same_parent][segment_cols].drop_duplicates()

        # check if segments's [file, i_event, i_mcp] is in denominator
        for kin in ["mcp_pt", "mcp_eta", "mcp_phi"]:
            for ((system, doublelayer), group) in segments.groupby(["ls_system",
                                                                    "ls_doublelayer",
                                                                    ]):
                layer = doublelayer * 2
                layers = range(layer, layer + 4)

                segment_keys = group[["file", "i_event", "i_mcp"]].drop_duplicates()
                merged = denom.merge(segment_keys, on=["file", "i_event", "i_mcp"], how="inner")

                n_denom, edges = np.histogram(denom[kin], bins=bins[kin])
                n_numer, edges = np.histogram(merged[kin], bins=bins[kin])
                efficiency = np.divide(n_numer, n_denom, out=np.zeros_like(n_numer, dtype=float), where=n_denom!=0)
                centers = 0.5 * (edges[1:] + edges[:-1])
                fig, ax = plt.subplots()
                ax.plot(
                    centers,
                    efficiency,
                    marker="o",
                    markersize=1,
                    linestyle="-",
                    color="dodgerblue",
                )
                ax.set_xlabel(xlabel[kin])
                ax.set_ylabel("T2 finding efficiency")
                ax.set_title(f"{NICKNAMES[system]}, layers {list(layers)}")
                ax.set_ylim(0.7, 1.03)
                pdf.savefig()
                plt.close()


    def plot_segment_quality_efficiency(self, pdf: PdfPages):

        bins = {
            "mcp_pt": np.linspace(0.0, 10.0, 201),
            "mcp_eta": np.linspace(-0.7, 0.7, 281),
            "mcp_phi": np.linspace(-3.2, 3.2, 321),
        }
        xlabel = {
            "mcp_pt": r"Muon $p_T$ [GeV]",
            "mcp_eta": r"Muon $\eta$",
            "mcp_phi": r"Muon $\phi$ [rad]",
        }

        # only consider truth-match doublets
        baseline = self.baseline_linesegment_mask()
        logger.info(f"Segment efficiency: total segments: {len(self.linesegments)}")
        logger.info(f"Segment efficiency: total segments in baseline: {baseline.sum()}")

        # consider efficiency vs kinematics
        for i_kin, kin in enumerate([
            "mcp_pt",
            "mcp_eta",
            "mcp_phi"
        ]):

            for ((system, doublelayer), group) in self.linesegments[baseline].groupby(["ls_system",
                                                                                       "ls_doublelayer",
            ]):

                logger.info(f"Plotting segment quality efficiency vs {kin}, system {system}, doublelayer {doublelayer} ...")
                layer = doublelayer * 2
                layers = range(layer, layer + 4)

                for req in LS_REQS:
                    req_text, req_mask = self.segment_requirements(group, req)
                    denom = group
                    numer = group[req_mask]
                    if i_kin == 0:
                        logger.info(f"Denom for system {system} layers {layers} {req}: {len(denom)} doublets")
                        logger.info(f"Numer for system {system} layers {layers} {req}: {len(numer)} doublets")

                    n_denom, edges = np.histogram(denom[kin], bins=bins[kin])
                    n_numer, edges = np.histogram(numer[kin], bins=bins[kin])
                    efficiency = np.divide(n_numer, n_denom, out=np.zeros_like(n_numer, dtype=float), where=n_denom!=0)
                    centers = 0.5 * (edges[1:] + edges[:-1])

                    fig, ax = plt.subplots()
                    ax.plot(
                        centers,
                        efficiency,
                        marker="o",
                        markersize=1,
                        linestyle="-",
                        color="dodgerblue",
                    )
                    ax.set_xlabel(xlabel[kin])
                    ax.set_ylabel("Segment quality efficiency")
                    ax.set_title(f"{NICKNAMES[system]} layers {layers}: {req_text}")
                    ax.set_ylim(0.965, 1.004)
                    pdf.savefig()
                    plt.close()


    def segment_requirements(self, df: pd.DataFrame, req: str) -> tuple[str, pd.DataFrame]:
        # return description and mask
        doublelayer = df["ls_doublelayer"]
        if len(doublelayer.unique()) != 1:
            raise ValueError(f"Multiple doublelayers found: {doublelayer.unique()}")
        doublelayer = doublelayer.iloc[0]
        if req == REQ_PASSTHROUGH:
            text = "No requirement"
            mask = np.ones(len(df), dtype=bool)
        elif req == LS_REQ_DR_POS:
            text = f"|dr| < {self.LS_DR_CUT[doublelayer]}mm"
            mask = np.abs(df["ls_dr"]) < self.LS_DR_CUT[doublelayer]
        elif req == LS_REQ_DZ_POS:
            text = f"|dz| < {self.LS_DZ_CUT[doublelayer]}mm"
            mask = np.abs(df["ls_dz"]) < self.LS_DZ_CUT[doublelayer]
        elif req == LS_REQ_RZ_ANG:
            text = f"|dtheta(rz)| < {self.LS_DTHETA_RZ_CUT[doublelayer]}rad"
            mask = np.abs(df["ls_dtheta_rz"]) < self.LS_DTHETA_RZ_CUT[doublelayer]
        elif req == LS_REQ_XY_CHI2:
            text = f"Chi2(xy,012) < {self.LS_CHI2_XY_CUT[doublelayer]}"
            mask = np.abs(df["ls_chi2_012"]) < self.LS_CHI2_XY_CUT[doublelayer]
        elif req == LS_REQ_ALL:
            text = f"All LS requirements"
            mask = df["ls_ok"]
        else:
            raise ValueError(f"Unknown segment requirement: {req}")
        return text, mask


    def baseline_t4_mask(self) -> pd.Series:
        return (
            (self.t4s["i_mcp"] != NO_MCP) &
            (self.t4s["t4_first_exit"]) &
            (np.abs(self.t4s["mcp_pdg"]) == MUON) &
            (self.t4s["mcp_q"] != 0) &
            (self.t4s["mcp_pt"] > ONE_POINT_FIVE_GEV) &
            (np.abs(self.t4s["mcp_eta"]) < BARREL_TRACKER_MAX_ETA) &
            (self.t4s["mcp_vertex_r"] < ZERO_POINT_ZERO_ONE_MM) &
            (np.abs(self.t4s["mcp_vertex_z"]) < ZERO_POINT_ZERO_ONE_MM) &
            (self.t4s["t4_ls_ok_lower"]) &
            (self.t4s["t4_ls_ok_upper"]) &
            np.ones(len(self.t4s), dtype=bool)
        )


    def plot_t4_features(self, pdf: PdfPages):
        logger.info("Plotting t4 features ...")
        if self.t4s is None or len(self.t4s) == 0:
            logger.info("No T4s to plot")
            return
        baseline = self.baseline_t4_mask() if self.signal else np.ones(len(self.t4s), dtype=bool)

        bins = {
            "t4_deta": np.linspace(-3.2, 3.2, 641) if not self.signal else np.linspace(-0.032, 0.032, 321),
            "t4_dphi": np.linspace(-3.2, 3.2, 321) if not self.signal else np.linspace(-0.32, 0.32, 321),
            "t4_dr": np.linspace(0, 1500, 501) if not self.signal else np.linspace(0, 1000, 401),
            "t4_dz": np.linspace(-30000, 30000, 201) if not self.signal else np.linspace(-200, 200, 201),
            "t4_dtheta_rz": np.linspace(-0.08, 0.08, 241),
            "t4_chi2_047": np.linspace(0, 0.5, 201),
        }
        xlabel = {
            "t4_deta": "upper T2 eta - lower T2 eta",
            "t4_dphi": "upper T2 phi - lower T2 phi [rad]",
            "t4_dr": "T4 dr [mm]",
            "t4_dz": "T4 dz [mm]",
            "t4_dtheta_rz": "upper T2 theta_rz - lower T2 theta_rz",
            "t4_chi2_047": "Diff^2 between circle(xy, 047) and 12356 [mm^2]",
        }
        formatting = {
            "t4_deta": ".5f",
            "t4_dphi": ".5f",
            "t4_dr": ".1f",
            "t4_dz": ".1f",
            "t4_dtheta_rz": ".5f",
            "t4_chi2_047": ".5f",
        }
        color = "cornflowerblue" if self.signal else "crimson"

        # 1d histograms
        for feature in [
            "t4_deta",
            "t4_dphi",
            "t4_dr",
            "t4_dz",
            "t4_dtheta_rz",
            "t4_chi2_047",
        ]:

            for semilogy in [
                False,
                # True,
            ]:

                groupby_cols = ["t4_system", "t4_doublelayer"]
                for ((system, doublelayer), group) in self.t4s[baseline].groupby(groupby_cols):

                        # logger.info(f"Plotting signal linesegment feature {feature}, system {system}, doublelayer {doublelayer} ...")

                        fig, ax = plt.subplots()
                        ax.hist(
                            group[feature],
                            bins=bins[feature],
                            histtype="stepfilled",
                            color=color,
                            edgecolor="black",
                            linewidth=1.0,
                            alpha=0.9,
                        )
                        if semilogy or feature in ["t4_chi2_047"]:
                            ax.semilogy()
                        num = len(group)
                        mean = np.mean(group[feature])
                        rms = np.sqrt(np.mean((group[feature] - mean) ** 2))
                        p997 = np.percentile(np.abs(group[feature]), 99.7)
                        fmt = formatting[feature]
                        ax.set_ylim(0.8 if ax.get_yscale() == "log" else 0, None)
                        ax.set_xlabel(xlabel[feature])
                        ax.set_ylabel("T4s")
                        ax.set_title(f"{NICKNAMES[system]}. DL={doublelayer}. N={num}, Mean={mean:{fmt}}, RMS={rms:{fmt}}")
                        ax.text(0.30, 0.92, f"99.7% in {p997:{fmt}}", transform=ax.transAxes, fontsize=16)
                        logger.info(f"{NICKNAMES[system]} doublelayer {doublelayer} {feature}: 99.7% in {p997:{fmt}}")
                        pdf.savefig()
                        plt.close()


    def plot_t4_efficiency_vs_kinematics(self, pdf: PdfPages):

        if self.t4s is None or len(self.t4s) == 0:
            logger.info("No T4s to plot")
            return

        bins = {
            "mcp_pt": np.linspace(0.0, 10.0, 201),
            "mcp_eta": np.linspace(-0.7, 0.7, 281),
            "mcp_phi": np.linspace(-3.2, 3.2, 321),
        }
        xlabel = {
            "mcp_pt": r"Muon $p_T$ [GeV]",
            "mcp_eta": r"Muon $\eta$",
            "mcp_phi": r"Muon $\phi$ [rad]",
        }

        # denominator
        dmask = self.get_denominator_mask()
        denom = self.mcps[dmask][["file", "i_event", "i_mcp", "mcp_pt", "mcp_eta", "mcp_phi"]]
        if denom.duplicated().any():
            raise ValueError("Denominator has duplicated rows!")

        # numerator
        numer_cols = [
            "file", # the file
            "i_event", # the event
            "i_mcp", # the parent mc particle
            "t4_system", # the system (IT, OT)
            "t4_doublelayer", # the first double layer
        ]

        # filter doublets to only those with same parent mcp
        same_parent = self.t4s["i_mcp"] != NO_MCP
        t4s = self.t4s[same_parent][numer_cols].drop_duplicates()

        # check if t4s's [file, i_event, i_mcp] is in denominator
        for kin in ["mcp_pt", "mcp_eta", "mcp_phi"]:
            for ((system, doublelayer), group) in t4s.groupby(["t4_system",
                                                               "t4_doublelayer",
                                                               ]):
                layer = doublelayer * 2
                layers = range(layer, layer + 8)

                keys = group[["file", "i_event", "i_mcp"]].drop_duplicates()
                merged = denom.merge(keys, on=["file", "i_event", "i_mcp"], how="inner")

                n_denom, edges = np.histogram(denom[kin], bins=bins[kin])
                n_numer, edges = np.histogram(merged[kin], bins=bins[kin])
                efficiency = np.divide(n_numer, n_denom, out=np.zeros_like(n_numer, dtype=float), where=n_denom!=0)
                centers = 0.5 * (edges[1:] + edges[:-1])
                fig, ax = plt.subplots()
                ax.plot(
                    centers,
                    efficiency,
                    marker="o",
                    markersize=1,
                    linestyle="-",
                    color="dodgerblue",
                )
                ax.set_xlabel(xlabel[kin])
                ax.set_ylabel("T4 finding efficiency")
                ax.set_title(f"{NICKNAMES[system]}, layers {list(layers)}")
                ax.set_ylim(0.7, 1.03)
                pdf.savefig()
                plt.close()


    def plot_t4_quality_efficiency(self, pdf: PdfPages):

        if self.t4s is None or len(self.t4s) == 0:
            logger.info("No T4s to plot")
            return

        bins = {
            "mcp_pt": np.linspace(0.0, 10.0, 201),
            "mcp_eta": np.linspace(-0.7, 0.7, 281),
            "mcp_phi": np.linspace(-3.2, 3.2, 321),
        }
        xlabel = {
            "mcp_pt": r"Muon $p_T$ [GeV]",
            "mcp_eta": r"Muon $\eta$",
            "mcp_phi": r"Muon $\phi$ [rad]",
        }

        # the cuts
        reqs = [col for col in self.t4s.columns if col.startswith("t4_ok")]

        # only consider truth-matched
        baseline = self.baseline_t4_mask()
        logger.info(f"T4 efficiency: total T4: {len(self.t4s)}")
        logger.info(f"T4 efficiency: total T4 in baseline: {baseline.sum()}")

        # consider efficiency vs kinematics
        for i_kin, kin in enumerate([
            "mcp_pt",
            "mcp_eta",
            "mcp_phi"
        ]):

            groupbys = ["t4_system", "t4_doublelayer"]
            for ((system, doublelayer), group) in self.t4s[baseline].groupby(groupbys):

                logger.info(f"Plotting T4 quality efficiency vs {kin}, system {system}, doublelayer {doublelayer} ...")
                layer = doublelayer * 2
                layers = range(layer, layer + 4)

                for req in reqs:
                    denom = group
                    numer = group[ group[req] ]
                    if i_kin == 0:
                        logger.info(f"Denom for system {system} layers {layers} {req}: {len(denom)} doublets")
                        logger.info(f"Numer for system {system} layers {layers} {req}: {len(numer)} doublets")

                    n_denom, edges = np.histogram(denom[kin], bins=bins[kin])
                    n_numer, edges = np.histogram(numer[kin], bins=bins[kin])
                    efficiency = np.divide(n_numer, n_denom, out=np.zeros_like(n_numer, dtype=float), where=n_denom!=0)
                    centers = 0.5 * (edges[1:] + edges[:-1])

                    fig, ax = plt.subplots()
                    ax.plot(
                        centers,
                        efficiency,
                        marker="o",
                        markersize=1,
                        linestyle="-",
                        color="dodgerblue",
                    )
                    ax.set_xlabel(xlabel[kin])
                    ax.set_ylabel("T4 quality efficiency")
                    ax.set_title(f"{NICKNAMES[system]} layers {layers}: {req}")
                    ax.set_ylim(0.965, 1.004)
                    pdf.savefig()
                    plt.close()

