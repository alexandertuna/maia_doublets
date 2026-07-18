import argparse
from glob import glob
import numpy as np
import pandas as pd
# import datashader as ds
import matplotlib.pyplot as plt
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

GAP = 2.0
RADII = [
    127, 167, 510, 550,
    819, 899, 1366, 1446,
]
AX_IDX = {
    "hits": (0, 1),
    "mds": (0, 2),
    "t2s": (1, 0),
    "t4s": (1, 1),
    "t8s": (1, 2),
}
RADII = [rad + GAP for rad in RADII]
R_MAX = max(RADII) * 1.1
N_LAYERS_IN_MDS = 2
N_LAYERS_IN_T2S = 4
N_LAYERS_IN_T4S = 8
N_LAYERS_IN_T8S = 16


def main():
    logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")
    args = options()
    check_options(args)
    scatter = HitsScatter2d(
        input_files_hits=get_input_files(args.hits),
        input_files_mds=get_input_files(args.mds),
        input_files_t2s=get_input_files(args.t2s),
        input_files_t4s=get_input_files(args.t4s),
        input_files_t8s=get_input_files(args.t8s),
        output_file=args.output,
        n_random=args.n,
    )
    scatter.plot()


def options():
    parser = argparse.ArgumentParser(description="Demo for 2D scatter plot")
    parser.add_argument("--hits", type=str, default="", help="Comma-separated list of input hits file paths")
    parser.add_argument("--mds", type=str, default="", help="Comma-separated list of input mds file paths")
    parser.add_argument("--t2s", type=str, default="", help="Comma-separated list of input t2s file paths")
    parser.add_argument("--t4s", type=str, default="", help="Comma-separated list of input t4s file paths")
    parser.add_argument("--t8s", type=str, default="", help="Comma-separated list of input t8s file paths")
    parser.add_argument("-n", type=int, help="Sample n random points to speed things up")
    parser.add_argument("--output", type=str, help="Output pdf file path")
    return parser.parse_args()


def check_options(args):
    if not args.hits:
        raise ValueError("Hits input file(s) must be specified")
    if not args.mds:
        raise ValueError("MDs input file(s) must be specified")
    if not args.output:
        raise ValueError("Output file must be specified")


def get_input_files(input_str):
    fnames = []
    for pattern in input_str.split(","):
        fnames.extend(glob(pattern))
    return fnames


class HitsScatter2d:

    def __init__(self,
                 input_files_hits: list[str],
                 input_files_mds: list[str],
                 input_files_t2s: list[str],
                 input_files_t4s: list[str],
                 input_files_t8s: list[str],
                 output_file: str,
                 n_random: int,
                 ):
        self.input_files_hits = input_files_hits
        self.input_files_mds = input_files_mds
        self.input_files_t2s = input_files_t2s
        self.input_files_t4s = input_files_t4s
        self.input_files_t8s = input_files_t8s
        self.objs = ["hits", "mds", "t2s", "t4s", "t8s"]
        if len(self.input_files_hits) == 0:
            self.objs.remove("hits")
        if len(self.input_files_mds) == 0:
            self.objs.remove("mds")
        if len(self.input_files_t2s) == 0:
            self.objs.remove("t2s")
        if len(self.input_files_t4s) == 0:
            self.objs.remove("t4s")
        if len(self.input_files_t8s) == 0:
            self.objs.remove("t8s")
        self.output_file = output_file
        self.hits_cols = {
            "simhit_x": "x",
            "simhit_y": "y",
            "simhit_z": "z",
        }
        self.n_random = n_random
        self.df = {}
        self.get_input_data()


    def get_input_data(self):
        self.df["hits"] = self.get_input_hits()
        self.df["mds"] = self.get_input_mds()
        self.df["t2s"] = self.get_input_t2s()
        self.df["t4s"] = self.get_input_t4s()
        self.df["t8s"] = self.get_input_t8s()


    def get_input_hits(self) -> pd.DataFrame:
        if len(self.input_files_hits) == 0:
            return pd.DataFrame()
        tmp = pd.concat([pd.read_pickle(fi) for fi in self.input_files_hits], ignore_index=True)
        df = tmp[list(self.hits_cols.keys())].rename(columns=self.hits_cols)
        df = self.sample_random(df)
        logger.info(f"Total hits read: {len(df)}")
        return df


    def get_input_mds(self) -> pd.DataFrame:
        if len(self.input_files_mds) == 0:
            return pd.DataFrame()
        tmp = pd.concat([pd.read_pickle(fi) for fi in self.input_files_mds], ignore_index=True)
        df = pd.concat([
            tmp[[f"doublet_x_{it}", f"doublet_y_{it}"]].rename(columns={f"doublet_x_{it}": "x", f"doublet_y_{it}": "y"})
            for it in range(N_LAYERS_IN_MDS)
        ])
        df = self.sample_random(df)
        logger.info(f"Total mds read: {len(df)}")
        return df


    def get_input_t2s(self) -> pd.DataFrame:
        if len(self.input_files_t2s) == 0:
            return pd.DataFrame()
        tmp = pd.concat([pd.read_pickle(fi) for fi in self.input_files_t2s], ignore_index=True)
        df = pd.concat([
            tmp[[f"t2_x_{it}", f"t2_y_{it}"]].rename(columns={f"t2_x_{it}": "x", f"t2_y_{it}": "y"})
            for it in range(N_LAYERS_IN_T2S)
        ])
        df = self.sample_random(df)
        logger.info(f"Total t2s read: {len(df)}")
        return df


    def get_input_t4s(self) -> pd.DataFrame:
        if len(self.input_files_t4s) == 0:
            return pd.DataFrame()
        tmp = pd.concat([pd.read_pickle(fi) for fi in self.input_files_t4s], ignore_index=True)
        df = pd.concat([
            tmp[[f"t4_x_{it}", f"t4_y_{it}"]].rename(columns={f"t4_x_{it}": "x", f"t4_y_{it}": "y"})
            for it in range(N_LAYERS_IN_T4S)
        ])
        df = self.sample_random(df)
        logger.info(f"Total t4s read: {len(df)}")
        return df


    def get_input_t8s(self) -> pd.DataFrame:
        if len(self.input_files_t8s) == 0:
            return pd.DataFrame()
        tmp = pd.concat([pd.read_pickle(fi) for fi in self.input_files_t8s], ignore_index=True)
        df = pd.concat([
            tmp[[f"t8_x_{it}", f"t8_y_{it}"]].rename(columns={f"t8_x_{it}": "x", f"t8_y_{it}": "y"})
            for it in range(N_LAYERS_IN_T8S)
        ])
        df = self.sample_random(df)
        logger.info(f"Total t8s read: {len(df)}")
        return df


    def sample_random(self, df: pd.DataFrame) -> pd.DataFrame:
        if self.n_random and self.n_random < len(df):
            logger.info(f"Sampling {self.n_random} random hits from {len(df)} total hits")
            random_n_indices = np.random.choice(df.index, size=self.n_random, replace=False)
            df = df.loc[random_n_indices]
        return df


    def plot(self):

        cargs = {
            "facecolor": None,
            "fill": False,
            "edgecolor": "black",
            "linestyle": "-",
            "linewidth": 0.5,
            "zorder": 0,
        }
        sargs = {
            "s": 10,
            "c": "blue",
            "alpha": 0.5,
            "edgecolors": "none",
        }

        # basic 2d
        logger.info(f"Making scatter plot ... ")
        fig, ax = plt.subplots(nrows=2, ncols=3)
        for obj in self.objs:
            logger.info(f"Plotting {obj} ... ")
            circles = [plt.Circle((0,0), rad, **cargs) for rad in RADII]
            row, col = AX_IDX[obj]
            ax[row, col].scatter(self.df[obj]["x"], self.df[obj]["y"], **sargs)
            for circle in circles:
                ax[row, col].add_patch(circle)
            ax[row, col].set_xlabel("x [mm]")
            ax[row, col].set_ylabel("y [mm]")
            ax[row, col].set_xlim(-R_MAX, R_MAX)
            ax[row, col].set_ylim(-R_MAX, R_MAX)
            ax[row, col].grid(True, alpha=0.3, linewidth=0.5)
        logger.info(f"Saving scatter plot ... ")
        fig.savefig(self.output_file, dpi=500)
        plt.close(fig)

        return

        # fancy 2d
        def project(xyz, tilt_x_deg=-5.0, tilt_y_deg=5.0):
            """Transverse x-y view tilted slightly for depth; z is the into-screen axis."""
            ax_, ay_ = np.radians(tilt_x_deg), np.radians(tilt_y_deg)
            Rx = np.array([[1, 0, 0],
                           [0, np.cos(ax_), -np.sin(ax_)],
                           [0, np.sin(ax_),  np.cos(ax_)]])
            Ry = np.array([[ np.cos(ay_), 0, np.sin(ay_)],
                           [0, 1, 0],
                           [-np.sin(ay_), 0, np.cos(ay_)]])
            p = xyz @ Rx.T @ Ry.T
            return p[:, 0], p[:, 1]        # screen = (x, y); rotated z is dropped

        FIG_KW = dict(figsize=(8, 8), dpi=500)
        LIMS   = dict(xlim=(-1600, 1600), ylim=(-1600, 1600))   # identical every stage
        ALPHA  = 0.2        # tune ONCE on the full-BIB frame
        MS     = 2

        def render(u, v, path, su=None, sv=None):
            fig, ax = plt.subplots(**FIG_KW)
            logger.info(f"Rendering {len(u)} points to {path} ... ")
            ax.plot(u, v, 'o', markersize=MS, alpha=ALPHA, mew=0, color='0.15', ls='none')
            if su is not None:
                ax.plot(su, sv, 'o', markersize=MS, alpha=1.0, mew=0, color='red', ls='none')
            ax.set(aspect='equal', **LIMS)
            ax.axis('off')
            logger.info(f"Saving figure to {path} ... ")
            fig.savefig(path, dpi=FIG_KW['dpi'], pad_inches=0)
            plt.close(fig)

        u, v = project(self.df["hits"][["x", "y", "z"]].to_numpy())
        render(u, v, self.output_file.replace(".png", "_tilted.png"))

        # logger.info(f"Projecting {len(self.hits_df)} hits to 2D for plotting ... ")
        # u, v = self.hits_df["x"], self.hits_df["y"]

        # logger.info(f"Creating ds.Canvas ... ")
        # cvs = ds.Canvas(plot_width=2400, plot_height=1800)

        # logger.info(f"Creating ds.Canvas.points ... ")
        # agg = cvs.points(pd.DataFrame({"u": u, "v": v}), "u", "v")

        # logger.info(f"Plotting {len(agg.values)} aggregated points to PDF ... ")
        # fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
        # ax.imshow(np.log1p(agg.values), origin="lower", cmap="Greys", extent=[u.min(), u.max(), v.min(), v.max()])

        # ax.set_aspect("equal")                       # keep detector proportions
        # ax.set_xlabel("projected u [mm]")
        # ax.set_ylabel("projected v [mm]")
        # ax.grid(True, alpha=0.3, linewidth=0.5)

        # fig.savefig("tmp.pdf", bbox_inches="tight")


        # xs = self.hits_df["x"]
        # ys = self.hits_df["y"]
        # zs = self.hits_df["z"]

        # fig = go.Figure(go.Scatter3d(
        #     x=xs, y=ys, z=zs,
        #     mode="markers",
        #     marker=dict(size=2, color=zs, colorscale="Viridis"),
        # ))
        # if self.output_file.endswith(".html"):
        #     fig.write_html(self.output_file)
        # elif self.output_file.endswith(".pdf"):
        #     fig.write_image(self.output_file)
        # else:
        #     raise ValueError(f"Unsupported output file format: {self.output_file}")



if __name__ == "__main__":
    main()
