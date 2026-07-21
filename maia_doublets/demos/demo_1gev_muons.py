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

import pyLCIO
MCPARTICLES = "MCParticle"
COLLECTIONS = [
    "InnerTrackerBarrelCollection",
    "OuterTrackerBarrelCollection",
]
MUON = 13
MAX_ETA = 0.65
N_PLOTS = 100
FNAME = "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_1p0_1p1/10um/muonGun_pT_1p0_1p1_sim_300.slcio"
RADII = [
    127, 167, 510, 550,
    819, 899, 1366, 1446,
]
R_MAX = max(RADII) * 1.02

def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    df = slcio_to_df(FNAME)
    with PdfPages("one_gev_muons.pdf") as pdf:
        plot_hits_per_mcp(df, pdf)


def plot_hits_per_mcp(df: pd.DataFrame, pdf: PdfPages):
    cargs = {
        "facecolor": None,
        "fill": False,
        "edgecolor": "black",
        "linestyle": "-",
        "linewidth": 0.2,
        "zorder": 0,
    }

    for i_group, (cols, group) in enumerate(df.groupby(["i_event", "i_mcp"])):
        i_event, i_mcp = cols
        n_hits = len(group)
        logger.info(f"Event {i_event}, MCP {i_mcp}: {n_hits} hits")
        fig, ax = plt.subplots()
        circles = [plt.Circle((0,0), rad, **cargs) for rad in RADII]
        for circle in circles:
            ax.add_patch(circle)
        ax.scatter(group["hit_x"],
                   group["hit_y"],
                   s=20,
                   c="red",
                   edgecolors="black",
                   )
        mcp_pt = group["mcp_pt"].iloc[0]
        mcp_eta = group["mcp_eta"].iloc[0]
        mcp_phi = group["mcp_phi"].iloc[0]
        mcp_props = f"pT: {mcp_pt:.2f} GeV\neta: {mcp_eta:.2f}\nphi: {mcp_phi:.2f}"
        ax.text(0.02, 0.98, mcp_props, transform=ax.transAxes, va="top", ha="left")
        ax.set_title(f"Event {i_event}, MCP {i_mcp}, : {n_hits} hits")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_xlim(-R_MAX, R_MAX)
        ax.set_ylim(-R_MAX, R_MAX)
        pdf.savefig(fig)
        plt.close(fig)
        if i_group >= N_PLOTS:
            logger.info(f"Stopping after {N_PLOTS} groups for demo purposes.")
            break


def slcio_to_df(slcio_file_path):

    logger.info(f"Processing file {slcio_file_path} ...")
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(slcio_file_path)

    hits = []
    for i_event, event in enumerate(reader):
        mcps = [mcp for mcp in event.getCollection(MCPARTICLES)]
        for i_col, collection in enumerate(COLLECTIONS):
            col = event.getCollection(collection)
            for i_hit, hit in enumerate(col):

                mcp = hit.getMCParticle()
                i_mcp = mcps.index(mcp) if mcp in mcps else -1
                if i_mcp == -1:
                    continue
                if np.abs(mcp.getPDG()) != MUON:
                    continue
                mcp_p = mcp.getMomentum()
                mcp_px, mcp_py, mcp_pz = mcp_p[0], mcp_p[1], mcp_p[2]
                mcp_pt = np.sqrt(mcp_px**2 + mcp_py**2)
                mcp_theta = np.arctan2(mcp_pt, mcp_pz)
                mcp_eta = -np.log(np.tan(mcp_theta / 2))
                mcp_phi = np.arctan2(mcp_py, mcp_px)
                if np.abs(mcp_eta) > MAX_ETA:
                    continue

                position = hit.getPosition()
                hits.append({
                    "i_event": i_event,
                    "i_col": i_col,
                    "i_hit": i_hit,
                    "i_mcp": i_mcp,
                    "hit_x": position[0],
                    "hit_y": position[1],
                    "hit_z": position[2],
                    "mcp_pt": mcp_pt,
                    "mcp_eta": mcp_eta,
                    "mcp_phi": mcp_phi,
                })

    logger.info(f"Processed {len(hits)} hits in {i_event + 1} events.")
    return pd.DataFrame(hits)


if __name__ == "__main__":
    main()
