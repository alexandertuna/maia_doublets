"""
Equal spaced detector aka v07
"""
import style
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import multiprocessing as mp

MIN_PT = 0.2 # GeV
MIN_T = -10 # ns
MAX_T = 20 # ns
MIN_Z = -4000 # mm
MAX_Z = 4000 # mm
MIN_R = 0
MAX_R = 1900 # mm
MCPARTICLE = "MCParticle"
GEN_STATUS = 1
N_BIB_FILES = 50
BIB = [
    f"/ceph/users/atuna/work/maia/maia_datasets/productions/bib.2026_08_14_17h50m00s/BIB10TeV/sim_{muon}/BIB_sim_{i+1}.slcio"
    for i in range(N_BIB_FILES)
    for muon in ["mp", "mm"]
]
MAX_PROC = 20


PKL = "fig-mechanism.pkl"

def main():

    if os.path.exists(PKL):
        df = pd.read_pickle(PKL)
    else:
        df = get_mcparticles(BIB)
        df.to_pickle(PKL)

    print(df)
    with PdfPages("fig-mechanism.pdf") as pdf:
        plot(df, pdf)


def get_mcparticles(fpaths: list[str]) -> pd.DataFrame:
    dfs = []
    with mp.Pool(processes=MAX_PROC) as pool:
        dfs = pool.map(get_mcparticles_one_file, fpaths)
    return pd.concat(dfs, ignore_index=True)


def get_mcparticles_one_file(fpath: str) -> pd.DataFrame:
    """
    https://github.com/MuonColliderSoft/LCIO/blob/master/src/cpp/include/IMPL/MCParticleImpl.h
    """
    import pyLCIO
    print(f"Processing file {fpath}")
    rows = []
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(fpath)
    for event in reader:
        for mcp in event.getCollection(MCPARTICLE):
            rows.append({
                "pdg": mcp.getPDG(),
                "t": mcp.getTime(),
                "x": mcp.getVertex()[0],
                "y": mcp.getVertex()[1],
                "z": mcp.getVertex()[2],
                "px": mcp.getMomentum()[0],
                "py": mcp.getMomentum()[1],
                "pz": mcp.getMomentum()[2],
                "gen_status": mcp.getGeneratorStatus(),
                "sim_status": mcp.getSimulatorStatus(),
                "nparents": len(mcp.getParents()),
                "nchildren": len(mcp.getDaughters()),
            })
    reader.close()
    return post_process(pd.DataFrame(rows))


def post_process(df: pd.DataFrame) -> pd.DataFrame:
    df["r"] = (df["x"]**2 + df["y"]**2)**0.5
    df["pt"] = (df["px"]**2 + df["py"]**2)**0.5
    df["p"] = (df["px"]**2 + df["py"]**2 + df["pz"]**2)**0.5
    df = df[df["pt"] > MIN_PT]
    df = df[(df["t"] > MIN_T) & (df["t"] < MAX_T)]
    df = df[df["gen_status"] == GEN_STATUS]
    return df


def plot(df: pd.DataFrame, pdf: PdfPages):
    plot_pt(df, pdf)
    plot_t(df, pdf)

    # one vector for each particle (z vs r)
    # vector points in the direction of the particle's momentum
    df["rz"] = np.arctan2(df["r"], df["z"])

    length = 200.0
    x, y, z, r = (df[c].to_numpy() for c in ("x", "y", "z", "r"))
    px, py, pz, pt = (df[c].to_numpy() for c in ("px", "py", "pz", "pt"))
    pr = pt
    norm = np.hypot(pz, pr)
    norm[norm == 0] = np.nan  # skip zero-momentum particles
    u = pz / norm * length
    v = pr / norm * length

    fig, ax = plt.subplots()
    ax.quiver(z, r, u, v,
              angles="xy", scale_units="xy", scale=1,  # arrow length = data units
              width=0.003)
    ax.scatter(z, r, s=5, color="k", zorder=3)  # mark the start points
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_xlim(MIN_Z, MAX_Z)
    ax.set_ylim(MIN_R, MAX_R)
    ax.set_aspect("equal")  # otherwise arrow angles look distorted
    pdf.savefig(fig)
    plt.close(fig)


def plot_pt(df: pd.DataFrame, pdf: PdfPages):
    # linear x-axis
    fig, ax = plt.subplots()
    ax.set_xlabel("pT [GeV]")
    ax.set_ylabel("Counts")
    bins = np.linspace(0, 1, 500)
    ax.hist(df["pt"], bins=bins)
    ax.semilogy()
    pdf.savefig(fig)
    plt.close(fig)

    # logarithmic x-axis
    fig, ax = plt.subplots()
    ax.set_xlabel("pT [GeV]")
    ax.set_ylabel("Counts")
    bins = np.logspace(-6, 0, 500)
    ax.hist(df["pt"], bins=bins)
    ax.semilogx()
    ax.semilogy()
    pdf.savefig(fig)
    plt.close(fig)

def plot_t(df: pd.DataFrame, pdf: PdfPages):
    fig, ax = plt.subplots()
    ax.set_xlabel("t [ns]")
    ax.set_ylabel("Counts")
    bins = np.linspace(-50, 100, 150)
    ax.hist(df["t"], bins=bins)
    ax.semilogy()
    pdf.savefig(fig)
    plt.close(fig)

if __name__ == "__main__":
    main()
