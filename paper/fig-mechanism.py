"""
Equal spaced detector aka v07
"""
import style
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

MIN_PT = 0.2 # GeV
MIN_T = -10 # ns
MAX_T = 20 # ns
MIN_Z = -4000 # mm
MAX_Z = 4000 # mm
MIN_R = 0
MAX_R = 1900 # mm
MCPARTICLE = "MCParticle"
N_BIB_FILES = 5
BIB = [
    f"/ceph/users/atuna/work/maia/maia_datasets/productions/bib.2026_08_14_17h50m00s/BIB10TeV/sim_{muon}/BIB_sim_{i+1}.slcio"
    for i in range(N_BIB_FILES)
    for muon in ["mp", "mm"]
]


PKL = "fig-mechanism.pkl"

def main():

    if os.path.exists(PKL):
        df = pd.read_pickle(PKL)
    else:
        df = get_mcparticles(BIB)
        df = post_process(df)
        df.to_pickle(PKL)

    print(df)
    with PdfPages("fig-mechanism.pdf") as pdf:
        plot(df, pdf)


def get_mcparticles(bib) -> pd.DataFrame:
    """
    https://github.com/MuonColliderSoft/LCIO/blob/master/src/cpp/include/IMPL/MCParticleImpl.h
    """
    import pyLCIO
    all_data = []
    for i_fpath, fpath in enumerate(bib):
        print(f"Processing file {i_fpath + 1}/{len(bib)}")
        reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
        reader.open(fpath)
        for event in reader:
            for mcparticle in event.getCollection(MCPARTICLE):
                all_data.append({
                    "pdg": mcparticle.getPDG(),
                    "t": mcparticle.getTime(),
                    "x": mcparticle.getVertex()[0],
                    "y": mcparticle.getVertex()[1],
                    "z": mcparticle.getVertex()[2],
                    "px": mcparticle.getMomentum()[0],
                    "py": mcparticle.getMomentum()[1],
                    "pz": mcparticle.getMomentum()[2],
                    "gen_status": mcparticle.getGeneratorStatus(),
                    "sim_status": mcparticle.getSimulatorStatus(),
                    "nparents": len(mcparticle.getParents()),
                    "nchildren": len(mcparticle.getDaughters()),
                })
        reader.close()
    return pd.DataFrame(all_data)


def post_process(df: pd.DataFrame) -> pd.DataFrame:
    df["r"] = (df["x"]**2 + df["y"]**2)**0.5
    df["pt"] = (df["px"]**2 + df["py"]**2)**0.5
    df["p"] = (df["px"]**2 + df["py"]**2 + df["pz"]**2)**0.5
    df = df[df["pt"] > MIN_PT]
    df = df[(df["t"] > MIN_T) & (df["t"] < MAX_T)]
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


    # ax.set_title("Mechanism")
    # ax.set_xlabel("z [mm]")
    # ax.set_ylabel("r [mm]")
    # ax.scatter(df["z"], df["r"], c=df["rz"], cmap="viridis", s=1)
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
