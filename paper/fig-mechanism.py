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

MIN_PT = 1.0 # GeV
MIN_T = -10 # ns
MAX_T = 20 # ns
MIN_Z = -4000 # mm
MAX_Z = 4000 # mm
MIN_R = 0
MAX_R = 1900 # mm
MCPARTICLE = "MCParticle"
GEN_STATUS = 1
N_BIB_FILES = 20
BIB = [
    f"/ceph/users/atuna/work/maia/maia_datasets/productions/bib.2026_08_14_17h50m00s/BIB10TeV/sim_{muon}/BIB_sim_{i+1}.slcio"
    for i in range(N_BIB_FILES)
    for muon in ["mp", "mm"]
]
TTBAR = [
    "/ceph/users/atuna/work/maia/maia_noodling/experiments/simulate_ttbar.2026_07_08_10h14m00s/ttbar_sim/ttbar_sim_10000.slcio", # event 1 is quite central
    # "/ceph/users/atuna/work/maia/maia_datasets/productions/ttbar.2026_09_19_13h43m00s/ttbar_sim_1.slcio",
]
MAX_PROC = 20
B_FIELD = 5 # T
T_STEP = 3 # ns
T_STEPS = 3


PKL = "fig-mechanism.pkl"

def main():

    # check_total_z_momentum(TTBAR)
    # return

    if os.path.exists(PKL):
        df = pd.read_pickle(PKL)
    else:
        df = get_mcparticles(BIB + TTBAR)
        df.to_pickle(PKL)

    print(df)
    print("Number of particles:", len(df))
    print("Number of BIB particles:", df["is_bib"].sum())
    print("Number of ttbar particles:", (~df["is_bib"]).sum())

    with PdfPages("fig-mechanism.pdf") as pdf:
        plot(df, pdf)


def check_total_z_momentum(fpaths: list[str]):
    info = []
    for fpath in fpaths:
        for event in range(10):
            df = get_mcparticles_one_file(fpath, events_of_interest=[event])
            num, mean, std = len(df), df["pz"].mean(), df["pz"].std()
            info.append([fpath, event, mean, std])

    for fpath, event, mean, std in info:
        print(f"File {os.path.basename(fpath)} event {event} num pz: {num:6d}, mean pz: {mean:6.1f}, std pz: {std:6.1f}")


def is_bib_file(fpath: str) -> bool:
    bname = os.path.basename(fpath)
    if "BIB" in bname:
        return True
    elif "ttbar" in bname:
        return False
    raise ValueError(f"Unknown file type: {fpath}")


def get_mcparticles(fpaths: list[str]) -> pd.DataFrame:
    dfs = []
    with mp.Pool(processes=MAX_PROC) as pool:
        dfs = pool.map(get_mcparticles_one_file, fpaths)
    return pd.concat(dfs, ignore_index=True)


def get_mcparticles_one_file(fpath: str, events_of_interest: list[int] = [0]) -> pd.DataFrame:
    """
    https://github.com/MuonColliderSoft/LCIO/blob/master/src/cpp/include/IMPL/MCParticleImpl.h
    """
    import pyLCIO
    is_bib = is_bib_file(fpath)

    # \HACK
    if not is_bib:
        events_of_interest = [1]
    # /HACK

    print(f"Processing file {fpath}, is_bib={is_bib}")
    rows = []
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(fpath)
    for i_event, event in enumerate(reader):
        if i_event not in events_of_interest:
            print(f"Skipping event {i_event} of {os.path.basename(fpath)}")
            if i_event > max(events_of_interest):
                break
            continue
        for mcp in event.getCollection(MCPARTICLE):
            rows.append({
                "pdg": mcp.getPDG(),
                "q": mcp.getCharge(),
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
                "is_bib": is_bib,
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

    sig_mask = df["is_bib"] == False
    bkg_mask = df["is_bib"] == True

    sig_trajs = []
    for mcp in df[sig_mask].itertuples():
        step_x, step_y, step_z, step_r = mcp.x, mcp.y, mcp.z, mcp.r
        if mcp.is_bib:
            continue
        traj = []
        for _ in range(T_STEPS):
            step_x, step_y, step_z = step_x + px * T_STEP, step_y + py * T_STEP, step_z + pz * T_STEP
            step_r = (step_x**2 + step_y**2)**0.5
            traj.append((step_z, step_r))
        sig_trajs.append(traj)

    fig, ax = plt.subplots()
    ax.quiver(z[bkg_mask], r[bkg_mask], u[bkg_mask], v[bkg_mask],
              angles="xy", scale_units="xy", scale=1,  # arrow length = data units
              width=0.003)

    print("Drawing BIB ...")
    ax.scatter(z[bkg_mask], r[bkg_mask], s=2, color="black", zorder=3)
    ax.scatter(z[sig_mask], r[sig_mask], s=2, color="red", zorder=3)
    print("Drawing ttbar ...")
    for i_traj, traj in enumerate(sig_trajs):
        traj = np.array(traj)
        print(f"Trajectory {i_traj}:")
        ax.plot(traj[:, 0], traj[:, 1], color="red", linewidth=0.5, zorder=2)
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
