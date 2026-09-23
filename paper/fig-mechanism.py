"""
Equal spaced detector aka v07
"""
import style
import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from concurrent.futures import ProcessPoolExecutor

MIN_PT = 0.5 # GeV
MIN_T = -25 # ns
MAX_T = 25 # ns
PROPAGATE_A_BIT = 0.3 # ns
PROPAGATE_A_LOT = 3.0 # ns
N_STEPS = 10
MIN_Z = -2.1 # m
MAX_Z = 2.1 # m
MAX_R = 0.9 # m
MM_TO_M = 1e-3
MCPARTICLE = "MCParticle"
GEN_STATUS = 1
N_PARENTS = 0
N_BIB_FILES_FULL = 1666
# N_BIB_FILES = 1666
# N_BIB_FILES = 833
# N_BIB_FILES = 304
N_BIB_FILES = 6665
BIB = [
    f"/ceph/users/atuna/work/maia/maia_datasets/productions/bib.2026_08_14_17h50m00s/BIB10TeV/sim_{muon}/BIB_sim_{i+1}.slcio"
    for i in range(N_BIB_FILES)
    for muon in ["mp", "mm"]
]
TTBAR = [
    "/ceph/users/atuna/work/maia/maia_datasets/productions/ttbar.2026_09_19_13h43m00s/ttbar_sim_10000.slcio",
]
PARALLEL = True
MAX_WORKERS = 10
N_FILES_PER_WORKER = 100
B_FIELD = 5 # T
SPEED_OF_LIGHT = 299.792458 # mm/ns
K_CONSTANT = 0.299792458e-3  # GeV / (T * mm) per unit charge
T_STEP = 3 # ns
T_STEPS = 3
SIGNAL_EVENT_OF_INTEREST = 5
BIB_EVENT_OF_INTEREST = 0
DOWNSAMPLE_BIB = 0.25
DPI = 1000

PKL = "fig-mechanism.pkl"

def main():

    # check_total_z_momentum(TTBAR)
    # return

    df = get_or_load_mcparticles()
    df = downsample_bib(df)
    df = remove_off_screen_mcparticles(df)

    print(df)
    print("Number of particles:", len(df))
    print("Number of BIB particles:", df["is_bib"].sum())
    print("Number of ttbar particles:", (~df["is_bib"]).sum())

    with PdfPages("fig-mechanism.pdf") as pdf:
        plot(df, pdf)


def get_or_load_mcparticles() -> pd.DataFrame:
    if os.path.exists(PKL):
        print(f"Loading {PKL} ...")
        df = pd.read_pickle(PKL)
    else:
        df = get_mcparticles(BIB + TTBAR)
        print(f"Writing to {PKL} ...")
        df.to_pickle(PKL)
    return df


def is_bib_file(fpath: str) -> bool:
    bname = os.path.basename(fpath)
    if "BIB" in bname:
        return True
    elif "ttbar" in bname:
        return False
    raise ValueError(f"Unknown file type: {fpath}")


def downsample_bib(df: pd.DataFrame) -> pd.DataFrame:
    if DOWNSAMPLE_BIB < 1.0:
        n_files = len(df[df["is_bib"]][["bib_num"]].drop_duplicates())
        print(f"Downsampling {n_files} BIB files by {DOWNSAMPLE_BIB} ...")
        df = pd.concat([
            df[df["is_bib"]].sample(frac=DOWNSAMPLE_BIB, random_state=42),
            df[~df["is_bib"]],
        ], ignore_index=True)
    return df


def remove_off_screen_mcparticles(df: pd.DataFrame) -> pd.DataFrame:
    min_z = MIN_Z / MM_TO_M
    max_z = MAX_Z / MM_TO_M
    max_r = MAX_R / MM_TO_M
    df = df[df["r"] < max_r]
    df = df[(df["z"] > min_z) & (df["z"] < max_z)]
    df = df[(df["t"] > MIN_T) & (df["t"] < MAX_T)]
    return df


def get_mcparticles(fpaths: list[str]) -> pd.DataFrame:
    chunks = [fpaths[i:i+N_FILES_PER_WORKER] for i in range(0, len(fpaths), N_FILES_PER_WORKER)]

    if PARALLEL:
        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
            futures = [pool.submit(get_mcparticles_worker, chunk, index) for index, chunk in enumerate(chunks)]
            for future in futures:
                _ = future.result()
    else:
        for index, chunk in enumerate(chunks):
            _ = get_mcparticles_worker(chunk, index)

    print(f"Reading pickled dataframes from {PKL}.* ...")
    dfs = [pd.read_pickle(f) for f in glob.glob(f"{PKL}.*")]

    print(f"Concatenating {len(dfs)} dfs ...")
    return pd.concat(dfs, ignore_index=True)


def get_mcparticles_worker(fpaths: list[str], index: int = None, event_of_interest: int = None) -> pd.DataFrame:
    """
    https://github.com/MuonColliderSoft/LCIO/blob/master/src/cpp/include/IMPL/MCParticleImpl.h
    """
    import pyLCIO

    if len(fpaths) == 0:
        raise ValueError("No file paths provided to get_mcparticles_worker.")
    
    rows = []
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

    for fpath in fpaths:

        if not os.path.exists(fpath):
            print(f"File {fpath} does not exist, skipping.")
            continue

        # filename parsing
        is_bib = is_bib_file(fpath)
        is_mm = "sim_mm" in fpath
        bib_num = int(re.search(r"BIB_sim_(\d+)\.slcio", fpath).group(1)) if is_bib else 0

        # choosing an event
        if event_of_interest is None:
            eoi = SIGNAL_EVENT_OF_INTEREST if not is_bib else BIB_EVENT_OF_INTEREST
        else:
            eoi = event_of_interest
        print(f"Processing {fpath}, is_bib={is_bib}, event_of_interest={eoi}, bib_num={bib_num}")

        counter = 0
        reader.open(fpath)

        # EVENT::LCEvent* evt = lcReader->readEvent(targetRun, targetEvent);
        reader.skipNEvents(eoi)
        for i_event, event in enumerate(reader):
            if i_event > 0:
                break
            mcps = event.getCollection(MCPARTICLE)
            if not mcps:
                raise ValueError(f"No MCParticles found in {fpath}, is_bib={is_bib}, event_of_interest={eoi}, bib_num={bib_num}")
            for mcp in mcps:
                # save memory
                if not mcp:
                    raise ValueError(f"Invalid MCParticle in {fpath}, is_bib={is_bib}, event_of_interest={eoi}, bib_num={bib_num}")
                if mcp.getGeneratorStatus() != GEN_STATUS:
                    continue
                if (mcp.getMomentum()[0]**2 + mcp.getMomentum()[1]**2) < MIN_PT**2:
                    continue
                rows.append({
                    "pdg": mcp.getPDG(),
                    "m": mcp.getMass(),
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
                    "is_mm": is_mm,
                    "bib_num": bib_num,
                })
                counter += 1

        print(f"Processing {fpath}, is_bib={is_bib}, event_of_interest={eoi}, bib_num={bib_num}, counter={counter}")
        reader.close()

    # post-processing
    df = post_process(pd.DataFrame(rows))

    # For big data processing,
    # save data to disk instead of holding it in memory
    if index is not None:
        print(f"Saving DataFrame to {PKL}.{index}")
        df.to_pickle(f"{PKL}.{index:04}")
        return

    return df


def post_process(df: pd.DataFrame) -> pd.DataFrame:
    df["r"] = (df["x"]**2 + df["y"]**2)**0.5
    df["pt"] = (df["px"]**2 + df["py"]**2)**0.5
    df["p"] = (df["px"]**2 + df["py"]**2 + df["pz"]**2)**0.5
    df["beta"] = df["p"] / np.sqrt(df["p"]**2 + df["m"]**2)
    # downscope column datatypes
    for col in ("m", "t", "x", "y", "z", "r", "px", "py", "pz", "pt", "p", "beta"):
        df[col] = df[col].astype(np.float32)
    for col in ("pdg", "nparents", "nchildren", "bib_num"):
        df[col] = df[col].astype(np.int32)
    return df


def plot(df: pd.DataFrame, pdf: PdfPages):
    # plot_pt(df, pdf)
    # plot_t(df, pdf)
    # plot_rz(df, pdf)
    # plot_xy(df, pdf)
    plot_xyz(df=df, t_max=PROPAGATE_A_BIT, pdf=pdf)
    plot_xyz(df=df, t_max=PROPAGATE_A_LOT, pdf=pdf)


def plot_rz(df: pd.DataFrame, pdf: PdfPages):

    print("Starting plot_rz...")

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

    # sig_trajs = []
    # for mcp in df[sig_mask].itertuples():
    #     step_x, step_y, step_z, step_r = mcp.x, mcp.y, mcp.z, mcp.r
    #     if mcp.is_bib:
    #         continue
    #     traj = []
    #     for _ in range(T_STEPS):
    #         step_x, step_y, step_z = step_x + px * T_STEP, step_y + py * T_STEP, step_z + pz * T_STEP
    #         step_r = (step_x**2 + step_y**2)**0.5
    #         traj.append((step_z, step_r))
    #     sig_trajs.append(traj)

    fig, ax = plt.subplots()
    # ax.quiver(z[bkg_mask], r[bkg_mask], u[bkg_mask], v[bkg_mask],
    #           angles="xy", scale_units="xy", scale=1,  # arrow length = data units
    #           width=0.003)

    print("Drawing BIB ...")
    ax.scatter(z[bkg_mask], r[bkg_mask], s=2, color="black", zorder=3)
    ax.scatter(z[sig_mask], r[sig_mask], s=2, color="red", zorder=3)
    print("Drawing ttbar ...")
    # for i_traj, traj in enumerate(sig_trajs):
    #     traj = np.array(traj)
    #     print(f"Trajectory {i_traj}:")
    #     ax.plot(traj[:, 0], traj[:, 1], color="red", linewidth=0.5, zorder=2)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_xlim(MIN_Z, MAX_Z)
    ax.set_ylim(0, MAX_R)
    # ax.set_aspect("equal")  # otherwise arrow angles look distorted
    pdf.savefig(fig)
    plt.close(fig)


def plot_xy(df: pd.DataFrame, pdf: PdfPages):
    pass


def plot_xyz(df: pd.DataFrame, t_max: float, pdf: PdfPages):

    is_signal = df["is_bib"] == False

    print("Propagating helices ...")
    x, y, z = (df[col].to_numpy() for col in ("x", "y", "z"))
    px, py, pz = (df[col].to_numpy() for col in ("px", "py", "pz"))
    q, beta = (df[col].to_numpy() for col in ("q", "beta"))
    trace_x, trace_y, trace_z = propagate_helix(x, y, z, px, py, pz, q, beta, t_max=t_max)
    trace_x *= MM_TO_M
    trace_y *= MM_TO_M
    trace_z *= MM_TO_M

    segs = np.stack([trace_z, trace_x, trace_y], axis=-1)  # (N, n_steps, 3); beam axis horizontal

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(projection="3d")
    ax.computed_zorder = False  # respect zorder so signal draws on top

    for mask, color, lw, alpha, zo, label in [
            ( is_signal, "red", 0.2, 0.1, 1, "Signal"),
            (~is_signal, "blue", 0.2, 0.1, 2, "Background"),
        ]:
        print(f"Drawing lines for {label} ...")
        lc = Line3DCollection(segs[mask],
                              colors=color,
                              linewidths=lw,
                              alpha=alpha,
                              zorder=zo,
                              label=label,
                              )
        lc.set_clip_on(False)
        lc.set_rasterized(True)  # small PDF, vector text
        ax.add_collection3d(lc)

    print("Drawing 3D figure ...")
    ax.set(xlabel="z [m]", ylabel="x [m]", zlabel="y [m]")
    physical_xlim, physical_ylim, physical_zlim = (-MAX_R, MAX_R), (-MAX_R, MAX_R), (MIN_Z, MAX_Z)
    ax.set_xlim(*physical_zlim)
    ax.set_ylim(*physical_xlim)
    ax.set_zlim(*physical_ylim)
    ax.set_box_aspect((np.ptp(physical_zlim), np.ptp(physical_xlim), np.ptp(physical_ylim)), zoom=1.2)

    # draw beamline as black line from z = MIN_Z to z = MAX_Z
    ax.plot([MIN_Z, MAX_Z], [0, 0], [0, 0], color="black", linewidth=0.5, alpha=0.5, zorder=0, clip_on=False)


    # add text: legend
    kwargs = dict(rotation=-4, fontsize=12, transform=ax.transAxes)
    ax.text2D(0.40, 0.22 if t_max > 1 else 0.28, f"Collision", **kwargs, color="red")
    ax.text2D(0.08, 0.28 if t_max > 1 else 0.35, f"Beam background", **kwargs, color="blue")

    # add text: minimum PT threshold, propagation time, bib scaling
    # angle text by 10 degrees
    rotation = -4
    fontsize = 10
    kwargs = dict(rotation=rotation, fontsize=fontsize, transform=ax.transAxes)
    bib_scaling = N_BIB_FILES / N_BIB_FILES_FULL * DOWNSAMPLE_BIB
    ax.text2D(0.08, 0.830, r"$p_{T} >$" + f"{MIN_PT} GeV", **kwargs)
    ax.text2D(0.44, 0.805, r"$\Delta t =$" + f"{t_max} ns", **kwargs)
    if bib_scaling < 1.0:
        ax.text2D(0.80, 0.780, f"BIB scaling: {bib_scaling:.1f}", **kwargs)

    # manage the grid and axes
    ax.grid(True, linewidth=0.1)
    for axis in (
        ax.xaxis,
        ax.yaxis,
        ax.zaxis,
    ):
        axis.set_tick_params(pad=-5 if axis == ax.xaxis else 0)
        axis.labelpad = -5
        axis.set_major_locator(MaxNLocator(5 if axis == ax.xaxis else 3))
        axis.label.set_size(10)
        axis.set_tick_params(labelsize=10)
        axis.pane.fill = False
        axis._axinfo['grid']['color'] = (0.9, 0.9, 0.9, 0.3)

    print("Adjusting figure layout and view ...")
    fig.subplots_adjust(left=0, right=0.95, bottom=0, top=1)
    ax.view_init(elev=25, azim=-80)
    # ax.legend(loc="upper left", frameon=False)

    print("Saving figure to PNG ...")
    tag = "early" if t_max < 1 else "later"
    fig.savefig(f"fig-mechanism-{tag}.png", dpi=DPI)

    print("Saving figure to PDF ...")
    pdf.savefig(fig, dpi=DPI)
    plt.close(fig)


def propagate_helix(x, y, z, px, py, pz, q, beta, t_max, B=B_FIELD, n_steps=N_STEPS):
    """Units: mm, GeV, ns; B along +z. Returns X, Y, Z of shape (N, n_steps)."""
    p  = np.sqrt(px**2 + py**2 + pz**2)
    pt = np.hypot(px, py)
    s = np.linspace(0, 1, n_steps)[None, :] * (beta * SPEED_OF_LIGHT * t_max)[:, None]
    phi0 = np.arctan2(py, px)[:, None]
    w = (-q * K_CONSTANT * B / p)[:, None] # signed d(phi)/ds
    ws = w * s
    small = np.abs(ws[:, -1:]) < 1e-6 # neutrals / stiff tracks -> straight line
    w_safe = np.where(small, 1.0, w)
    fx = np.where(small, s * np.cos(phi0),  (np.sin(phi0 + ws) - np.sin(phi0)) / w_safe)
    fy = np.where(small, s * np.sin(phi0), -(np.cos(phi0 + ws) - np.cos(phi0)) / w_safe)
    trace_x = x[:, None] + (pt / p)[:, None] * fx
    trace_y = y[:, None] + (pt / p)[:, None] * fy
    trace_z = z[:, None] + (pz / p)[:, None] * s
    return trace_x, trace_y, trace_z


def check_total_z_momentum(fpaths: list[str]):
    """
    This is helpful for picking a more-central event
    """
    info = []
    for fpath in fpaths:
        for event in range(10):
            df = get_mcparticles_one_file(fpath, event_of_interest=event)
            num, mean, std = len(df), df["pz"].mean(), df["pz"].std()
            info.append([fpath, event, num, mean, std])

    for fpath, event, num, mean, std in info:
        print(f"File {os.path.basename(fpath)} event {event} num: {num:6d}, mean pz: {mean:6.1f}, std pz: {std:6.1f}")


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
    for is_bib in (True, False):
        fig, ax = plt.subplots()
        ax.set_xlabel("t [ns]")
        ax.set_ylabel("Counts")
        bins = np.linspace(-25, 25, 200)
        ax.hist(df[df["is_bib"] == is_bib]["t"], bins=bins)
        ax.semilogy()
        ax.set_title("BIB" if is_bib else "Non-BIB")
        pdf.savefig(fig)
        plt.close(fig)

if __name__ == "__main__":
    main()
