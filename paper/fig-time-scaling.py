"""
Plot of the C++ DoubletTracking time per event as a function of BIB percentage.

Each BIB file holds one event. At each percentage, the time is the average over
files 0, 1, ..., NUM_FILES-1. The times come from the timing table that Gaudi's
AlgTimingAuditor prints at the end of each log written by cpp/run_scan.sh.
"""
import style

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

RESOLUTION = "10um"
DETECTOR = "v06"
OUTPUT_DIR = "/ceph/users/atuna/work/maia/maia_doublets/output"
NEUTRINOGUNS = {
    10: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun10_digi_{RESOLUTION}",
    20: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun20_digi_{RESOLUTION}",
    30: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun30_digi_{RESOLUTION}",
    40: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun40_digi_{RESOLUTION}",
    50: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun50_digi_{RESOLUTION}",
    60: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun60_digi_{RESOLUTION}",
    70: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun70_digi_{RESOLUTION}",
    80: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun80_digi_{RESOLUTION}",
    90: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun90_digi_{RESOLUTION}",
    100: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun100_digi_{RESOLUTION}",
    110: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun110_digi_{RESOLUTION}",
    120: f"{OUTPUT_DIR}/{DETECTOR}_neutrinoGun120_digi_{RESOLUTION}",
}
PERCENTAGES = sorted(NEUTRINOGUNS.keys())
TIMEWINDOW = "none"
FILENAME_TEMPLATE = f"doublets_{TIMEWINDOW}_XXX.log"
NUM_FILES = 10
ALGORITHM = "DoubletTracking"
OUTPUT_FILE = f"fig-time-scaling-{DETECTOR}-{TIMEWINDOW}.pdf"


def main():
    times = load_times()
    with PdfPages(OUTPUT_FILE) as pdf:
        plot_times(times, pdf)


def load_times() -> dict[int, float]:
    """Average DoubletTracking time per event (s), for each BIB percentage."""
    times = {}
    for percentage in PERCENTAGES:
        times[percentage] = load_average_time(NEUTRINOGUNS[percentage])
    return times


def load_average_time(dirpath: str) -> float:
    """Average DoubletTracking time per event (s) over files 0..NUM_FILES-1 in one directory."""
    times = []
    for file_number in range(NUM_FILES):
        filename = FILENAME_TEMPLATE.replace("XXX", str(file_number))
        times.append(read_time(f"{dirpath}/{filename}"))
    return np.mean(times)


def read_time(logpath: str) -> float:
    """
    DoubletTracking time per event (s) from one log, i.e. total (s) / count from the row:
    AlgTimingAuditor     INFO  DoubletTracking               | 1.778e+05 |         1 |     177.8
    Raises an error naming the file if the row is missing (job failed or still running).
    """
    print(f"Reading {logpath} ...")
    with open(logpath) as log:
        for line in log:
            if line.startswith("AlgTimingAuditor") and f" {ALGORITHM} " in line:
                name, exec_ms, count, total_s = line.split("|")
                return float(total_s) / int(count)
    raise ValueError(f"No {ALGORITHM} timing row in {logpath} (job failed or still running?)")


def plot_times(times: dict[int, float], pdf: PdfPages) -> None:
    """Time vs BIB percentage with a power-law fit, as in fig-scaling.py."""
    percentages = np.array(PERCENTAGES)
    seconds = np.array([times[percentage] for percentage in PERCENTAGES])

    fig, ax = plt.subplots()
    ax.plot(percentages, seconds, marker="o", linestyle="None", color="black")

    # fit to power law: a straight line in log(time) vs log(percentage)
    slope, intercept = np.polyfit(np.log(percentages), np.log(seconds), 1)
    fit = np.exp(intercept) * percentages ** slope
    ax.plot(percentages, fit, linestyle="--", color="black")
    ax.text(0.05, 0.92, f"$t \\propto f^{{{slope:.1f}}}$", transform=ax.transAxes)

    padding_lo, padding_hi = 10, 9
    ax.set_xlim([min(percentages) - padding_lo, max(percentages) + padding_hi])
    ax.set_xlabel("Percentage of BIB considered")
    ax.set_ylabel("Time per event [s]")

    # lin x, lin y
    pdf.savefig(fig)
    plt.close(fig)


if __name__ == "__main__":
    main()
