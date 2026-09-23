"""
Plot of neutrinoGun object multiplicity as a function of BIB percentage
"""
import style

import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

NEUTRINOGUNS = {
    10: "/ceph/users/atuna/work/maia/maia_doublets/output/v06_neutrinoGun10_digi_10um",
    20: "/ceph/users/atuna/work/maia/maia_doublets/output/v06_neutrinoGun20_digi_10um",
    30: "/ceph/users/atuna/work/maia/maia_doublets/output/v06_neutrinoGun30_digi_10um",
    40: "/ceph/users/atuna/work/maia/maia_doublets/output/v06_neutrinoGun40_digi_10um",
    50: "/ceph/users/atuna/work/maia/maia_doublets/output/v06_neutrinoGun50_digi_10um",
    70: "/ceph/users/atuna/work/maia/maia_doublets/output/v06_neutrinoGun70_digi_10um",
    80: "/ceph/users/atuna/work/maia/maia_doublets/output/v06_neutrinoGun80_digi_10um",
    100: "/ceph/users/atuna/work/maia/maia_doublets/output/v06_neutrinoGun_digi_10um",
}
PERCENTAGES = sorted(NEUTRINOGUNS.keys())
MAX_FILES = 10
OBJECTS = [
    "hits",
    "mds",
    "t2s",
    "t4s",
    # "t8s",
]
NICKNAME = {
    "hits": "Hits",
    "mds": "MDs",
    "t2s": "T2s",
    "t4s": "T4s",
    "t8s": "T8s",
}
OKAY = {
    "hits": None,
    "mds": "md_ok",
    "t2s": "t2_ok",
    "t4s": "t4_ok",
    "t8s": "t8_ok",
}
COLOR = {
    "hits": "red",
    "mds": "blue",
    "t2s": "green",
    "t4s": "orange",
    "t8s": "purple",
}
POWER_LAW_COORDS = {
    "hits": (0.06, 0.92),
    "mds": (0.14, 0.72),
    "t2s": (0.22, 0.47),
    "t4s": (0.60, 0.16),
    "t8s": (0.50, 0.50),
}


def main():
    plot = ScalingPlot()
    with PdfPages("fig-scaling.pdf") as pdf:
        plot.plot(pdf)


class ScalingPlot:

    def __init__(self):
        self.load_data()


    def load_data(self) -> None:
        self.data = {}
        for percentage, dirpath in NEUTRINOGUNS.items():
            self.data[percentage] = self.load_data_one_percentage(dirpath)


    def load_data_one_percentage(self, dirpath) -> dict:
        return {
            obj: self.load_data_one_object(dirpath, obj)
            for obj in OBJECTS
        }


    def load_data_one_object(self, dirpath: str, obj: str) -> float:
        """
        For example, load all the mds which have `md_ok=True`.
        And find the average number of mds per event like this.
        """
        okay = OKAY[obj]
        pklpath = f"{dirpath}/{obj}_*.pkl"
        fpaths = sorted(glob.glob(pklpath))[:MAX_FILES]
        if not fpaths:
            raise FileNotFoundError(f"No files found for {pklpath}")
        total = 0
        for fpath in fpaths:
            print(f"Loading {fpath} ...")
            df = pd.read_pickle(fpath)
            total += (df[okay].sum() if (okay is not None and len(df) > 0) else len(df))
        return total / len(fpaths)


    def plot(self, pdf):
        fig, ax = plt.subplots()
        for obj in OBJECTS:
            yields = [self.data[percentage][obj] for percentage in PERCENTAGES]
            ax.plot(PERCENTAGES, yields, label=obj, marker="o", linestyle="None", color=COLOR[obj])

            # fit to power law
            text_x, text_y = POWER_LAW_COORDS[obj]
            name = NICKNAME[obj]
            coeffs = np.polyfit(np.log(PERCENTAGES), np.log(yields), 1)
            fit = np.exp(coeffs[1]) * np.array(PERCENTAGES) ** coeffs[0]
            ax.plot(PERCENTAGES, fit, linestyle="--", color=COLOR[obj])
            kwargs = dict(color=COLOR[obj], transform=ax.transAxes)
            ax.text(text_x, text_y, f"{name}: $y = x^{{{coeffs[0]:.1f}}}$", **kwargs)

        ax.set_xlabel("BIB percentage")
        ax.set_ylabel("Average multiplicity")

        # lin x, lin y
        pdf.savefig(fig)

        # lin x, log y
        ax.semilogy()
        pdf.savefig(fig)

        # log x, log y
        ax.semilogx()
        ax.semilogy()
        pdf.savefig(fig)

        # fin
        plt.close(fig)



if __name__ == "__main__":
    main()
