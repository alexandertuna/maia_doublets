"""
All units are mm!
Layers are treated as infinitely thin.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams

NZS = np.array([
    42 // 2,
    42 // 2,
    42 // 2,
    42 // 2,
    42 // 2,
    42 // 2,
    42 // 2,
    42 // 2,
])
GAP = 2.0 # mm
LAYER_RADII = {
    "v01": np.array([
        819,
        819 + GAP,
        899,
        899 + GAP,
        1366,
        1366 + GAP,
        1446,
        1446 + GAP,
    ])
}
LENGTH_IT = 30.1 # mm
LENGTH_OT = 60.2 # mm
LENGTHS = np.array([
    LENGTH_OT,
    LENGTH_OT,
    LENGTH_OT,
    LENGTH_OT,
    LENGTH_OT,
    LENGTH_OT,
    LENGTH_OT,
    LENGTH_OT,
])
LAYER_PAIRS = [
    (0, 1),
    (2, 3),
    (4, 5),
    (6, 7),
]


def main():
    zstag = ZStaggering()
    zstag.optimize()
    zstag.evaluate()
    zstag.plot()


class ZStaggering:

    def __init__(self):
        # i/o
        self.pdf = "z_stagger.pdf"

        # building the detector
        self.radii = LAYER_RADII["v01"]
        self.n_layers = len(self.radii)
        self.angles = np.deg2rad(np.linspace(90.0, 30.0, int(1e4))) # from eta=0 to eta=1.3ish
        self.etas = -np.log(np.tan(self.angles / 2.0))
        self.module_starts = np.array([[zsensor*LENGTHS[layer] for zsensor in range(NZS[layer])] for layer in range(self.n_layers)])
        self.module_stops = np.array([[(zsensor+1)*LENGTHS[layer] for zsensor in range(NZS[layer])] for layer in range(self.n_layers)])


    def optimize(self):
        pass


    def evaluate(self):
        self.passes_thru = []
        for radius, starts, stops in zip(self.radii,
                                         self.module_starts,
                                         self.module_stops,
                                         ):
            passes = []
            for (start, stop) in zip(starts, stops):
                # tan(angle) = r/z
                z = radius / np.tan(self.angles)
                passes.append( (z > start) & (z < stop) )
            self.passes_thru.append(passes)


    def plot(self):
        with PdfPages(self.pdf) as pdf:
            self.plot_passing(pdf)


    def plot_passing(self, pdf: PdfPages):

        for (lower, upper) in LAYER_PAIRS:
            lower_pass = self.passes_thru[lower]
            upper_pass = self.passes_thru[upper]
            if len(lower_pass) != len(upper_pass):
                raise ValueError("There are different n(z-sensors) for lower and upper")

            passes_thru_both = np.zeros_like(self.angles).astype(bool)
            for module, (lower_module_mask, upper_module_mask) in enumerate(zip(lower_pass, upper_pass)):
                passes_thru_both |= (lower_module_mask & upper_module_mask)

            passes_thru_both = passes_thru_both.astype(float)

            fig, ax = plt.subplots()
            ax.plot(self.etas, passes_thru_both, marker="o", markersize=1, linewidth=3, linestyle="-", color="blue")
            ax.set_xlabel("Eta")
            ax.set_ylabel("Passes through double-layer")
            ax.set_title(f"Layer pair {lower}, {upper}")
            ax.set_ylim(0.0, 1.03)
            pdf.savefig()
            plt.close()



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


if __name__ == "__main__":
    main()
