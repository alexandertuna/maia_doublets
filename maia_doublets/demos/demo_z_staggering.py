"""
All units are mm!
Layers are treated as infinitely thin.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams

ZPOS = np.array([
    [0.0, 60.2, 120.4, 180.6, 240.8, 301.0, 361.2, 421.4, 481.6, 541.8, 602.0, 662.2, 722.4, 782.6, 842.8, 903.0, 963.2, 1023.4, 1083.6, 1143.8, 1204.0],
    [0.0, 60.2, 120.4, 180.6, 240.8, 301.0, 361.2, 421.4, 481.6, 541.8, 602.0, 662.2, 722.4, 782.6, 842.8, 903.0, 963.2, 1023.4, 1083.6, 1143.8, 1204.0],
])

NZS = np.array([
    42 // 2,
    42 // 2,
    # 42 // 2,
    # 42 // 2,
    # 42 // 2,
    # 42 // 2,
    # 42 // 2,
    # 42 // 2,
])
PHI_STAGGER_DR = 4.0
GAP = 2.0 # mm
LAYER_RADII_Z_MOD_2_EQ_0 = {
    "v01": np.array([
        819,
        819 + GAP,
    ])
}
LAYER_RADII_Z_MOD_2_EQ_1 = {
    "v01": np.array([
        819,
        819 + GAP,
    ])
}
LENGTH_IT = 30.1 # mm
LENGTH_OT = 60.2 # mm
LENGTHS = np.array([
    LENGTH_OT,
    LENGTH_OT,
    # LENGTH_OT,
    # LENGTH_OT,
    # LENGTH_OT,
    # LENGTH_OT,
    # LENGTH_OT,
    # LENGTH_OT,
])
LAYER_PAIRS = [
    (0, 1),
    # (2, 3),
    # (4, 5),
    # (6, 7),
]


def main():
    version = "v01"
    phi_mod_2 = 0
    zstag = ZStaggering(version=version,
                        phi_mod_2=phi_mod_2,
                        )
    zstag.optimize()
    zstag.evaluate()
    zstag.plot()


class Module:
    def __init__(self,
                 layer: int,
                 rad: float,
                 z_min: float,
                 z_max: float,
                 z_module: int,
                 ):
        self.layer = layer
        self.rad = rad
        self.z_min = z_min
        self.z_max = z_max
        self.z_module = z_module
        self.color = "blue" if (z_module % 2 == 0) else "red"
        self.thetamin = np.arctan2(self.rad, self.z_max)
        self.thetamax = np.arctan2(self.rad, self.z_min)
        self.etamin = -np.log(np.tan(self.thetamax / 2.0))
        self.etamax = -np.log(np.tan(self.thetamin / 2.0))

    def contains(self, eta: float) -> bool:
        return (eta >= self.etamin) & (eta <= self.etamax)


class ZStaggering:

    def __init__(self,
                 version: str,
                 phi_mod_2: int,
                 ):
        self.version = version
        self.phi_mod_2 = phi_mod_2

        # i/o
        self.pdf = "z_stagger.pdf"

        # scanning parameters
        angles = np.deg2rad(np.linspace(90.0, 55.0, int(1e4))) # from eta=0 to eta=0.6ish
        self.etas = -np.log(np.tan(angles / 2.0))

        # detector parameters
        self.radii_z_mod_2_eq_0 = LAYER_RADII_Z_MOD_2_EQ_0[self.version] + (PHI_STAGGER_DR * self.phi_mod_2)
        self.radii_z_mod_2_eq_1 = LAYER_RADII_Z_MOD_2_EQ_1[self.version] + (PHI_STAGGER_DR * self.phi_mod_2)
        self.n_layers = len(self.radii_z_mod_2_eq_0)
        self.z_mins = ZPOS
        self.z_maxs = ZPOS + LENGTHS[:, np.newaxis]

        # build the detector
        self.modules = []
        for layer in range(self.n_layers):
            self.modules.append([])
            for z_module, (z_min, z_max) in enumerate(zip(self.z_mins[layer],
                                                        self.z_maxs[layer])):
                if z_module % 2 == 0:
                    rad = self.radii_z_mod_2_eq_0[layer]
                else:
                    rad = self.radii_z_mod_2_eq_1[layer]
                self.modules[-1].append(Module(layer=layer,
                                               rad=rad,
                                               z_min=z_min,
                                               z_max=z_max,
                                               z_module=z_module
                                               ))


    def optimize(self):
        pass


    def evaluate(self):
        # evaluate individual layers
        passes_thru = []
        for layer in range(self.n_layers):
            modules = self.modules[layer]
            passes_thru.append(np.array([module.contains(self.etas) for module in modules]))

        # evaluate layer pairs
        self.passes_thru_both = {}
        for (lower, upper) in LAYER_PAIRS:
            lower_pass = passes_thru[lower]
            upper_pass = passes_thru[upper]
            if len(lower_pass) != len(upper_pass):
                raise ValueError("There are different n(z-sensors) for lower and upper")

            passes_thru_both = np.zeros_like(self.etas).astype(bool)
            for lower_module_mask, upper_module_mask in zip(lower_pass, upper_pass):
                passes_thru_both |= (lower_module_mask & upper_module_mask)

            passes_thru_both = passes_thru_both.astype(float)
            self.passes_thru_both[(lower, upper)] = passes_thru_both


    def plot(self):
        with PdfPages(self.pdf) as pdf:
            self.draw_detector(pdf)
            self.plot_passing(pdf)


    def draw_detector(self, pdf: PdfPages):
        fig, ax = plt.subplots()
        for layer in range(self.n_layers):
            for module in self.modules[layer]:
                ax.plot([module.z_min, module.z_max], [module.rad, module.rad], color=module.color, linewidth=2)
        ax.set_xlabel("Z [mm]")
        ax.set_ylabel("R [mm]")
        ax.set_title(f"Detector layout. phi % 2 = {self.phi_mod_2}")
        pdf.savefig()
        plt.close()


    def plot_passing(self, pdf: PdfPages):
        for (lower, upper), passes_thru_both in self.passes_thru_both.items():
            percent_ok = np.mean(passes_thru_both)
            fig, ax = plt.subplots()
            ax.plot(self.etas, passes_thru_both, marker="o", markersize=1, linewidth=2, linestyle="-", color="blue")
            ax.set_xlabel("Eta")
            ax.set_ylabel("Passes through double-layer")
            ax.set_title(f"Layer pair {lower}, {upper}. phi % 2 = {self.phi_mod_2}. Passing: {percent_ok:.2%}")
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
