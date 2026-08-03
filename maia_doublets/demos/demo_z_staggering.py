"""
All units are mm!
Layers are treated as infinitely thin.
Modules are optimized sequentially because we only change the z-position of odd-numbered modules, and the even-numbered modules are fixed.
This means that the optimization of one module does not affect the optimization of another module.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams

DO_STAGGER = True
DO_OPTIMIZE = True
DO_OFFSETS = False
DO_PLOT = True
GAP = 2.0 # mm
MODULE_THICKNESS = 1.663 # mm
PHI_STAGGER_DR = 4.0 # mm

L_IT_0 = 30.1 # mm
L_OT_0 = 60.2 # mm
NZ_OT = 42 // 2
ZPOS_OT = [it*L_OT_0 for it in range(NZ_OT)]
ZPOS = np.array([
    ZPOS_OT,
    ZPOS_OT,
    ZPOS_OT,
    ZPOS_OT,
    ZPOS_OT,
    ZPOS_OT,
    ZPOS_OT,
    ZPOS_OT,
])

if DO_STAGGER:
    Z_STAGGER_DR = 4.0 # mm
    L_IT = 31.3 # mm
    L_OT = 62.6 # mm
else:
    Z_STAGGER_DR = 0.0 # mm
    L_IT = L_IT_0 # mm
    L_OT = L_OT_0 # mm

CHOSEN_OFFSETS_OT = np.array([0, 1.2, 0, 1.3, 0, 1.6, 0, 2.1, 0, 2.6, 0, 2.8, 0, 2.8, 0, 2.8, 0, 3.1, 0, 3.3, 0])
LENGTHS = np.array([
    L_OT,
    L_OT,
    L_OT,
    L_OT,
    L_OT,
    L_OT,
    L_OT,
    L_OT,
])
LAYER_PAIRS = [
    (0, 1),
    (2, 3),
    (4, 5),
    (6, 7),
]
LAYER_RADII_Z_MOD_2_EQ_0 = {
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
LAYER_RADII_Z_MOD_2_EQ_1 = {
    "v01": LAYER_RADII_Z_MOD_2_EQ_0["v01"] + Z_STAGGER_DR
}


def main():
    check_constants()
    version = "v01"
    phi_mod_2 = 0
    zstag = ZStaggering(version=version,
                        phi_mod_2=phi_mod_2,
                        )
    if DO_OPTIMIZE:
        zstag.optimize()
    else:
        zstag.evaluate()
        zstag.announce()
        if DO_PLOT:
            zstag.plot()


def check_constants():
    assert DO_STAGGER in [True, False], "DO_STAGGER must be a boolean"
    assert len(ZPOS) == len(LENGTHS), "ZPOS and LENGTHS must have the same length"
    assert len(ZPOS) == len(LAYER_PAIRS) * 2, "ZPOS must have twice the length of LAYER_PAIRS"
    assert len(ZPOS) == len(LAYER_RADII_Z_MOD_2_EQ_0["v01"]), "ZPOS and LAYER_RADII_Z_MOD_2_EQ_0 must have the same length"
    assert len(ZPOS) == len(LAYER_RADII_Z_MOD_2_EQ_1["v01"]), "ZPOS and LAYER_RADII_Z_MOD_2_EQ_1 must have the same length"
    assert len(ZPOS_OT) == NZ_OT, "ZPOS_OT must have length NZ_OT"
    assert len(CHOSEN_OFFSETS_OT) == NZ_OT, "CHOSEN_OFFSETS_OT must have length NZ_OT"


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
        self.theta_min = np.arctan2(self.rad, self.z_max)
        self.theta_max = np.arctan2(self.rad, self.z_min)
        self.eta_min = -np.log(np.tan(self.theta_max / 2.0))
        self.eta_max = -np.log(np.tan(self.theta_min / 2.0))

    def contains(self, eta: float) -> bool:
        return (eta >= self.eta_min) & (eta <= self.eta_max)


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
        if DO_OPTIMIZE:
            self.etas = np.linspace(1e-5, 0.80, int(1e5))
        else:
            self.etas = np.linspace(1e-5, 0.65, int(1e5))

        # detector parameters
        self.radii_z_mod_2_eq_0 = LAYER_RADII_Z_MOD_2_EQ_0[self.version] + (PHI_STAGGER_DR * self.phi_mod_2)
        self.radii_z_mod_2_eq_1 = LAYER_RADII_Z_MOD_2_EQ_1[self.version] + (PHI_STAGGER_DR * self.phi_mod_2)
        self.n_layers = len(self.radii_z_mod_2_eq_0)
        self.z_mins = ZPOS
        if DO_OFFSETS:
            for layer in range(self.n_layers):
                self.z_mins[layer] += np.array(CHOSEN_OFFSETS_OT)
        self.z_maxs = self.z_mins + LENGTHS[:, np.newaxis]

        # build detector
        self.build_detector()


    def build_detector(self):
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
        def choose_best(params, results):
            arr = np.array(results)
            plateau = np.flatnonzero(arr >= arr.max())
            choice = plateau[len(plateau) // 2]
            return params[choice], arr[choice]

        modules_under_test = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
        chosen_offsets = [0] * NZ_OT
        zpos_ot_offset = ZPOS_OT.copy()
        offsets = np.arange(0, 4, 0.1)

        # optimize each module sequentially
        for module in modules_under_test:

            # optimize module
            efficiencies = []
            
            for offset in offsets:
                zpos_ot_offset[module] = ZPOS_OT[module] + offset
                self.z_mins = np.array([
                    zpos_ot_offset, zpos_ot_offset,
                    zpos_ot_offset, zpos_ot_offset,
                    zpos_ot_offset, zpos_ot_offset,
                    zpos_ot_offset, zpos_ot_offset,
                ])
                self.z_maxs = self.z_mins + LENGTHS[:, np.newaxis]
                self.build_detector()
                self.evaluate()
                efficiencies.append(self.efficiency_thru_all)
                # print(f"Efficiency, offset={offset:.1f}: {self.efficiency_thru_all:.5f}")

            # save the best offset
            offset, eff = choose_best(offsets, efficiencies)
            print(f"Best offset, module {module}: {offset:.2f}, eff: {eff:.5f}, effs = {' '.join([f'{eff:.5f}' for eff in efficiencies])}")
            chosen_offsets[module] = float(offset)
            zpos_ot_offset[module] = ZPOS_OT[module] + offset

        # let the people know
        print(f"Chosen offsets: {chosen_offsets}")


    def evaluate(self):
        # evaluate individual layers
        passes_thru = []
        for layer in range(self.n_layers):
            modules = self.modules[layer]
            passes_thru.append(np.array([module.contains(self.etas) for module in modules]))

        # evaluate layer pairs
        self.passes_thru_both = {}
        self.efficiency_thru_both = {}
        for (lower, upper) in LAYER_PAIRS:
            lower_pass = passes_thru[lower]
            upper_pass = passes_thru[upper]
            if len(lower_pass) != len(upper_pass):
                raise ValueError("There are different n(z-sensors) for lower and upper")

            # does a trajectory pass thru any pair of modules?
            passes_thru_both = np.zeros_like(self.etas).astype(bool)
            for lower_module_mask, upper_module_mask in zip(lower_pass, upper_pass):
                passes_thru_both |= (lower_module_mask & upper_module_mask)

            passes_thru_both = passes_thru_both.astype(float)
            self.passes_thru_both[(lower, upper)] = passes_thru_both
            self.efficiency_thru_both[(lower, upper)] = passes_thru_both.mean()

        # does a trajectory pass thru all layer pairs?
        self.passes_thru_all = np.ones_like(self.etas).astype(bool)
        for passes_thru_both in self.passes_thru_both.values():
            self.passes_thru_all &= passes_thru_both.astype(bool)
        self.efficiency_thru_all = self.passes_thru_all.mean()


    def announce(self, verbose: bool = True):
        if verbose:
            for (lower, upper), efficiency in self.efficiency_thru_both.items():
                print(f"Efficiency thru layer pair ({lower}, {upper}): {efficiency:.5f}")
        print(f"Efficiency thru all layer pairs: {self.efficiency_thru_all:.5f}")

    def plot(self):
        with PdfPages(self.pdf) as pdf:
            self.draw_detector(pdf, eta_lines=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
            self.plot_passing(pdf)


    def draw_detector(self, pdf: PdfPages, eta_lines: list[float] = []):
        fig, ax = plt.subplots()
        for layer in range(self.n_layers):
            for module in self.modules[layer]:
                # plot entire module, including thickness in the r-direction, as a rectangle
                rect = patches.Rectangle((module.z_min, module.rad - MODULE_THICKNESS / 2),
                                          module.z_max - module.z_min,
                                          MODULE_THICKNESS,
                                          linewidth=1,
                                          edgecolor=module.color,
                                          facecolor=module.color,
                                          alpha=0.5,
                                          )
                ax.add_patch(rect)
                ax.plot([module.z_min, module.z_max], [module.rad, module.rad], color=module.color, linewidth=1)
        ax.set_xlabel("Z [mm]")
        ax.set_ylabel("R [mm]")
        z_min, z_max = ax.get_xlim()
        r_min, r_max = ax.get_ylim()
        ax.set_xlim(z_min - 10, z_max + 10)
        ax.set_ylim(r_min - 10, r_max + 10)
        for i_eta, eta in enumerate(eta_lines):
            x0 = y0 = 0
            y1 = r_max + 1 + i_eta
            x1 = y1 / np.tan(eta_to_theta(eta))
            ax.plot([x0, x1], [y0, y1], color="black", linestyle="-", alpha=0.1)
            ax.text(x1, y1, f"eta = {eta:.3f}", fontsize=10, color="lightgray", ha="right")
        ax.set_title(f"Detector layout. phi % 2 = {self.phi_mod_2}. Efficiency = {self.efficiency_thru_all:.3%}")
        pdf.savefig()
        plt.close()


    def plot_passing(self, pdf: PdfPages):
        for (lower, upper), passes_thru_both in self.passes_thru_both.items():
            percent_ok = self.efficiency_thru_both[(lower, upper)]
            print(f"Layer pair {lower}, {upper}. phi % 2 = {self.phi_mod_2}. Passing: {percent_ok:.3%}")
            fig, ax = plt.subplots()
            ax.plot(self.etas, passes_thru_both, marker="o", markersize=1, linewidth=2, linestyle="-", color="blue")
            ax.set_xlabel("Eta")
            ax.set_ylabel("Passes through double-layer")
            ax.set_title(f"Layer pair {lower}, {upper}. phi % 2 = {self.phi_mod_2}. Passing: {percent_ok:.3%}")
            ax.set_ylim(0.0, 1.03)
            pdf.savefig()
            plt.close()


def eta_to_theta(eta: float) -> float:
    return 2 * np.arctan(np.exp(-eta))


class DetectorParameters:

    def __init__(
        self,
        version: str,
        tracker: str,
        phi_mod_2: int,
        do_z_stagger: bool,
        do_z_offsets: bool,
    ):
        # check inputs
        if not tracker in ["IT", "OT"]:
            raise ValueError(f"Invalid tracker: {tracker}. Must be 'IT' or 'OT'")
        if not version in ["v01"]:
            raise ValueError(f"Invalid version: {version}. Must be 'v01'")

        self.version = version
        self.tracker = tracker
        self.phi_mod_2 = phi_mod_2
        self.do_z_stagger = do_z_stagger
        self.do_z_offsets = do_z_offsets

        # set the sensor length
        if tracker == "IT":
            raise NotImplementedError("Staggering not implemented for IT")
        else:
            self.original_length = 60.2
            self.length = 62.6 if self.do_z_stagger else 60.2

        # set the number of z-sensors to consider
        if self.tracker == "IT":
            raise NotImplementedError("nz not implemented for IT")
        else:
            self.nz = 42 // 2

        # set the dr-staggering for z-sensors
        self.z_stagger_dr = 4.0 if self.do_z_stagger else 0.0

        # set the layer radii
        if version == "v01":
            self.layer_radii_z_mod_2_eq_0 = np.array([
                819,
                819 + GAP,
                899,
                899 + GAP,
                1366,
                1366 + GAP,
                1446,
                1446 + GAP,
            ])
            self.layer_radii_z_mod_2_eq_1 = self.layer_radii_z_mod_2_eq_0 + self.z_stagger_dr
        else:
            raise ValueError(f"Invalid version: {version}. Must be 'v01'")

        # set the double-layers
        self.layer_pairs = [
            (0, 1),
            (2, 3),
            (4, 5),
            (6, 7),
        ]

        # set the z-positions of the z-sensors
        self.zpos = [it*self.original_length for it in range(self.nz)]

        # set the z-offsets of the z-sensors
        if self.do_z_offsets:
            self.chosen_offsets = np.array([0, 1.2, 0, 1.3, 0, 1.6, 0, 2.1, 0, 2.6, 0, 2.8, 0, 2.8, 0, 2.8, 0, 3.1, 0, 3.3, 0])
        else:
            self.chosen_offsets = np.array([0] * self.nz)


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
