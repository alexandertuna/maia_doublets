"""
All units are mm!
Layers are treated as infinitely thin.
Modules are optimized sequentially because we only change the z-position of odd-numbered modules, and the even-numbered modules are fixed.
This means that the optimization of one module does not (strongly?) affect the optimization of another module.

NB: Its a little complicated to choose z-offsets for each module *considering consistency between other modules*.
This matters because the optimal z-offset for each module usually has a few equal choices, i.e. a broad maximum in the efficiency vs. z-offset curve.
Right now, I just take the middle of the plateau. Then I eyeball small adjustments afterward for consistency.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
DetectorParameters = None

GAP = 2.0 # mm
MODULE_THICKNESS = 1.663 # mm
PHI_STAGGER_DR_WITH_Z_STAGGER = 8.0 # mm
PHI_STAGGER_DR_WITHOUT_Z_STAGGER = 4.0 # mm
Z_STAGGER_DR = 4.0 # mm

def main():
    ops = options()
    if ops.optimize and not ops.stagger:
        raise ValueError("Cannot optimize without staggering")
    if ops.optimize and ops.offsets:
        raise ValueError("Cannot optimize and use offsets at the same time")

    params = DetectorParameters(version=ops.version,
                                tracker=ops.tracker,
                                phi_mod_2=ops.phi_mod_2,
                                do_z_stagger=ops.stagger,
                                do_z_offsets=ops.offsets,
                                )
    params.describe()

    zstag = ZStaggering(params=params)
    if ops.optimize:
        zstag.optimize()
        return
    zstag.evaluate()
    zstag.announce()
    if ops.plot:
        zstag.plot()
    zstag.write_xml(ops.xml)


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
                 params: DetectorParameters,
                #  version: str,
                #  phi_mod_2: int,
                 ):

        # i/o
        self.pdf = "z_stagger.pdf"

        # detector parameters
        self.params = params

        # build detector
        self.build_detector()

        # scanning parameters
        self.etas = np.linspace(1e-5, 0.65, int(1e5))
        self.etas_for_optimization = np.linspace(1e-5, 0.80, int(1e5))


    def build_detector(self):
        self.modules = []
        for layer in range(self.params.n_layers):
            self.modules.append([])
            for z_module, (z_min, z_max) in enumerate(zip(self.params.z_mins[layer],
                                                          self.params.z_maxs[layer])):
                rad = (self.params.layer_radii_z_mod_2_eq_0[layer] if (z_module % 2 == 0) else
                       self.params.layer_radii_z_mod_2_eq_1[layer])
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

        # set up an un-optimized detector to start with
        if np.any(self.params.chosen_offsets != 0):
            raise ValueError("Cannot optimize if offsets are already set")

        # optimize each module sequentially
        for module in self.params.modules_under_test:

            # optimize module
            efficiencies = []
            for offset in self.params.possible_offsets:
                self.params.chosen_offsets[module] = float(offset)
                self.params.set_z_positions()
                self.build_detector()
                self.evaluate(during_optimization=True)
                efficiencies.append(self.efficiency_thru_all)

            # save the best offset
            best_offset, eff = choose_best(self.params.possible_offsets, efficiencies)
            breakdown = ""
            for (offset, eff) in zip(self.params.possible_offsets, efficiencies):
                breakdown += f"{offset:.2f}:{eff:.5f}, "
            print(f"Best offset, module {module}: {best_offset:.2f}, eff: {eff:.5f}. Breakdown: {breakdown}")
            self.params.chosen_offsets[module] = float(best_offset)

        # finalize
        print(f"Chosen offsets: {self.params.chosen_offsets}")
        self.params.set_z_positions()


    def evaluate(self, during_optimization: bool = False):

        # over-cover eta during optimization
        etas = self.etas if not during_optimization else self.etas_for_optimization

        # evaluate individual layers
        passes_thru = []
        for layer in range(self.params.n_layers):
            modules = self.modules[layer]
            passes_thru.append(np.array([module.contains(etas) for module in modules]))

        # evaluate layer pairs
        self.passes_thru_both = {}
        self.efficiency_thru_both = {}
        for (lower, upper) in self.params.layer_pairs:
            lower_pass = passes_thru[lower]
            upper_pass = passes_thru[upper]
            if len(lower_pass) != len(upper_pass):
                raise ValueError("There are different n(z-sensors) for lower and upper")

            # does a trajectory pass thru any pair of modules?
            passes_thru_both = np.zeros_like(etas).astype(bool)
            for lower_module_mask, upper_module_mask in zip(lower_pass, upper_pass):
                passes_thru_both |= (lower_module_mask & upper_module_mask)

            passes_thru_both = passes_thru_both.astype(float)
            self.passes_thru_both[(lower, upper)] = passes_thru_both
            self.efficiency_thru_both[(lower, upper)] = passes_thru_both.mean()

        # does a trajectory pass thru all layer pairs?
        self.passes_thru_all = np.ones_like(etas).astype(bool)
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
        for layer in range(self.params.n_layers):
            for module in self.modules[layer]:
                # plot entire module, including thickness in the r-direction, as a rectangle
                rect = patches.Rectangle((module.z_min, module.rad - MODULE_THICKNESS / 2),
                                          module.z_max - module.z_min,
                                          MODULE_THICKNESS,
                                          linewidth=0.1,
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
        # ax.set_xlim(0, 60)
        # ax.set_ylim(0, 200)
        # ax.set_xlim(70, 120)
        # ax.set_ylim(500, 580)
        for i_eta, eta in enumerate(eta_lines):
            x0 = y0 = 0
            y1 = r_max + 1 + i_eta
            x1 = y1 / np.tan(eta_to_theta(eta))
            ax.plot([x0, x1], [y0, y1], color="black", linestyle="-", alpha=0.1)
            ax.text(x1, y1, f"eta = {eta:.3f}", fontsize=10, color="lightgray", ha="right")
        ax.set_title(f"Detector layout. phi % 2 = {self.params.phi_mod_2}. Efficiency = {self.efficiency_thru_all:.3%}")
        pdf.savefig()
        plt.close()


    def plot_passing(self, pdf: PdfPages):
        for (lower, upper), passes_thru_both in self.passes_thru_both.items():
            percent_ok = self.efficiency_thru_both[(lower, upper)]
            print(f"Layer pair {lower}, {upper}. phi % 2 = {self.params.phi_mod_2}. Passing: {percent_ok:.3%}")
            fig, ax = plt.subplots()
            ax.plot(self.etas, passes_thru_both, marker="o", markersize=1, linewidth=2, linestyle="-", color="blue")
            ax.set_xlabel("Eta")
            ax.set_ylabel("Passes through double-layer")
            ax.set_title(f"Layer pair {lower}, {upper}. phi % 2 = {self.params.phi_mod_2}. Passing: {percent_ok:.3%}")
            ax.set_ylim(0.0, 1.03)
            pdf.savefig()
            plt.close()


    def write_xml(self, filename: str):
        self.params.make_xml_blurbs()
        with open(filename, "w") as fi:
            for blurb in self.params.xml_blurbs:
                fi.write(blurb)


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
        if not version in ["v01", "v05"]:
            raise ValueError(f"Invalid version: {version}. Must be 'v01' or 'v05'")
        if not phi_mod_2 in [0, 1]:
            raise ValueError(f"Invalid phi_mod_2: {phi_mod_2}. Must be 0 or 1")

        self.version = version
        self.tracker = tracker
        self.phi_mod_2 = phi_mod_2
        self.do_z_stagger = do_z_stagger
        self.do_z_offsets = do_z_offsets

        # set the sensor length
        if tracker == "IT":
            self.original_length = 30.1
            self.length = 32.5 if self.do_z_stagger else 30.1
        else:
            self.original_length = 60.2
            self.length = 62.6 if self.do_z_stagger else 60.2

        # set the number of z-sensors to consider
        if self.tracker == "IT":
            self.nz = 46 // 2
        else:
            self.nz = 42 // 2

        # set the dr-staggering for z-sensors
        self.z_stagger_dr = Z_STAGGER_DR if self.do_z_stagger else 0.0

        # set the layer radii
        if self.do_z_stagger:
            if version == "v01":
                if self.tracker == "IT":
                    self.layer_radii_z_mod_2_eq_0 = np.array([
                        127,
                        127 + GAP,
                        167,
                        167 + GAP,
                        502,
                        502 + GAP,
                        542,
                        542 + GAP,
                    ])
                else:
                    self.layer_radii_z_mod_2_eq_0 = np.array([
                        819,
                        819 + GAP,
                        899,
                        899 + GAP,
                        1358,
                        1358 + GAP,
                        1438,
                        1438 + GAP,
                    ])
            elif version == "v05":
                if self.tracker == "IT":
                    self.layer_radii_z_mod_2_eq_0 = np.array([
                        127,
                        127 + GAP,
                        265.333,
                        265.333 + GAP,
                        403.666,
                        403.666 + GAP,
                        542,
                        542 + GAP,
                    ])
                else:
                    self.layer_radii_z_mod_2_eq_0 = np.array([
                        819,
                        819 + GAP,
                        1025.333,
                        1025.333 + GAP,
                        1231.666,
                        1231.666 + GAP,
                        1438,
                        1438 + GAP,
                    ])
            else:
                raise ValueError(f"Invalid version: {version}. Must be 'v01' or 'v05'")
        else:
            # no z-stagger
            if version == "v01":
                if self.tracker == "IT":
                    self.layer_radii_z_mod_2_eq_0 = np.array([
                        127,
                        127 + GAP,
                        167,
                        167 + GAP,
                        510,
                        510 + GAP,
                        550,
                        550 + GAP,
                    ])
                else:
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
            elif version == "v05":
                if self.tracker == "IT":
                    self.layer_radii_z_mod_2_eq_0 = np.array([
                        127,
                        127 + GAP,
                        268,
                        268 + GAP,
                        409,
                        409 + GAP,
                        550,
                        550 + GAP,
                    ])
                else:
                    self.layer_radii_z_mod_2_eq_0 = np.array([
                        819,
                        819 + GAP,
                        1028,
                        1028 + GAP,
                        1237,
                        1237 + GAP,
                        1446,
                        1446 + GAP,
                    ])
            else:
                raise ValueError(f"Invalid version: {version}. Must be 'v01' or 'v05'")

        # set the layer radii which depends on phi_mod_2 and z_mod_2
        self.phi_stagger_dr = PHI_STAGGER_DR_WITH_Z_STAGGER if self.do_z_stagger else PHI_STAGGER_DR_WITHOUT_Z_STAGGER
        self.layer_radii_z_mod_2_eq_0 += (self.phi_stagger_dr * self.phi_mod_2)
        self.layer_radii_z_mod_2_eq_1 = self.layer_radii_z_mod_2_eq_0 + self.z_stagger_dr

        # set the number of layers
        self.n_layers = len(self.layer_radii_z_mod_2_eq_0)

        # set the double-layers
        self.layer_pairs = [
            (0, 1),
            (2, 3),
            (4, 5),
            (6, 7),
        ]

        # set the z-positions of the z-sensors [n_z_sensors]
        self.zpos = np.array([it*self.original_length for it in range(self.nz)])

        # set the z-offsets of the z-sensors
        if self.do_z_offsets:
            if self.tracker == "IT":
                if self.version in ["v01", "v05"]:
                    self.chosen_offsets = np.array([0, 1.5, 0, 1.8, 0, 1.8, 0, 1.8, 0, 2.3, 0, 2.8, 0, 3.2, 0, 3.4, 0, 3.4, 0, 3.4, 0, 3.4, 0])
                else:
                    raise ValueError(f"Invalid version: {self.version}. Must be 'v01' or 'v05'")
            else:
                if self.version in ["v01", "v05"]:
                    self.chosen_offsets = np.array([0, 1.2, 0, 1.3, 0, 1.6, 0, 2.1, 0, 2.6, 0, 2.8, 0, 2.8, 0, 2.8, 0, 3.1, 0, 3.3, 0])
                else:
                    raise ValueError(f"Invalid version: {self.version}. Must be 'v01' or 'v05'")
        else:
            self.chosen_offsets = np.array([0.0] * self.nz)

        # optimization settings
        self.possible_offsets = np.arange(0, 4, 0.1)
        if self.tracker == "IT":
            self.modules_under_test = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
        else:
            self.modules_under_test = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]

        # given all that, set the z-positions of the z-sensors
        self.set_z_positions()
        self.set_max_z_and_r()


    def set_z_positions(self):
        # zmins and zmaxs of the z-sensors [n_layers, n_z_sensors]
        self.z_mins = np.array([self.zpos + self.chosen_offsets for _ in range(self.n_layers)])
        self.z_maxs = self.z_mins + self.length


    def set_max_z_and_r(self):
        self.max_z = np.max(self.z_maxs)
        self.max_r = np.max(np.concatenate([self.layer_radii_z_mod_2_eq_0,
                                            self.layer_radii_z_mod_2_eq_1]))


    def describe(self):
        self.set_max_z_and_r()
        print(f"DetectorParameters:")
        print(f"  version: {self.version}")
        print(f"  tracker: {self.tracker}")
        print(f"  phi_mod_2: {self.phi_mod_2}")
        print(f"  do_z_stagger: {self.do_z_stagger}")
        print(f"  do_z_offsets: {self.do_z_offsets}")
        print(f"  original_length: {self.original_length}")
        print(f"  length: {self.length}")
        print(f"  nz: {self.nz}")
        print(f"  z_stagger_dr: {self.z_stagger_dr}")
        print(f"  phi_stagger_dr: {self.phi_stagger_dr}")
        print(f"  max_z: {self.max_z}")
        print(f"  max_r: {self.max_r}")
        # print(f"  layer_radii_z_mod_2_eq_0: {self.layer_radii_z_mod_2_eq_0}")
        # print(f"  layer_radii_z_mod_2_eq_1: {self.layer_radii_z_mod_2_eq_1}")
        # print(f"  layer_pairs: {self.layer_pairs}")
        # print(f"  zpos: {self.zpos}")
        # print(f"  chosen_offsets: {self.chosen_offsets}")
        # print(f"  z_mins: {self.z_mins}")
        # print(f"  z_maxs: {self.z_maxs}")


    def make_xml_blurbs(self):
        self.xml_blurbs = []
        self.set_xml_params()
        self.make_xml_constants_blurbs()
        self.make_xml_layers_blurbs()


    def set_xml_params(self):
        self.xml = {}
        if self.version == "v01":
            if self.tracker == "IT":
                self.xml["nzs"] = [32, 32, 32, 32, 46, 46, 46, 46]
                self.xml["nphis"] = [15*2, 15*2, 20*2, 20*2, 58*2, 58*2, 62*2, 62*2]
                self.xml["half_length"] = [481.6, 481.6, 481.6, 481.6, 692.3, 692.3, 692.3, 692.3]
            elif self.tracker == "OT":
                self.xml["nzs"] = [42] * self.n_layers
                self.xml["nphis"] = [48*2, 48*2, 52*2, 52*2, 80*2, 80*2, 84*2, 84*2]
                self.xml["half_length"] = [1264.2] * self.n_layers
        elif self.version == "v05":
            if self.tracker == "IT":
                self.xml["nzs"] = [32, 32, 32, 32, 32, 32, 46, 46]
                self.xml["nphis"] = [15*2, 15*2, 30*2, 30*2, 46*2, 46*2, 62*2, 62*2]
                self.xml["half_length"] = [481.6, 481.6, 481.6, 481.6, 481.6, 481.6, 692.3, 692.3]
            elif self.tracker == "OT":
                self.xml["nzs"] = [42] * self.n_layers
                self.xml["nphis"] = [48*2, 48*2, 60*2, 60*2, 72*2, 72*2, 84*2, 84*2]
                self.xml["half_length"] = [1264.2] * self.n_layers

        # +z sensors are mirrored to -z, so we only need to consider half of them
        self.xml["nzs"] = [nz // 2 for nz in self.xml["nzs"]]


    def make_xml_constants_blurbs(self):
        """
        <constant name="InnerTracker_Barrel_radius_0" value="127*mm"/>
        <constant name="InnerTracker_Barrel_radius_1" value="167*mm"/>
        <constant name="InnerTracker_Barrel_radius_2" value="510*mm"/>
        <constant name="InnerTracker_Barrel_radius_3" value="550*mm"/>
        """
        tracker = "InnerTracker" if self.tracker == "IT" else "OuterTracker"
        self.xml_names = {}
        self.xml_names["gap"] = tracker + "_Barrel_DoubleLayer_Gap"

        # layer radius (without z-stagger dr)
        # layer radius (with z-stagger dr)
        for lower, upper in self.layer_pairs:
            doubler_layer = lower // 2
            self.xml_names[f"layer_radii_z_mod_2_eq_0_{doubler_layer}"] = f"{tracker}_Barrel_radius_{doubler_layer}"
            self.xml_names[f"layer_radii_z_mod_2_eq_1_{doubler_layer}"] = f"{tracker}_Barrel_radius_{doubler_layer}_z_staggered"
            self.xml_blurbs.append(xml_constant_mm_template().format(name=self.xml_names[f"layer_radii_z_mod_2_eq_0_{doubler_layer}"],
                                                                     value=self.layer_radii_z_mod_2_eq_0[lower]))
            self.xml_blurbs.append(xml_constant_mm_template().format(name=self.xml_names[f"layer_radii_z_mod_2_eq_1_{doubler_layer}"],
                                                                     value=self.layer_radii_z_mod_2_eq_1[lower]))

        # original sensor length
        # sensor length
        self.xml_names["sensor_length_original"] = tracker + "_Barrel_OriginalSensorLength"
        self.xml_names["sensor_length"] = tracker + "_Barrel_SensorLength"
        self.xml_blurbs.append(xml_constant_mm_template().format(name=self.xml_names["sensor_length_original"],
                                                                 value=self.original_length))
        self.xml_blurbs.append(xml_constant_mm_template().format(name=self.xml_names["sensor_length"],
                                                                 value=self.length))


    def make_xml_layers_blurbs(self):

        bits_zsensor = 8
        bits_layer = 5
        bits_total = bits_zsensor + bits_layer
        id_max = 2**bits_total - 1

        layer_module = ("InnerTrackerBarrelModule_01"
                        if self.tracker == "IT" else
                        "OuterTrackerBarrelModule_In")

        for layer in range(self.n_layers):

            layer_dicts = []

            # step 1: make a dict for each barrel ring
            for iz in range(self.xml["nzs"][layer]):
                z_mod_2 = int(iz % 2)
                offset = self.chosen_offsets[iz]
                layer_id = layer
                nphi = self.xml["nphis"][layer]
                radius = self.xml_names[f"layer_radii_z_mod_2_eq_{z_mod_2}_{int(layer // 2)}"]
                if layer % 2 == 1:
                    radius += f" + {self.xml_names['gap']}"
                dr = self.phi_stagger_dr
                z0 = iz * self.original_length + offset + self.length / 2
                nz = 1
                dic = dict(
                    layer_module=layer_module,
                    layer_id=layer_id,
                    nphi=nphi,
                    radius=radius,
                    dr=dr,
                    z0=z0,
                    nz=nz,
                )

                # positive z
                layer_dicts.append(dic)

                # negative z
                ndic = dic.copy()
                ndic["z0"] = -1.0 * dic["z0"]
                layer_dicts.append(ndic)

            # step 2: sort dicts by z0 and label them by increasing z0
            layer_dicts.sort(key=lambda dic: dic["z0"], reverse=True)
            for it in range(len(layer_dicts)):
                layer_id = (layer << bits_zsensor) + it
                if layer_id > id_max:
                    raise ValueError(f"Layer ID {layer_id} exceeds maximum {id_max}")
                layer_dicts[it]["layer_id"] = layer_id

            # step 3: fin
            for dic in layer_dicts:
                self.xml_blurbs.append(xml_layer_template().format(**dic))


def xml_constant_mm_template():
    """
    <constant name="InnerTracker_Barrel_radius_0" value="127*mm"/>
    <constant name="InnerTracker_Barrel_radius_1" value="167*mm"/>
    <constant name="InnerTracker_Barrel_radius_2" value="510*mm"/>
    <constant name="InnerTracker_Barrel_radius_3" value="550*mm"/>
    """
    return '<constant name="{name}" value="{value}*mm"/>\n'


def xml_layer_template():
    return (
"""
<layer module="{layer_module}" id="{layer_id}">
    <rphi_layout phi_tilt="0*deg" nphi="{nphi}" phi0="0" rc="{radius}" dr="{dr}*mm"/>
    <z_layout dr="0" z0="{z0:.3f}*mm" nz="{nz}"/>
</layer>
""")


def eta_to_theta(eta: float) -> float:
    return 2 * np.arctan(np.exp(-eta))


def options():
    parser = argparse.ArgumentParser(description="Demo for 2D scatter plot")
    parser.add_argument("--stagger", action="store_true", help="Enable z-staggering for the demo")
    parser.add_argument("--optimize", action="store_true", help="Enable optimization for the demo")
    parser.add_argument("--offsets", action="store_true", help="Enable offsets for the demo")
    parser.add_argument("--plot", action="store_true", help="Enable plotting for the demo")
    parser.add_argument("--phi-mod-2", type=int, default=None, choices=[0, 1], help="Specify the phi_mod_2 for the demo")
    parser.add_argument("--tracker", type=str, default=None, choices=["IT", "OT"], help="Specify the tracker for the demo")
    parser.add_argument("--version", type=str, default=None, choices=["v01", "v05"], help="Specify the version for the demo")
    parser.add_argument("--xml", type=str, default="tmp.xml", help="Filename of output xml")
    return parser.parse_args()


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
