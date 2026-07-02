"""
Steering file for counting doublets in a LST-friendly MAIA detector
"""
import argparse
from glob import glob
import os
import pandas as pd
import time
import logging
logger = logging.getLogger(__name__)

from maia_doublets.constants import SIGNAL, NICKNAME_TO_SYSTEM
from maia_doublets.datasets import get_filepaths, parse_filepaths
from maia_doublets.slcio import HitMaker
from maia_doublets.md import MDMaker
from maia_doublets.t2 import T2Maker
from maia_doublets.t4 import T4Maker
from maia_doublets.t8 import T8Maker
from maia_doublets.plot import Plotter
from maia_doublets.calib import CalibConstants
from maia_doublets.calib import MDCalibrator
from maia_doublets.calib import T2Calibrator
from maia_doublets.calib import T4Calibrator
from maia_doublets.calib import T8Calibrator


def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    # parse options
    ops = options()
    valid_geos = ["v01", "v04", "v05"]
    valid_smears = ["00um", "10um"]
    if ops.geo not in valid_geos:
        raise ValueError(f"Invalid geometry version specified, must be one of {valid_geos}")
    if ops.smear not in valid_smears:
        raise ValueError(f"Invalid smear value specified, must be one of {valid_smears}")
    if ops.i:
        fnames = parse_filepaths(ops.i)
    else:
        fnames = get_filepaths(
            geometry_version=ops.geo,
            signal=ops.signal,
            background10=ops.background10,
            background100=ops.background100,
            sim=ops.sim,
            digi=ops.digi,
            smear=ops.smear,
        )
    if not fnames:
        raise ValueError("No input files found")
    if not ops.layers:
        raise ValueError("At least one layer must be specified")
    layers = parse_layers(ops.layers)
    geometry = ops.geometry
    signal = ops.signal or any(SIGNAL in os.path.basename(fname) for fname in fnames)
    cut_mds = ops.cut_mds or not signal
    cut_t2s = ops.cut_t2s or not signal
    cut_t4s = ops.cut_t4s or not signal
    cut_t8s = ops.cut_t8s or not signal
    if not ops.sim and not ops.digi:
        raise ValueError("At least one of --sim or --digi must be specified")
    if ops.sim and ops.digi:
        raise ValueError("Only one of --sim or --digi can be specified, not both")
    if ops.calibrate and (cut_mds or cut_t2s or cut_t4s or cut_t8s):
        raise ValueError("Cannot use --calibrate with any of --cut-mds, --cut-t2s, --cut-t4s, or --cut-t8s")

    # log some info
    logger.info(f"Detected {'signal' if signal else 'background'} files")
    logger.info(f"Found {len(fnames)} files")
    logger.info(f"Layers provided: {ops.layers}")
    logger.info(f"Layers decoded: {layers}")
    logger.info(f"Do calibration: {ops.calibrate}")
    logger.info(f"Calib json: {calib_json(ops)}")
    logger.info(f"Cut MDs: {cut_mds}")
    logger.info(f"Cut T2s: {cut_t2s}")
    logger.info(f"Cut T4s: {cut_t4s}")
    logger.info(f"Cut T8s: {cut_t8s}")
    logger.info(f"Fast MDs: {ops.fast_mds}")
    logger.info(f"Geometry version: {ops.geo}")
    logger.info(f"Using sim hits: {ops.sim}")
    logger.info(f"Using digi hits: {ops.digi}")
    if ops.digi:
        logger.info(f"Smear value for digi hits: {ops.smear}")

    # calib constants
    calibs = CalibConstants(calib_json(ops)).calibs

    # simhits and mcparticles
    simhits, mcps, hit_time = get_simhits_and_mcps(ops, fnames, geometry, signal, layers)
    write_simhits_and_mcps(ops, simhits, mcps)

    # mini-doublets (mds)
    doublets, md_time = get_mds(ops, simhits, signal, cut_mds, calibs)
    write_mds(ops, doublets)
    if ops.calibrate:
        calib_mds(ops, doublets)
        calibs = CalibConstants(calib_json(ops)).calibs
        doublets, md_time = get_mds(ops, simhits, signal, cut_mds, calibs)

    # t2s
    t2s, t2_time = get_t2s(ops, doublets, signal, cut_t2s, calibs)
    write_t2s(ops, t2s)
    if ops.calibrate:
        calib_t2s(ops, t2s)
        calibs = CalibConstants(calib_json(ops)).calibs
        t2s, t2_time = get_t2s(ops, doublets, signal, cut_t2s, calibs)


    # t4s
    t4s, t4_time = get_t4s(ops, t2s, signal, cut_t4s, calibs)
    write_t4s(ops, t4s)
    if ops.calibrate:
        calib_t4s(ops, t4s)
        calibs = CalibConstants(calib_json(ops)).calibs
        t4s, t4_time = get_t4s(ops, t2s, signal, cut_t4s, calibs)

    # t8s
    t8s, t8_time = get_t8s(ops, t4s, signal, cut_t8s, calibs)
    write_t8s(ops, t8s)
    if ops.calibrate:
        calib_t8s(ops, t8s)
        calibs = CalibConstants(calib_json(ops)).calibs
        t8s, t8_time = get_t8s(ops, t4s, signal, cut_t8s, calibs)

    # plot stuff
    with Timer() as plot_time:
        if ops.plot:
            logger.info("Creating plots ...")
            plotter = Plotter(
                signal=signal,
                mcps=mcps,
                simhits=simhits,
                doublets=doublets,
                linesegments=t2s,
                t4s=t4s,
                t8s=t8s,
                calibs=calibs,
                pdf="doublets.pdf",
            )
            plotter.plot()

    # log timing info
    logger.info(f"Timing info (in seconds):")
    logger.info(f"  Hit making: {hit_time:.2f}")
    logger.info(f"  MD making: {md_time:.2f}")
    logger.info(f"  T2 making: {t2_time:.2f}")
    logger.info(f"  T4 making: {t4_time:.2f}")
    logger.info(f"  T8 making: {t8_time:.2f}")
    logger.info(f"  Plotting: {plot_time.duration:.2f}")



def get_simhits_and_mcps(
    ops: argparse.Namespace,
    fnames: list[str],
    geometry: bool,
    signal: bool,
    layers: dict[int, set[int]]
) -> tuple[pd.DataFrame, pd.DataFrame, float]:

    with Timer() as hit_time:
        if ops.read_mcps and ops.read_simhits:
            logger.info(f"Reading simhits {ops.read_simhits} and mcps {ops.read_mcps} ...")
            mcps = pd.read_pickle(ops.read_mcps)
            simhits = pd.read_pickle(ops.read_simhits)
        elif any([
            ops.read_mcps and not ops.read_simhits,
            ops.read_simhits and not ops.read_mcps,
        ]):
            raise ValueError("Both --read-mcps and --read-simhits must be specified together")
        else:
            # convert slcio to hits dataframe
            converter = HitMaker(slcio_file_paths=fnames,
                                load_geometry=geometry,
                                signal=signal,
                                sim=ops.sim,
                                layers=layers,
                                )
            mcps, simhits = converter.convert()

    return simhits, mcps, hit_time.duration


def write_simhits_and_mcps(ops: argparse.Namespace, simhits: pd.DataFrame, mcps: pd.DataFrame) -> None:
    if ops.write_mcps:
        logger.info(f"Saving mcps to {ops.write_mcps} ...")
        mcps.to_pickle(ops.write_mcps)
    if ops.write_simhits:
        logger.info(f"Saving simhits to {ops.write_simhits} ...")
        simhits.to_pickle(ops.write_simhits)


def get_mds(ops: argparse.Namespace, simhits: pd.DataFrame, signal: bool, cut_mds: bool, calibs: dict) -> tuple[pd.DataFrame, float]:
    with Timer() as md_time:
        if ops.read_mds:
            logger.info(f"Reading mini-doublets from {ops.read_mds} ...")
            doublets = pd.read_pickle(ops.read_mds)
        else:
            # make mini-doublets from hits
            doublets = None
            doublets = MDMaker(
                signal=signal,
                cut_mds=cut_mds,
                fast_merge=ops.fast_mds,
                calibs=calibs,
                simhits=simhits,
            ).df

    return doublets, md_time.duration


def write_mds(ops: argparse.Namespace, doublets: pd.DataFrame) -> None:
    if ops.write_mds:
        logger.info(f"Saving mini-doublets to {ops.write_mds} ...")
        doublets.to_pickle(ops.write_mds)


def calib_mds(ops: argparse.Namespace, doublets: pd.DataFrame) -> None:
    if not ops.calibrate:
        return
    logger.info("Calibrating MDs ...")
    calib = MDCalibrator(doublets, calib_json=calib_json(ops))
    calib.calibrate()


def get_t2s(ops: argparse.Namespace, doublets: pd.DataFrame, signal: bool, cut_t2s: bool, calibs: dict) -> tuple[pd.DataFrame, float]:
    with Timer() as t2_time:
        if ops.read_t2s:
            logger.info(f"Reading T2s (line segments) from {ops.read_t2s} ...")
            t2s = pd.read_pickle(ops.read_t2s)
        else:
            # make T2s (line segments) from mini-doublets
            t2s = None
            t2s = T2Maker(
                signal=signal,
                cut_t2s=cut_t2s,
                calibs=calibs,
                doublets=doublets,
            ).df

    return t2s, t2_time.duration


def write_t2s(ops: argparse.Namespace, t2s: pd.DataFrame) -> None:
    if ops.write_t2s:
        logger.info(f"Saving T2s (line segments) to {ops.write_t2s} ...")
        t2s.to_pickle(ops.write_t2s)


def calib_t2s(ops: argparse.Namespace, t2s: pd.DataFrame) -> None:
    if not ops.calibrate:
        return
    logger.info("Calibrating T2s ...")
    calib = T2Calibrator(t2s, calib_json=calib_json(ops))
    calib.calibrate()


def get_t4s(ops: argparse.Namespace, t2s: pd.DataFrame, signal: bool, cut_t4s: bool, calibs: dict) -> tuple[pd.DataFrame, float]:
    with Timer() as t4_time:
        if ops.read_t4s:
            logger.info(f"Reading T4s from {ops.read_t4s} ...")
            t4s = pd.read_pickle(ops.read_t4s)
        else:
            # make T4s from T2s (line segments)
            t4s = None
            t4s = T4Maker(
                signal=signal,
                cut_t4s=cut_t4s,
                calibs=calibs,
                t2s=t2s,
            ).df

    return t4s, t4_time.duration


def write_t4s(ops: argparse.Namespace, t4s: pd.DataFrame) -> None:
    if ops.write_t4s:
        logger.info(f"Saving T4s to {ops.write_t4s} ...")
        t4s.to_pickle(ops.write_t4s)


def calib_t4s(ops: argparse.Namespace, t4s: pd.DataFrame) -> None:
    if not ops.calibrate:
        return
    logger.info("Calibrating T4s ...")
    calib = T4Calibrator(t4s, calib_json=calib_json(ops))
    calib.calibrate()


def get_t8s(ops: argparse.Namespace, t4s: pd.DataFrame, signal: bool, cut_t8s: bool, calibs: dict) -> tuple[pd.DataFrame, float]:
    with Timer() as t8_time:
        if ops.read_t8s:
            logger.info(f"Reading T8s from {ops.read_t8s} ...")
            t8s = pd.read_pickle(ops.read_t8s)
        else:
            # make T8s from T4s
            t8s = None
            t8s = T8Maker(
                signal=signal,
                cut_t8s=cut_t8s,
                calibs=calibs,
                t4s=t4s,
            ).df

    return t8s, t8_time.duration


def write_t8s(ops: argparse.Namespace, t8s: pd.DataFrame) -> None:
    if ops.write_t8s:
        logger.info(f"Saving T8s to {ops.write_t8s} ...")
        t8s.to_pickle(ops.write_t8s)


def calib_t8s(ops: argparse.Namespace, t8s: pd.DataFrame) -> None:
    if not ops.calibrate:
        return
    logger.info("Calibrating T8s ...")
    calib = T8Calibrator(t8s, calib_json=calib_json(ops))
    calib.calibrate()


def calib_json(ops: argparse.Namespace) -> str:
    key = (ops.geo, "sim") if ops.sim else (ops.geo, "digi", ops.smear)
    key = "_".join(key)
    return os.path.join(ops.calib_dir, f"{key}.json")


def guess_calib_dir() -> str:
    calib_guess = glob("../*/calibs/") or [""]
    return calib_guess[0]


def options():
    preset = [
        "ITB0", "ITB1", "ITB2", "ITB3",
        "ITB4", "ITB5", "ITB6", "ITB7",
        "OTB0", "OTB1", "OTB2", "OTB3",
        "OTB4", "OTB5", "OTB6", "OTB7",
    ]
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", default=[], help="Input slcio file or glob pattern")
    parser.add_argument("--calibrate", action="store_true", help="Measure and write calibration constants (signal intervals) to file")
    parser.add_argument("--calib-dir", type=str, default=guess_calib_dir(), help="Directory of calibration constants")
    parser.add_argument("--geometry", action="store_true", help="Load compact geometry from xml")
    parser.add_argument("--layers", nargs="+", type=str, default=preset, help="List of layers to consider")
    parser.add_argument("--sim", action="store_true", help="Use sim hits in the analysis")
    parser.add_argument("--digi", action="store_true", help="Use digi hits in the analysis")
    parser.add_argument("--plot", action="store_true", help="Include plots in the analysis")
    parser.add_argument("--cut-mds", action="store_true", help="Cut MDs based on MD_DZ_CUT and MD_DR_CUT")
    parser.add_argument("--cut-t2s", action="store_true", help="Cut T2s (line segments) based on [[ something ]]")
    parser.add_argument("--cut-t4s", action="store_true", help="Cut T4s based on [[ something ]]")
    parser.add_argument("--cut-t8s", action="store_true", help="Cut T8s based on [[ something ]]")
    parser.add_argument("--read-mcps", type=str, help="Read mcps from pickle file")
    parser.add_argument("--write-mcps", type=str, help="Write mcps to pickle file")
    parser.add_argument("--read-simhits", type=str, help="Read simhits from pickle file")
    parser.add_argument("--write-simhits", type=str, help="Write simhits to pickle file")
    parser.add_argument("--read-mds", type=str, help="Read mini-doublets from pickle file")
    parser.add_argument("--write-mds", type=str, help="Write mini-doublets to pickle file")
    parser.add_argument("--fast-mds", action="store_true", help="Use fast binned merge for mini-doublets")
    parser.add_argument("--read-t2s", type=str, help="Read T2s (line segments) from pickle file")
    parser.add_argument("--write-t2s", type=str, help="Write T2s (line segments) to pickle file")
    parser.add_argument("--read-t4s", type=str, help="Read T4s from pickle file")
    parser.add_argument("--write-t4s", type=str, help="Write T4s to pickle file")
    parser.add_argument("--read-t8s", type=str, help="Read T8s from pickle file")
    parser.add_argument("--write-t8s", type=str, help="Write T8s to pickle file")
    parser.add_argument("--geo", type=str, help="Version of geometry to use for cuts (e.g. v01, v04)", required=True)
    parser.add_argument("--smear", type=str, default="00um", help="Smear value to use for digi hits (e.g. 10um)")
    parser.add_argument("--signal", action="store_true", help="Use signal files in the analysis")
    parser.add_argument("--background10", action="store_true", help="Use background files (10 percent) in the analysis")
    parser.add_argument("--background100", action="store_true", help="Use background files (100 percent) in the analysis")
    parser.add_argument("--debug", action="store_true", help="Print some debug information")
    return parser.parse_args()


def parse_layers(layers_str_list: list[str]) -> dict[int, set[int]]:
    """
    Parse layers like ITB4, OTB3, etc. into a dict of {system: set of layers}
    e.g. ["ITB4", "OTB3"] -> {INNER_TRACKER_BARREL: {4}, OUTER_TRACKER_BARREL: {3}}
    """
    dict_of_system_layer_pairs = {}
    for layer_str in layers_str_list:
        if len(layer_str) != len("ITB4"):
            raise ValueError(f"Invalid layer specified: {layer_str}")
        layer = int(layer_str[-1])
        system_str = layer_str[:-1]
        system = parse_system(system_str)
        if system not in dict_of_system_layer_pairs:
            dict_of_system_layer_pairs[system] = set()
        dict_of_system_layer_pairs[system].add(layer)
    return dict_of_system_layer_pairs


def parse_system(system_str: str) -> int:
    if len(system_str) != len("OTB"):
        raise ValueError(f"Invalid system specified: {system_str}")
    return NICKNAME_TO_SYSTEM[system_str]


class Timer:
    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, *args):
        self.end = time.perf_counter()
        self.duration = self.end - self.start


if __name__ == "__main__":
    main()
