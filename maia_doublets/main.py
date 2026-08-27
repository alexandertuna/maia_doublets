"""
Steering file for counting doublets in a LST-friendly MAIA detector
"""
import argparse
from glob import glob
import ndjson
import os
import pandas as pd
import time
import logging
logger = logging.getLogger(__name__)

from maia_doublets.constants import MUONGUN, PIONGUN, NICKNAME_TO_SYSTEM
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
    check_options(ops)
    if ops.i:
        fnames = parse_filepaths(ops.i)
    else:
        fnames = get_filepaths(
            geometry_version=ops.geo,
            dataset=ops.dataset,
            sim=ops.sim,
            digi=ops.digi,
            smear=ops.smear,
        )
    if not fnames:
        raise ValueError("No input files found")
    layers = parse_layers(ops.layers)
    data_source = ops.dataset
    signal = (any(MUONGUN in os.path.basename(fname) for fname in fnames) or
              any(PIONGUN in os.path.basename(fname) for fname in fnames))
    pdf = ops.pdf or f"{calib_key(ops)}_{data_source}.pdf"
    geo_version = ops.geo
    cut_mds = ops.cut_mds or not signal
    cut_t2s = ops.cut_t2s or not signal
    cut_t4s = ops.cut_t4s or not signal
    cut_t8s = ops.cut_t8s or not signal
    if ops.no_cuts:
        cut_mds = cut_t2s = cut_t4s = cut_t8s = False
    if ops.calibrate and (cut_mds or cut_t2s or cut_t4s or cut_t8s):
        raise ValueError("Cannot use --calibrate with any of --cut-mds, --cut-t2s, --cut-t4s, or --cut-t8s")

    # log some info
    logger.info(f"Detected {data_source} files")
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
    logger.info(f"Writing ndjson: {ops.cutflow}")
    logger.info(f"Writing pdf: {pdf}")
    if ops.digi:
        logger.info(f"Smear value for digi hits: {ops.smear}")

    # calib constants
    calibs = CalibConstants(calib_json(ops)).calibs

    # hits and mcparticles
    hits, mcps, hit_cutflow, hit_time = get_hits_and_mcps(ops, fnames, geo_version, signal, layers)
    write_hits_and_mcps(ops, hits, mcps, hit_cutflow)
    if ops.stop_after_hits:
        logger.info("Stopping after hits, as requested")
        return

    # mini-doublets (mds)
    mds, md_cutflow, md_time = get_mds(ops, hits, signal, cut_mds, calibs)
    write_mds(ops, mds, md_cutflow)
    if ops.calibrate:
        calib_mds(ops, mds)
        calibs = CalibConstants(calib_json(ops)).calibs
        mds, md_cutflow, md_time = get_mds(ops, hits, signal, cut_mds, calibs)
    if ops.stop_after_mds:
        logger.info("Stopping after MDs, as requested")
        return

    # t2s
    t2s, t2_cutflow, t2_time = get_t2s(ops, mds, signal, cut_t2s, calibs)
    write_t2s(ops, t2s, t2_cutflow)
    if ops.calibrate:
        calib_t2s(ops, t2s)
        calibs = CalibConstants(calib_json(ops)).calibs
        t2s, t2_cutflow, t2_time = get_t2s(ops, mds, signal, cut_t2s, calibs)

    # t4s
    t4s, t4_cutflow, t4_time = get_t4s(ops, t2s, signal, cut_t4s, calibs)
    write_t4s(ops, t4s, t4_cutflow)
    if ops.calibrate:
        calib_t4s(ops, t4s)
        calibs = CalibConstants(calib_json(ops)).calibs
        t4s, t4_cutflow, t4_time = get_t4s(ops, t2s, signal, cut_t4s, calibs)

    # t8s
    t8s, t8_cutflow, t8_time = get_t8s(ops, t4s, signal, cut_t8s, calibs)
    write_t8s(ops, t8s, t8_cutflow)
    if ops.calibrate:
        calib_t8s(ops, t8s)
        calibs = CalibConstants(calib_json(ops)).calibs
        t8s, t8_cutflow, t8_time = get_t8s(ops, t4s, signal, cut_t8s, calibs)

    # plot stuff
    with Timer() as plot_time:
        if ops.plot:
            logger.info("Creating plots ...")
            plotter = Plotter(
                signal=signal,
                mcps=mcps,
                hits=hits,
                mds=mds,
                t2s=t2s,
                t4s=t4s,
                t8s=t8s,
                calibs=calibs,
                pdf=pdf,
            )
            plotter.plot()

    # write cutflows
    if ops.cutflow:
        logger.info(f"Writing cutflows to {ops.cutflow} ...")
        with open(ops.cutflow, "w") as fi:
            ndjson.dump([
                hit_cutflow.to_dict(orient="records"),
                md_cutflow.to_dict(orient="records"),
                t2_cutflow.to_dict(orient="records"),
                t4_cutflow.to_dict(orient="records"),
                t8_cutflow.to_dict(orient="records"),
            ], fi)

    # log timing info
    logger.info(f"Timing info (in seconds):")
    logger.info(f"  Hit making: {hit_time:.2f}")
    logger.info(f"  MD making: {md_time:.2f}")
    logger.info(f"  T2 making: {t2_time:.2f}")
    logger.info(f"  T4 making: {t4_time:.2f}")
    logger.info(f"  T8 making: {t8_time:.2f}")
    logger.info(f"  Plotting: {plot_time.duration:.2f}")


def check_options(ops: argparse.Namespace) -> None:
    valid_geos = ["v01", "v04", "v05", "v06", "v07"]
    valid_smears = ["00um", "05um", "10um", "20um"]
    if ops.geo not in valid_geos:
        raise ValueError(f"Invalid geometry version specified, must be one of {valid_geos}")
    if ops.smear not in valid_smears:
        raise ValueError(f"Invalid smear value specified, must be one of {valid_smears}")
    if not ops.sim and not ops.digi:
        raise ValueError("At least one of --sim or --digi must be specified")
    if ops.sim and ops.digi:
        raise ValueError("Only one of --sim or --digi can be specified, not both")
    if not ops.layers:
        raise ValueError("At least one layer must be specified")
    if ops.write_hits and not ops.write_hits.endswith(".pkl"):
        raise ValueError("Output file for --write-hits must end with .pkl")
    if ops.write_mcps and not ops.write_mcps.endswith(".pkl"):
        raise ValueError("Output file for --write-mcps must end with .pkl")
    if ops.write_mds and not ops.write_mds.endswith(".pkl"):
        raise ValueError("Output file for --write-mds must end with .pkl")
    if ops.write_t2s and not ops.write_t2s.endswith(".pkl"):
        raise ValueError("Output file for --write-t2s must end with .pkl")
    if ops.write_t4s and not ops.write_t4s.endswith(".pkl"):
        raise ValueError("Output file for --write-t4s must end with .pkl")


def cutflow_path(df_path: str) -> str:
    return df_path.replace(".pkl", ".json")


def get_hits_and_mcps(
    ops: argparse.Namespace,
    fnames: list[str],
    geo_version: str,
    signal: bool,
    layers: dict[int, set[int]]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, float]:

    with Timer() as hit_time:
        if ops.read_mcps and ops.read_hits:
            logger.info(f"Reading hits {ops.read_hits} and mcps {ops.read_mcps} ...")
            mcps = pd.read_pickle(ops.read_mcps)
            hits = pd.read_pickle(ops.read_hits)
            cutflow = pd.read_json(cutflow_path(ops.read_hits))
        elif any([
            ops.read_mcps and not ops.read_hits,
            ops.read_hits and not ops.read_mcps,
        ]):
            raise ValueError("Both --read-mcps and --read-hits must be specified together")
        else:
            # convert slcio to hits dataframe
            converter = HitMaker(slcio_file_paths=fnames,
                                geo_version=geo_version,
                                signal=signal,
                                sim=ops.sim,
                                layers=layers,
                                )
            mcps, hits, cutflow = converter.convert()

    return hits, mcps, cutflow, hit_time.duration


def write_hits_and_mcps(ops: argparse.Namespace, hits: pd.DataFrame, mcps: pd.DataFrame, cutflow: pd.DataFrame) -> None:
    if ops.write_mcps:
        logger.info(f"Saving mcps to {ops.write_mcps} ...")
        mcps.to_pickle(ops.write_mcps)
    if ops.write_hits:
        jname = cutflow_path(ops.write_hits)
        logger.info(f"Saving hits to {ops.write_hits} and cutflow to {jname} ...")
        hits.to_pickle(ops.write_hits)
        cutflow.to_json(jname, indent=4)


def get_mds(ops: argparse.Namespace, hits: pd.DataFrame, signal: bool, cut_mds: bool, calibs: dict) -> tuple[pd.DataFrame,
                                                                                                                pd.DataFrame,
                                                                                                                float]:
    with Timer() as md_time:
        if ops.read_mds:
            cpath = cutflow_path(ops.read_mds)
            logger.info(f"Reading mini-doublets from {ops.read_mds} and cutflow from {cpath} ...")
            doublets = pd.read_pickle(ops.read_mds)
            cutflow = pd.read_json(cpath)
        else:
            # make mini-doublets from hits
            doublets = None
            maker = MDMaker(
                signal=signal,
                cut_mds=cut_mds,
                fast_merge=ops.fast_mds,
                calibs=calibs,
                hits=hits,
            )
            doublets = maker.df
            cutflow = maker.cutflow

    return doublets, cutflow, md_time.duration


def write_mds(ops: argparse.Namespace, doublets: pd.DataFrame, cutflow: pd.DataFrame) -> None:
    if not ops.write_mds:
        return
    jname = ops.write_mds.replace(".pkl", ".json")
    logger.info(f"Saving mini-doublets to {ops.write_mds} and cutflow to {jname} ...")
    doublets.to_pickle(ops.write_mds)
    cutflow.to_json(jname, indent=4)


def calib_mds(ops: argparse.Namespace, doublets: pd.DataFrame) -> None:
    if not ops.calibrate:
        return
    logger.info("Calibrating MDs ...")
    calib = MDCalibrator(doublets, calib_json=calib_json(ops))
    calib.calibrate()


def get_t2s(ops: argparse.Namespace, mds: pd.DataFrame, signal: bool, cut_t2s: bool, calibs: dict) -> tuple[pd.DataFrame,
                                                                                                            pd.DataFrame,
                                                                                                            float]:
    with Timer() as t2_time:
        if ops.read_t2s:
            cpath = cutflow_path(ops.read_t2s)
            logger.info(f"Reading T2s from {ops.read_t2s} and cutflow from {cpath} ...")
            t2s = pd.read_pickle(ops.read_t2s)
        else:
            # make T2s from mini-doublets
            t2s = None
            maker = T2Maker(
                signal=signal,
                cut_t2s=cut_t2s,
                calibs=calibs,
                mds=mds,
            )
            t2s = maker.df
            cutflow = maker.cutflow

    return t2s, cutflow, t2_time.duration


def write_t2s(ops: argparse.Namespace, t2s: pd.DataFrame, cutflow: pd.DataFrame) -> None:
    if not ops.write_t2s:
        return
    cpath = cutflow_path(ops.write_t2s)
    logger.info(f"Saving T2s to {ops.write_t2s} and cutflow to {cpath} ...")
    t2s.to_pickle(ops.write_t2s)
    cutflow.to_json(cpath, indent=4)


def calib_t2s(ops: argparse.Namespace, t2s: pd.DataFrame) -> None:
    if not ops.calibrate:
        return
    logger.info("Calibrating T2s ...")
    calib = T2Calibrator(t2s, calib_json=calib_json(ops))
    calib.calibrate()


def get_t4s(ops: argparse.Namespace, t2s: pd.DataFrame, signal: bool, cut_t4s: bool, calibs: dict) -> tuple[pd.DataFrame,
                                                                                                            pd.DataFrame,
                                                                                                            float]:
    with Timer() as t4_time:
        if ops.read_t4s:
            logger.info(f"Reading T4s from {ops.read_t4s} ...")
            t4s = pd.read_pickle(ops.read_t4s)
        else:
            # make T4s from T2s
            t4s = None
            maker = T4Maker(
                signal=signal,
                cut_t4s=cut_t4s,
                calibs=calibs,
                t2s=t2s,
            )
            t4s = maker.df
            cutflow = maker.cutflow

    return t4s, cutflow, t4_time.duration


def write_t4s(ops: argparse.Namespace, t4s: pd.DataFrame, cutflow: pd.DataFrame) -> None:
    if not ops.write_t4s:
        return
    cpath = cutflow_path(ops.write_t4s)
    logger.info(f"Saving T4s to {ops.write_t4s} and cutflow to {cpath} ...")
    t4s.to_pickle(ops.write_t4s)
    cutflow.to_json(cpath, indent=4)


def calib_t4s(ops: argparse.Namespace, t4s: pd.DataFrame) -> None:
    if not ops.calibrate:
        return
    logger.info("Calibrating T4s ...")
    calib = T4Calibrator(t4s, calib_json=calib_json(ops))
    calib.calibrate()


def get_t8s(ops: argparse.Namespace, t4s: pd.DataFrame, signal: bool, cut_t8s: bool, calibs: dict) -> tuple[pd.DataFrame,
                                                                                                            pd.DataFrame,
                                                                                                            float]:
    with Timer() as t8_time:
        if ops.read_t8s:
            logger.info(f"Reading T8s from {ops.read_t8s} ...")
            t8s = pd.read_pickle(ops.read_t8s)
        else:
            # make T8s from T4s
            t8s = None
            maker = T8Maker(
                signal=signal,
                cut_t8s=cut_t8s,
                calibs=calibs,
                t4s=t4s,
            )
            t8s = maker.df
            cutflow = maker.cutflow

    return t8s, cutflow, t8_time.duration


def write_t8s(ops: argparse.Namespace, t8s: pd.DataFrame, cutflow: pd.DataFrame) -> None:
    if not ops.write_t8s:
        return
    cpath = cutflow_path(ops.write_t8s)
    logger.info(f"Saving T8s to {ops.write_t8s} and cutflow to {cpath} ...")
    t8s.to_pickle(ops.write_t8s)
    cutflow.to_json(cpath, indent=4)


def calib_t8s(ops: argparse.Namespace, t8s: pd.DataFrame) -> None:
    if not ops.calibrate:
        return
    logger.info("Calibrating T8s ...")
    calib = T8Calibrator(t8s, calib_json=calib_json(ops))
    calib.calibrate()


def calib_key(ops: argparse.Namespace) -> str:
    key = (ops.geo, "sim") if ops.sim else (ops.geo, "digi", ops.smear)
    key = "_".join(key)
    return key


def calib_json(ops: argparse.Namespace) -> str:
    key = calib_key(ops)
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
    parser.add_argument("--calibrate", action="store_true", help="Measure and write calibration constants (signal intervals) to file")
    parser.add_argument("--no-cuts", action="store_true", help="Dont cut anything (overrides default behavior)")
    parser.add_argument("--stop-after-hits", action="store_true", help="Stop the analysis after processing hits")
    parser.add_argument("--stop-after-mds", action="store_true", help="Stop the analysis after processing mini-doublets")
    parser.add_argument("--calib-dir", type=str, default=guess_calib_dir(), help="Directory of calibration constants")
    parser.add_argument("--layers", nargs="+", type=str, default=preset, help="List of layers to consider")
    parser.add_argument("--sim", action="store_true", help="Use sim hits in the analysis")
    parser.add_argument("--digi", action="store_true", help="Use digi hits in the analysis")
    parser.add_argument("--plot", action="store_true", help="Include plots in the analysis")
    parser.add_argument("--cut-mds", action="store_true", help="Cut MDs based on MD_DZ_CUT and MD_DR_CUT")
    parser.add_argument("--cut-t2s", action="store_true", help="Cut T2s based on [[ something ]]")
    parser.add_argument("--cut-t4s", action="store_true", help="Cut T4s based on [[ something ]]")
    parser.add_argument("--cut-t8s", action="store_true", help="Cut T8s based on [[ something ]]")
    parser.add_argument("--read-mcps", type=str, help="Read mcps from pickle file")
    parser.add_argument("--write-mcps", type=str, help="Write mcps to pickle file")
    parser.add_argument("--read-hits", type=str, help="Read hits from pickle file")
    parser.add_argument("--write-hits", type=str, help="Write hits to pickle file")
    parser.add_argument("--read-mds", type=str, help="Read mini-doublets from pickle file")
    parser.add_argument("--write-mds", type=str, help="Write mini-doublets to pickle file")
    parser.add_argument("--fast-mds", action="store_true", help="Use fast binned merge for mini-doublets")
    parser.add_argument("--read-t2s", type=str, help="Read T2s from pickle file")
    parser.add_argument("--write-t2s", type=str, help="Write T2s to pickle file")
    parser.add_argument("--read-t4s", type=str, help="Read T4s from pickle file")
    parser.add_argument("--write-t4s", type=str, help="Write T4s to pickle file")
    parser.add_argument("--read-t8s", type=str, help="Read T8s from pickle file")
    parser.add_argument("--write-t8s", type=str, help="Write T8s to pickle file")
    parser.add_argument("--geo", type=str, help="Version of geometry to use for cuts (e.g. v01, v04)", required=True)
    parser.add_argument("--smear", type=str, default="00um", help="Smear value to use for digi hits (e.g. 10um)")
    parser.add_argument("--dataset", type=str, help="Specify the dataset to use in the analysis")
    parser.add_argument("--cutflow", type=str, default="cutflow.ndjson", help="Path to output newline-delimited JSON for cutflows file")
    parser.add_argument("--debug", action="store_true", help="Print some debug information")
    parser.add_argument("--pdf", type=str, default="", help="Path to output PDF file")
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
