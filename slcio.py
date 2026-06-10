import contextlib
import os
import sys
import numpy as np
import pandas as pd
import multiprocessing as mp
import logging
logger = logging.getLogger(__name__)

from constants import OUTSIDE_BOUNDS, INSIDE_BOUNDS, UNDEFINED_BOUNDS, BOUNDS
from constants import EPSILON, MCPARTICLE, PARTICLES_OF_INTEREST, SPEED_OF_LIGHT
from constants import MM_TO_CM, CM_TO_MM
from constants import XML
from constants import INNER_TRACKER_BARREL_COLLECTION, OUTER_TRACKER_BARREL_COLLECTION
from constants import INNER_TRACKER_BARREL_HITS, OUTER_TRACKER_BARREL_HITS
from constants import INNER_TRACKER_BARREL_RELATIONS, OUTER_TRACKER_BARREL_RELATIONS
from constants import BYTE_TO_MB, NO_MCP
from constants import MIN_COSTHETA, MIN_SIMHIT_PT_FRACTION, MAX_TIME
from constants import INNER_TRACKER_BARREL, OUTER_TRACKER_BARREL
from constants import NICKNAMES

_detector = None
_surfman = None
_maps = None

class HitMaker:

    def __init__(
            self,
            slcio_file_paths: list[str],
            load_geometry: bool,
            signal: bool,
            sim: bool,
            inner: bool,
            outer: bool,
            layers: list[int],
        ):
        self.slcio_file_paths = slcio_file_paths
        self.load_geometry = load_geometry
        self.signal = signal
        self.sim = sim
        self.inner = inner
        self.outer = outer
        self.layers = layers


    def convert(self) -> pd.DataFrame:
        mcps, simhits = self.convert_all_files()
        mcps = sort_mcps(mcps)
        simhits = sort_simhits(simhits)
        announce_inside_bounds(simhits)
        memory_usage = simhits.memory_usage(deep=True).sum() * BYTE_TO_MB
        logger.info(f"simhits.memory_usage: {memory_usage:.1f} MB")
        counts = simhits.groupby([
            "simhit_system",
            "simhit_layer",
        ]).size()
        for (system, layer), total in counts.items():
            logger.info(f"N(simhits) in system {system} layer {layer}: {total}")
        return mcps, simhits


    def convert_all_files(self) -> pd.DataFrame:
        logger.info(f"Converting {len(self.slcio_file_paths)} slcio files to a DataFrame ...")
        initializer = init_worker if self.load_geometry else init_dummy
        processes = min(mp.cpu_count(), len(self.slcio_file_paths))
        logger.info(f"Using {processes} processes for conversion ...")
        with mp.Pool(processes=processes, initializer=initializer) as pool:
            n_map = len(self.slcio_file_paths)
            file_numbers = list(range(n_map))
            load_geometry = [self.load_geometry]*n_map
            signal = [self.signal]*n_map
            sim = [self.sim]*n_map
            inner = [self.inner]*n_map
            outer = [self.outer]*n_map
            layers = [self.layers]*n_map
            results = pool.starmap(
                convert_one_file,
                zip(self.slcio_file_paths,
                    file_numbers,
                    load_geometry,
                    signal,
                    sim,
                    inner,
                    outer,
                    layers,
                )
            )
        logger.info("Merging DataFrames ...")
        return [
            pd.concat([mcps for (mcps, simhits) in results], ignore_index=True),
            pd.concat([simhits for (mcps, simhits) in results], ignore_index=True),
        ]


def init_dummy():
    pass


def init_worker():
    # Sorry for these global variables. They are needed for multiprocessing
    global _detector, _surfman, _maps
    import dd4hep, DDRec
    dd4hep.setPrintLevel(dd4hep.PrintLevel.WARNING)
    with silence_c_stdout_stderr():
        # Sorry for this context manager. dd4hep can be very noisy
        _detector = dd4hep.Detector.getInstance()
        _detector.fromCompact(XML)
        _surfman = DDRec.SurfaceManager(_detector)
        dets = {
            INNER_TRACKER_BARREL_COLLECTION: _detector.detector("InnerTrackerBarrel"),
            OUTER_TRACKER_BARREL_COLLECTION: _detector.detector("OuterTrackerBarrel"),
        }
        _maps = {name: _surfman.map(det.name()) for name, det in dets.items()}


def convert_one_file(
        slcio_file_path: str,
        file_number: int,
        load_geometry: bool,
        signal: bool,
        use_sim: bool,
        inner: bool,
        outer: bool,
        layers: list[int],
    ) -> tuple[pd.DataFrame, pd.DataFrame]:

    # import here to avoid:
    #  - unnecessary imports if not used
    #  - issues with multiprocessing
    with silence_c_stdout_stderr():
        import pyLCIO
    if load_geometry:
        import dd4hep
        import DDRec

    # check for file existence
    if not os.path.isfile(slcio_file_path):
        msg = f"File {slcio_file_path} does not exist"
        logger.error(msg)
        raise FileNotFoundError(msg)

    # open the SLCIO file
    logger.info(f"Processing file {slcio_file_path} ...")
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(slcio_file_path)

    # list for holding all hits
    mcps = []
    simhits = []

    # loop over all events in the slcio file
    for i_event, event in enumerate(reader):

        # if i_event > 0:
        #     break

        # inspect MCParticles
        mcparticles = list(event.getCollection(MCPARTICLE))
        mcp_px = [mcp.getMomentum()[0] for mcp in mcparticles]
        mcp_py = [mcp.getMomentum()[1] for mcp in mcparticles]
        mcp_pz = [mcp.getMomentum()[2] for mcp in mcparticles]
        mcp_m = [mcp.getMass() for mcp in mcparticles]
        mcp_q = [mcp.getCharge() for mcp in mcparticles]
        mcp_pdg = [mcp.getPDG() for mcp in mcparticles]
        mcp_vertex_x = [mcp.getVertex()[0] for mcp in mcparticles]
        mcp_vertex_y = [mcp.getVertex()[1] for mcp in mcparticles]
        mcp_vertex_z = [mcp.getVertex()[2] for mcp in mcparticles]
        mcp_endpoint_x = [mcp.getEndpoint()[0] for mcp in mcparticles]
        mcp_endpoint_y = [mcp.getEndpoint()[1] for mcp in mcparticles]
        mcp_endpoint_z = [mcp.getEndpoint()[2] for mcp in mcparticles]
        for i_mcp in range(len(mcparticles)):
            if abs(mcp_pdg[i_mcp]) not in PARTICLES_OF_INTEREST:
                continue
            mcps.append({
                'file': file_number,
                'i_event': i_event,
                'i_mcp': i_mcp,
                'mcp_px': mcp_px[i_mcp],
                'mcp_py': mcp_py[i_mcp],
                'mcp_pz': mcp_pz[i_mcp],
                'mcp_m': mcp_m[i_mcp],
                'mcp_q': mcp_q[i_mcp],
                'mcp_pdg': mcp_pdg[i_mcp],
                'mcp_vertex_x': mcp_vertex_x[i_mcp],
                'mcp_vertex_y': mcp_vertex_y[i_mcp],
                'mcp_vertex_z': mcp_vertex_z[i_mcp],
                'mcp_endpoint_x': mcp_endpoint_x[i_mcp],
                'mcp_endpoint_y': mcp_endpoint_y[i_mcp],
                'mcp_endpoint_z': mcp_endpoint_z[i_mcp],
            })

        # choose trackers
        collections = []
        if inner:
            if use_sim:
                collections.append(INNER_TRACKER_BARREL_COLLECTION)
            else:
                if signal:
                    collections.append(INNER_TRACKER_BARREL_RELATIONS)
                else:
                    collections.append(INNER_TRACKER_BARREL_HITS)
        if outer:
            if use_sim:
                collections.append(OUTER_TRACKER_BARREL_COLLECTION)
            else:
                if signal:
                    collections.append(OUTER_TRACKER_BARREL_RELATIONS)
                else:
                    collections.append(OUTER_TRACKER_BARREL_HITS)
        if len(collections) == 0:
            msg = f"No tracker collections selected for file {os.path.basename(slcio_file_path)}"
            logger.error(msg)
            raise ValueError(msg)

        # inspect tracking detectors
        for collection in collections:

            col = event.getCollection(collection)
            n_obj = len(col)

            for i_obj, obj in enumerate(col):

                # define which objects are available
                # be careful
                if use_sim:
                    simhit, hit = obj, None
                else:
                    if signal:
                        simhit, hit = obj.getTo(), obj.getFrom()
                    else:
                        simhit, hit = None, obj

                if i_obj > 0 and i_obj % 1000000 == 0:
                    logger.info(f"Processing file {os.path.basename(slcio_file_path)} "
                                f"event {i_event} collection {collection} "
                                f"hit {i_obj}/{n_obj} ...")

                # consider a particular set of layers
                cellid0 = simhit.getCellID0() if use_sim else hit.getCellID0()
                layer = np.right_shift(cellid0, 7) & 0b11_1111
                if layer not in layers:
                    continue
                # module 0 only, sensor 20 only?
                # if (np.right_shift(hit.getCellID0(), 13) & 0b111_1111_1111) != 0:
                #     continue
                # if (np.right_shift(hit.getCellID0(), 24) & 0b1111_1111) != 20:
                #     continue

                # associated MCParticle
                if signal or use_sim:
                    mcp = simhit.getMCParticle()
                    i_mcp = mcparticles.index(mcp) if mcp in mcparticles else NO_MCP
                else:
                    i_mcp = NO_MCP

                # more simhit or hit attributes
                position = simhit.getPosition() if use_sim else hit.getPosition()
                time = simhit.getTime() if use_sim else hit.getTime()
                energy = simhit.getEDep() if use_sim else hit.getEDep()
                momentum = simhit.getMomentum() if (use_sim or signal) else [0, 0, 0]
                pathlength = simhit.getPathLength() if (use_sim or signal) else 0
                correction = (np.sqrt(position[0]**2 + position[1]**2 + position[2]**2) / SPEED_OF_LIGHT) if use_sim else 0.0

                # hit/surface relations
                if load_geometry:
                    surf = _maps[collection].find(simhit.getCellID0()).second
                    pos = dd4hep.rec.Vector3D(position[0] * MM_TO_CM,
                                              position[1] * MM_TO_CM,
                                              position[2] * MM_TO_CM)
                    inside_bounds = INSIDE_BOUNDS if surf.insideBounds(pos) else OUTSIDE_BOUNDS
                    distance = surf.distance(pos) * CM_TO_MM
                else:
                    inside_bounds = UNDEFINED_BOUNDS
                    distance = -1

                # ignore hits outside bounds
                if inside_bounds == OUTSIDE_BOUNDS:
                    continue

                # record the hit info
                simhits.append({
                    'file': file_number,
                    'i_event': i_event,
                    'i_mcp': i_mcp,
                    'simhit_x': position[0],
                    'simhit_y': position[1],
                    'simhit_z': position[2],
                    'simhit_cellid0': cellid0,
                    'simhit_inside_bounds': inside_bounds,
                    'simhit_t_corrected': time - correction,
                })
                if signal:
                    mcp_ok = i_mcp != NO_MCP
                    simhits[-1].update({
                        'simhit_px': momentum[0],
                        'simhit_py': momentum[1],
                        'simhit_pz': momentum[2],
                        'simhit_pathlength': pathlength,
                        'simhit_distance': distance,
                        'simhit_t': time,
                        'simhit_e': energy,
                        'mcp_px': mcp_px[i_mcp] if mcp_ok else 0,
                        'mcp_py': mcp_py[i_mcp] if mcp_ok else 0,
                        'mcp_pz': mcp_pz[i_mcp] if mcp_ok else 0,
                        "mcp_pdg": mcp_pdg[i_mcp] if mcp_ok else 0,
                        "mcp_q": mcp_q[i_mcp] if mcp_ok else 0,
                        "mcp_vertex_x": mcp_vertex_x[i_mcp] if mcp_ok else 0,
                        "mcp_vertex_y": mcp_vertex_y[i_mcp] if mcp_ok else 0,
                        "mcp_vertex_z": mcp_vertex_z[i_mcp] if mcp_ok else 0,
                        "mcp_endpoint_x": mcp_endpoint_x[i_mcp] if mcp_ok else 0,
                        "mcp_endpoint_y": mcp_endpoint_y[i_mcp] if mcp_ok else 0,
                        "mcp_endpoint_z": mcp_endpoint_z[i_mcp] if mcp_ok else 0,
                    })

    # Close
    reader.close()

    # Convert the list of hits to a pandas DataFrame and postprocess
    logger.info("Creating DataFrames ...")
    mcps = pd.DataFrame(mcps)
    simhits = pd.DataFrame(simhits)

    # sanity check
    if len(mcps) == 0:
        msg = f"No MCParticles found in file {os.path.basename(slcio_file_path)}"
        logger.error(msg)
        raise RuntimeError(msg)
    if len(simhits) == 0:
        msg = f"No simhits found in file {os.path.basename(slcio_file_path)}"
        logger.error(msg)
        raise RuntimeError(msg)

    # And postprocess
    logger.info("Postprocessing DataFrames ...")
    mcps = postprocess_mcps(mcps)
    if signal:
        simhits = postprocess_mcps(simhits)
    simhits = postprocess_simhits(simhits, signal)

    # Bonus features: define if a mcp is "detectable" or not
    if signal:
        mcps = add_detectable_columns(mcps, simhits)

    return mcps, simhits


def postprocess_mcps(df: pd.DataFrame) -> pd.DataFrame:
    df["mcp_p"] = np.sqrt(df["mcp_px"]**2 + df["mcp_py"]**2 + df["mcp_pz"]**2)
    df["mcp_pt"] = np.sqrt(df["mcp_px"]**2 + df["mcp_py"]**2)
    df["mcp_theta"] = np.arctan2(df["mcp_pt"], df["mcp_pz"])
    df["mcp_eta"] = -np.log(np.tan(df["mcp_theta"] / 2))
    df["mcp_phi"] = np.arctan2(df["mcp_py"], df["mcp_px"])
    df["mcp_qoverpt"] = df["mcp_q"] / df["mcp_pt"]
    df["mcp_vertex_r"] = np.sqrt(df["mcp_vertex_x"]**2 + df["mcp_vertex_y"]**2)
    df["mcp_endpoint_r"] = np.sqrt(df["mcp_endpoint_x"]**2 + df["mcp_endpoint_y"]**2)

    # remove redundant columns
    df.drop(columns=[
        "mcp_px",
        "mcp_py",
        "mcp_theta",
        "mcp_vertex_x",
        "mcp_vertex_y",
        "mcp_endpoint_x",
        "mcp_endpoint_y",
    ], inplace=True)

    # downcast to save memory
    df["file"] = df["file"].astype(np.uint32)
    df["i_event"] = df["i_event"].astype(np.uint32)
    df["i_mcp"] = df["i_mcp"].astype(np.uint32)

    # sort columns alphabetically
    return df[sorted(df.columns)]


def postprocess_simhits(df: pd.DataFrame, signal: bool) -> pd.DataFrame:
    logger.info(f"Postprocessing DataFrame, signal={signal} ...")
    df["simhit_r"] = np.sqrt(df["simhit_x"]**2 + df["simhit_y"]**2)
    df["simhit_system"] = np.right_shift(df["simhit_cellid0"], 0) & 0b1_1111
    df["simhit_side"] = np.right_shift(df["simhit_cellid0"], 5) & 0b11
    df["simhit_layer"] = np.right_shift(df["simhit_cellid0"], 7) & 0b11_1111
    df["simhit_module"] = np.right_shift(df["simhit_cellid0"], 13) & 0b111_1111_1111
    df["simhit_sensor"] = np.right_shift(df["simhit_cellid0"], 24) & 0b1111_1111
    df["simhit_layer_div_2"] = df["simhit_layer"] // 2
    df["simhit_layer_mod_2"] = df["simhit_layer"] % 2
    # df["simhit_theta"] = np.maximum(np.arctan2(df["simhit_r"], df["simhit_z"]), EPSILON)
    # df["simhit_eta"] = -np.log(np.tan(df["simhit_theta"] / 2))
    # df["simhit_phi"] = np.arctan2(df["simhit_y"], df["simhit_x"])
    # df["simhit_pt"] = np.sqrt(df["simhit_px"]**2 + df["simhit_py"]**2)
    # df["simhit_R"] = np.sqrt(df["simhit_x"]**2 + df["simhit_y"]**2 + df["simhit_z"]**2)
    # df["simhit_t_corrected"] = df["simhit_t"] - (df["simhit_R"] / SPEED_OF_LIGHT)
    if signal:
        df["simhit_R"] = np.sqrt(df["simhit_x"]**2 + df["simhit_y"]**2 + df["simhit_z"]**2)
        df["simhit_p"] = np.sqrt(df["simhit_px"]**2 + df["simhit_py"]**2 + df["simhit_pz"]**2)
        df["simhit_costheta"] = (df["simhit_x"] * df["simhit_px"] +
                                 df["simhit_y"] * df["simhit_py"] +
                                 df["simhit_z"] * df["simhit_pz"]) / (df["simhit_R"] * df["simhit_p"])
        df["simhit_first_exit"] = (
            (df["simhit_t_corrected"] < MAX_TIME) &
            (df["simhit_costheta"] > MIN_COSTHETA) &
            (df["simhit_p"] / df["mcp_p"] > MIN_SIMHIT_PT_FRACTION)
        )

    # remove unused columns
    if signal:
        drop_cols = [
            "simhit_px",
            "simhit_py",
            "simhit_pz",
            "simhit_R",
        ]
        df.drop(columns=drop_cols, inplace=True)

    # downcast to save memory
    df["file"] = df["file"].astype(np.uint32)
    df["i_event"] = df["i_event"].astype(np.uint32)
    df["i_mcp"] = df["i_mcp"].astype(np.uint32)
    df["simhit_cellid0"] = df["simhit_cellid0"].astype(np.uint32)
    df["simhit_inside_bounds"] = df["simhit_inside_bounds"].astype(np.uint8)
    df["simhit_side"] = df["simhit_side"].astype(np.uint8)
    df["simhit_system"] = df["simhit_system"].astype(np.uint8)
    df["simhit_layer"] = df["simhit_layer"].astype(np.uint8)
    df["simhit_layer_div_2"] = df["simhit_layer_div_2"].astype(np.uint8)
    df["simhit_layer_mod_2"] = df["simhit_layer_mod_2"].astype(np.uint8)
    df["simhit_module"] = df["simhit_module"].astype(np.uint16)
    df["simhit_sensor"] = df["simhit_sensor"].astype(np.uint16)
    if signal:
        df["simhit_first_exit"] = df["simhit_first_exit"].astype(bool)

    # sort columns alphabetically
    return df[sorted(df.columns)]


def sort_mcps(df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Sorting DataFrame ...")
    columns = [
        "file",
        "i_event",
        "i_mcp",
    ]
    return df.sort_values(by=columns).reset_index(drop=True)


def sort_simhits(df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Sorting DataFrame ...")
    columns = [
        "file",
        "i_event",
        "i_mcp",
        "simhit_system",
        "simhit_layer",
        "simhit_module",
        "simhit_sensor",
    ]
    return df.sort_values(by=columns).reset_index(drop=True)


def add_detectable_columns(mcps: pd.DataFrame, simhits: pd.DataFrame) -> pd.DataFrame:

    GROUP_COLS = ["file", "i_event", "i_mcp"]
    MATCH_COLS = GROUP_COLS + ["simhit_module", "simhit_sensor"]
    SYSTEMS = [INNER_TRACKER_BARREL, OUTER_TRACKER_BARREL]
    LAYER_PAIRS = [(0, 1), (2, 3), (4, 5), (6, 7)]

    def get_detectable_mcps(
        simhits_system: pd.DataFrame,
        layer_lower: int,
        layer_upper: int,
    ) -> pd.DataFrame:
        """
        Returns a DataFrame of GROUP_COLS rows for every MCP that has at least
        one (module, sensor) pair with hits in both layer_lower and layer_upper.
        """
        lo = simhits_system[simhits_system["simhit_layer"] == layer_lower][MATCH_COLS]
        hi = simhits_system[simhits_system["simhit_layer"] == layer_upper][MATCH_COLS]
        return lo.merge(hi, on=MATCH_COLS)[GROUP_COLS].drop_duplicates()

    for system in SYSTEMS:

        nickname = NICKNAMES[system]
        sysdf = simhits[ simhits["simhit_first_exit"] & (simhits["simhit_system"] == system) ]

        for lo, hi in LAYER_PAIRS:

            column = f"mcp_detectable_{nickname}_{lo}{hi}"

            # get mcps that are detectable in each double layer
            detectable_mcps = get_detectable_mcps(sysdf, lo, hi).assign(**{column: True})

            # merge detectable mcps with the main mcp df
            mcps = mcps.merge(detectable_mcps, on=GROUP_COLS, how="left")

            # fill nan with False (not detectable)
            mcps[column] = mcps[column].fillna(False).astype(bool)

        # combine all double layers to get overall detectable status
        mcps[f"mcp_detectable_{nickname}"] = np.logical_and.reduce([
            mcps[f"mcp_detectable_{nickname}_{lo}{hi}"] for (lo, hi) in LAYER_PAIRS
        ])

    return mcps


def announce_inside_bounds(df: pd.DataFrame):
    for bounds in [OUTSIDE_BOUNDS, INSIDE_BOUNDS, UNDEFINED_BOUNDS]:
        n_bounds = len(df[df["simhit_inside_bounds"] == bounds])
        logger.info(f"N(simhits) with bounds == {BOUNDS[bounds]}: {n_bounds}")


@contextlib.contextmanager
def silence_c_stdout_stderr():
    # Flush Python buffers
    sys.stdout.flush()
    sys.stderr.flush()

    # Save original FDs
    old_stdout_fd = os.dup(1)
    old_stderr_fd = os.dup(2)

    # Redirect to /dev/null
    with open(os.devnull, "w") as devnull:
        os.dup2(devnull.fileno(), 1)
        os.dup2(devnull.fileno(), 2)
        try:
            yield
        finally:
            # Restore original FDs
            os.dup2(old_stdout_fd, 1)
            os.dup2(old_stderr_fd, 2)
            os.close(old_stdout_fd)
            os.close(old_stderr_fd)

