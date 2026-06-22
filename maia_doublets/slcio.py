import contextlib
import os
import sys
import numpy as np
import pandas as pd
import uproot
import multiprocessing as mp
import logging
logger = logging.getLogger(__name__)

from maia_doublets.constants import OUTSIDE_BOUNDS, INSIDE_BOUNDS, UNDEFINED_BOUNDS, BOUNDS
from maia_doublets.constants import EPSILON, MCPARTICLE, PARTICLES_OF_INTEREST, SPEED_OF_LIGHT
from maia_doublets.constants import MM_TO_CM, CM_TO_MM
from maia_doublets.constants import XML
from maia_doublets.constants import INNER_TRACKER_BARREL_COLLECTION, OUTER_TRACKER_BARREL_COLLECTION
from maia_doublets.constants import INNER_TRACKER_BARREL_HITS, OUTER_TRACKER_BARREL_HITS
from maia_doublets.constants import INNER_TRACKER_BARREL_RELATIONS, OUTER_TRACKER_BARREL_RELATIONS
from maia_doublets.constants import BYTE_TO_MB, NO_MCP, PODIO_NO_MCP
from maia_doublets.constants import MIN_COSTHETA, MIN_SIMHIT_PT_FRACTION, MAX_TIME
from maia_doublets.constants import INNER_TRACKER_BARREL, OUTER_TRACKER_BARREL
from maia_doublets.constants import NICKNAMES, LAYER_OFFSET
from maia_doublets.constants import MUON, ONE_POINT_FIVE_GEV, BARREL_TRACKER_MAX_ETA, ZERO_POINT_ZERO_ONE_MM

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
            layers: dict[int, set[int]],
        ):
        self.slcio_file_paths = slcio_file_paths
        self.load_geometry = load_geometry
        self.signal = signal
        self.sim = sim
        self.layers = layers


    def convert_debug(self):
        # tmp debug
        logger.info("Opening pkl ...")
        # lcio_mcps = pd.read_pickle("/ceph/users/atuna/work/maia/maia_doublets/run/mcps.pkl")
        # lcio_hits = pd.read_pickle("/ceph/users/atuna/work/maia/maia_doublets/run/simhits.pkl")
        # lcio_mcps = pd.read_pickle("../data/signal/mcps.pkl")
        # lcio_hits = pd.read_pickle("../data/signal/simhits.pkl")
        # lcio_mcps = pd.read_pickle("../data/signal_10um/mcps.pkl")
        # lcio_hits = pd.read_pickle("../data/signal_10um/simhits.pkl")
        lcio_mcps = pd.read_pickle("mcps.pkl")
        lcio_hits = pd.read_pickle("simhits.pkl")
        logger.info("Converting ROOT files ...")
        fnames = [pat + ".root" for pat in self.slcio_file_paths]
        results = [
            convert_one_root_file(fnames[i], i, self.load_geometry, self.signal, self.sim, self.layers)
            for i in range(len(fnames))
        ]
        root_mcps = pd.concat([mcps for (mcps, simhits) in results], ignore_index=True)
        root_hits = pd.concat([simhits for (mcps, simhits) in results], ignore_index=True)
        logger.info("ROOT files converted ...")
        # root_mcps, root_hits = convert_one_root_file(fnames[0], 0, self.load_geometry, self.signal, self.sim, self.layers)

        print("root_hits\n", root_hits)
        print("lcio_hits\n", lcio_hits)
        print("x"*40)
        print("root_hits.equals(lcio_hits):", root_hits.equals(lcio_hits))
        print("root_mcps.equals(lcio_mcps):", root_mcps.equals(lcio_mcps))
        print("x"*40)
        for col in root_hits.columns:
            if not col in lcio_hits.columns:
                print(f"Column {col} is not in lcio_hits")
        for col in lcio_hits.columns:
            if not col in root_hits.columns:
                print(f"Column {col} is not in root_hits")
        print("*"*40)
        from pandas.testing import assert_frame_equal
        try:
            assert_frame_equal(root_hits, lcio_hits, check_dtype=False, check_column_type=False)
            logger.info("Hits dataframes are equal according to assert_frame_equal")
        except AssertionError as e:
            logger.error(f"Hits dataframes are not equal: {e}")
        try:
            assert_frame_equal(root_mcps, lcio_mcps, check_dtype=False, check_column_type=False)
            logger.info("MCP dataframes are equal according to assert_frame_equal")
        except AssertionError as e:
            logger.error(f"MCP dataframes are not equal: {e}")
        print("x"*40)

        # check for any non-equal columns in hits
        for col in root_hits.columns:
            if not root_hits[col].equals(lcio_hits[col]):
                if np.allclose(root_hits[col], lcio_hits[col]):
                    print(f"Column {col} is not exactly equal but close. dtype: root {root_hits[col].dtype}, lcio {lcio_hits[col].dtype}")
                    continue
                print(f"Column {col} is not equal or close. dtype: root {root_hits[col].dtype}, lcio {lcio_hits[col].dtype}")
                break
        else:
            logger.info(f"All columns in hits are either equal or close")
        print("x"*40)

        # check for any non-equal columns in mcps
        for col in root_mcps.columns:
            if not root_mcps[col].equals(lcio_mcps[col]):
                if np.allclose(root_mcps[col], lcio_mcps[col]):
                    print(f"Column {col} is not exactly equal but close. dtype: root {root_mcps[col].dtype}, lcio {lcio_mcps[col].dtype}")
                    continue
                print(f"Column {col} is not equal or close. dtype: root {root_mcps[col].dtype}, lcio {lcio_mcps[col].dtype}")
                break
        else:
            logger.info(f"All columns in mcps are either equal or close")
        print("x"*40)

        import sys
        sys.exit("Testing purposes")


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
            layers = [self.layers]*n_map
            results = pool.starmap(
                convert_one_file,
                zip(self.slcio_file_paths,
                    file_numbers,
                    load_geometry,
                    signal,
                    sim,
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


def convert_one_root_file(
        root_file_path: str,
        file_number: int,
        load_geometry: bool,
        signal: bool,
        use_sim: bool,
        layers: dict[int, set[int]],
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse a ROOT file for MCParticles and simhits, returning two DataFrames.
    This is an alternative to convert_one_file that uses uproot instead of pyLCIO.
    """
    rf = uproot.open(root_file_path)
    evs = rf["events"]

    logger.info(f"Converting ROOT file {root_file_path} to mcps ...")
    mcps = convert_one_root_file_to_mcps(evs, file_number)

    logger.info(f"Converting ROOT file {root_file_path} to hits ...")
    hits = convert_one_root_file_to_hits(evs, file_number, signal, use_sim, layers)

    # connect hits and mcps
    if signal:
        hits = merge_mcp_info_into_hits(hits, mcps)

    # keep interesting mcps only
    mask = np.abs(mcps["mcp_pdg"]).isin(PARTICLES_OF_INTEREST)
    mcps = mcps[mask].reset_index(drop=True)

    # post-process
    logger.info(f"Post-processing DataFrames for ROOT file {root_file_path} ...")
    mcps = postprocess_mcps(mcps)
    mcps = sort_mcps(mcps)
    if signal:
        hits = postprocess_mcps(hits)
    hits = postprocess_simhits(hits, signal)
    hits = sort_simhits(hits)

    # Bonus features: define if a mcp is "detectable" or not
    if signal:
        mcps = add_detectable_columns(mcps, hits)

    return mcps, hits


def convert_one_root_file_to_hits(evs: uproot.TTree,
                                  file_number: int,
                                  signal: bool,
                                  use_sim: bool,
                                  layers: dict[int, set[int]],
                                  ) -> pd.DataFrame:

    sim_cols = [
        INNER_TRACKER_BARREL_COLLECTION,
        OUTER_TRACKER_BARREL_COLLECTION,
    ]
    digi_cols = [
        INNER_TRACKER_BARREL_HITS,
        OUTER_TRACKER_BARREL_HITS,
    ]
    rel_cols = [
        INNER_TRACKER_BARREL_RELATIONS,
        OUTER_TRACKER_BARREL_RELATIONS
    ]

    datas = {}
    for (sim_col, digi_col, rel_col) in zip(sim_cols, digi_cols, rel_cols):
        datas[sim_col] = convert_one_root_file_to_hits_per_system(evs=evs,
                                                                  file_number=file_number,
                                                                  signal=signal,
                                                                  use_sim=use_sim,
                                                                  layers=layers,
                                                                  sim_col=sim_col,
                                                                  digi_col=digi_col,
                                                                  rel_col=rel_col,
                                                                  )

    hits = pd.concat(datas.values(), ignore_index=True)
    return hits


def convert_one_root_file_to_hits_per_system(
        evs: uproot.TTree,
        file_number: int,
        signal: bool,
        use_sim: bool,
        layers: dict[int, set[int]],
        sim_col: str,
        digi_col: str,
        rel_col: str
        ) -> pd.DataFrame:
    """
    If use_sim:
        We need the sim collections
        We need the link from sim hit to MCParticle
    Else:
        If signal:
            We need the digi collections
            We need the sim collections
            We need the relation from digi hit to sim hit
            We need the link from sim hit to MCParticle
        Else:
            We need the digi collections
    """
    logger.info(f"Converting ROOT file to hits for {sim_col} ...")

    # names = {}
    # ancil = {}

    # todo: add comment
    sim_basics = {
        f"{sim_col}.position.x": "simhit_x",
        f"{sim_col}.position.y": "simhit_y",
        f"{sim_col}.position.z": "simhit_z",
        f"{sim_col}.time": "simhit_t",
        f"{sim_col}.cellID": "simhit_cellid0",
    }

    digi_basics = {
        f"{digi_col}.position.x": "simhit_x",
        f"{digi_col}.position.y": "simhit_y",
        f"{digi_col}.position.z": "simhit_z",
        f"{digi_col}.time": "simhit_t",
        f"{digi_col}.cellID": "simhit_cellid0",
    }

    sim_extra = {
        f"{sim_col}.momentum.x": "simhit_px",
        f"{sim_col}.momentum.y": "simhit_py",
        f"{sim_col}.momentum.z": "simhit_pz",
        f"{sim_col}.eDep": "simhit_e",
        f"{sim_col}.pathLength": "simhit_pathlength",
        f"_{sim_col}_particle.index": "i_mcp",
    }

    digi_extra = {
        f"_{rel_col}_to.index": "i_to",
        f"_{rel_col}_from.index": "i_from",
    }

    # # basic info
    # if use_sim:
    #     names |= {
    #         f"{sim_col}.position.x": "simhit_x",
    #         f"{sim_col}.position.y": "simhit_y",
    #         f"{sim_col}.position.z": "simhit_z",
    #         f"{sim_col}.time": "simhit_t",
    #         f"{sim_col}.cellID": "simhit_cellid0",
    #     }
    # else:
    #     names |= {
    #         f"{digi_col}.position.x": "simhit_x",
    #         f"{digi_col}.position.y": "simhit_y",
    #         f"{digi_col}.position.z": "simhit_z",
    #         f"{digi_col}.time": "simhit_t",
    #         f"{digi_col}.cellID": "simhit_cellid0",
    #     }
    # # more sim hit info
    # if signal or use_sim:
    #     sim_features = {
    #         f"{sim_col}.momentum.x": "simhit_px",
    #         f"{sim_col}.momentum.y": "simhit_py",
    #         f"{sim_col}.momentum.z": "simhit_pz",
    #         f"{sim_col}.eDep": "simhit_e",
    #         f"{sim_col}.pathLength": "simhit_pathlength",
    #         f"_{sim_col}_particle.index": "i_mcp",
    #     }
    #     if use_sim:
    #         names |= sim_features
    #     elif signal:
    #         ancil |= sim_features

    # the relevant data in each circumstance
    sim_columns = {}
    digi_columns = {}
    if use_sim or signal:
        sim_columns |= sim_basics
        if signal:
            sim_columns |= sim_extra
    if not use_sim:
        digi_columns |= digi_basics
        if signal:
            digi_columns |= digi_extra

    # fetch the data
    data = {}
    hits = {}
    for (col, columns) in [
        (sim_col, sim_columns),
        (digi_col, digi_columns),
    ]:

        if not columns:
            continue

        # fetch data
        data[col] = evs.arrays(list(columns.keys()), library="np")

        # metadata: event number
        n_events = len(data[col][f"{col}.cellID"])
        i_events = np.repeat(np.arange(n_events), [len(arr) for arr in data[col][f"{col}.cellID"]])
        n_rows = len(i_events)

        # check
        if n_rows == 0:
            msg = f"No hits found in {col}"
            logger.error(msg)
            raise RuntimeError(msg)

        # flatten the data
        for key, value in columns.items():
            data[col][value] = np.concatenate(data[col].pop(key))

        # engineer a few columns
        data[col]["file"] = np.array([file_number]*n_rows)
        data[col]["i_event"] = np.array(i_events)
        data[col]["simhit_inside_bounds"] = np.array([UNDEFINED_BOUNDS]*n_rows)
        if col == sim_col:
            correction = (np.sqrt(data[col]["simhit_x"]**2 + \
                                  data[col]["simhit_y"]**2 + \
                                  data[col]["simhit_z"]**2) / SPEED_OF_LIGHT)
            data[col]["simhit_t_corrected"] = data[col]["simhit_t"] - correction
        if signal:
            data[col]["simhit_distance"] = np.array([-1]*n_rows)
        else:
            data[col]["i_mcp"] = np.array([NO_MCP]*n_rows)

        # convert to DataFrame
        hits[col] = pd.DataFrame(data[col])

        # adjust the "no MCParticle" value for i_mcp
        hits[col]["i_mcp"] = hits[col]["i_mcp"].replace(PODIO_NO_MCP, NO_MCP)

        # sanity check
        hit_system = np.right_shift(hits[col]["simhit_cellid0"], 0) & 0b1_1111
        if hit_system.nunique() != 1:
            msg = f"Expected one system, found {hit_system.unique()}"
            logger.error(msg)
            raise RuntimeError(msg)

        # filter by layer
        hit_system = hit_system.iloc[0]
        hit_layer = np.right_shift(hits[col]["simhit_cellid0"], 7) & 0b11_1111
        mask = hit_layer.isin(layers[hit_system])
        hits[col] = hits[col][mask]


    # if sim_columns:
    #     sim_data = evs.arrays(list(sim_columns.keys()), library="np")
    # if digi_columns:
    #     digi_data = evs.arrays(list(digi_columns.keys()), library="np")

    # # metadata: event number
    # meta_col = sim_col if use_sim else digi_col
    # i_events = []
    # for i_event, branch in enumerate(sim_data[f"{meta_col}.cellID"]):
    #     i_events.extend([i_event]*len(branch))

    # # count them up
    # n_rows = len(i_events)
    # if n_rows == 0:
    #     msg = f"No hits found in {sim_col}"
    #     logger.error(msg)
    #     raise RuntimeError(msg)

    # # numpify the data
    # for key, value in sim_columns.items():
    #     sim_data[value] = np.concatenate(sim_data.pop(key))

    # # engineer a few columns
    # sim_data["file"] = np.array([file_number]*n_rows)
    # sim_data["i_event"] = np.array(i_events)
    # correction = (np.sqrt(sim_data["simhit_x"]**2 + \
    #                       sim_data["simhit_y"]**2 + \
    #                       sim_data["simhit_z"]**2) / SPEED_OF_LIGHT) if use_sim else 0.0
    # sim_data["simhit_t_corrected"] = sim_data["simhit_t"] - correction
    # sim_data["simhit_inside_bounds"] = np.array([UNDEFINED_BOUNDS]*n_rows)
    # if signal:
    #     sim_data["simhit_distance"] = np.array([-1]*n_rows)
    # else:
    #     sim_data["i_mcp"] = np.array([NO_MCP]*n_rows)

    # engineer a tricky column:
    # navigate from digi hits to sim hits to mcp index
    if False and signal and not use_sim:
        # simhit, hit = getTo(), getFrom()
        key_to = f"_{rel_col}_to.index"
        key_from = f"_{rel_col}_from.index"
        key_mcp = f"_{sim_col}_particle.index"
        rels = evs.arrays([key_to, key_from, key_mcp], library="np")
        rel_to = np.concatenate(rels.pop(key_to))
        rel_from = np.concatenate(rels.pop(key_from))
        sim_mcp = np.concatenate(rels.pop(key_mcp))
        # subset = 1880
        # print("allclose", np.allclose(rel_to[subset-10:subset], rel_from[subset-10:subset]))
        print("rel_to.shape", rel_to.shape)
        print("rel_from.shape", rel_from.shape)
        print("sim_mcp.shape", sim_mcp.shape)
        print("sim_mcp[rel_to].shape", sim_mcp[rel_to].shape)
        i_mcp = np.ones_like(rel_to) * NO_MCP

        rel_df = pd.DataFrame({
            "i_event": data["i_event"],
            "i_digi":  np.concatenate(rel_from),
            "i_sim":   np.concatenate(rel_to),
        })

        for (r_to, r_from) in zip(rel_to, rel_from):
            i_mcp[r_from] = sim_mcp[r_to]
        print("i_mcp.shape", i_mcp.shape)
        print("i_mcp", i_mcp)
        # print("rel_to[:10]", rel_to[subset-10:subset])
        # print("rel_from[:10]", rel_from[subset-10:subset])


    # # convert to DataFrame
    # sim_hits = pd.DataFrame(sim_data)

    # # adjust the "no MCParticle" value for i_mcp
    # sim_hits["i_mcp"] = sim_hits["i_mcp"].replace(PODIO_NO_MCP, NO_MCP)

    # # sanity check
    # hit_system = np.right_shift(sim_hits["simhit_cellid0"], 0) & 0b1_1111
    # if hit_system.nunique() != 1:
    #     msg = f"Expected one system, found {hit_system.unique()}"
    #     logger.error(msg)
    #     raise RuntimeError(msg)

    # # filter by layer
    # hit_system = hit_system.iloc[0]
    # hit_layer = np.right_shift(sim_hits["simhit_cellid0"], 7) & 0b11_1111
    # mask = hit_layer.isin(layers[hit_system])
    # sim_hits = sim_hits[mask]

    if use_sim:
        return hits[sim_col]
        # return sim_hits
    else:
        raise NotImplementedError("Digi hits are not yet supported in this function")

    # return sim_hits


def convert_one_root_file_to_mcps(evs: uproot.TTree,
                                  file_number: int,
                                  ) -> pd.DataFrame:
    names = {
        "MCParticle.momentum.x": "mcp_px",
        "MCParticle.momentum.y": "mcp_py",
        "MCParticle.momentum.z": "mcp_pz",
        "MCParticle.mass": "mcp_m",
        "MCParticle.charge": "mcp_q",
        "MCParticle.PDG": "mcp_pdg",
        "MCParticle.vertex.x": "mcp_vertex_x",
        "MCParticle.vertex.y": "mcp_vertex_y",
        "MCParticle.vertex.z": "mcp_vertex_z",
        "MCParticle.endpoint.x": "mcp_endpoint_x",
        "MCParticle.endpoint.y": "mcp_endpoint_y",
        "MCParticle.endpoint.z": "mcp_endpoint_z",
    }
    data = evs.arrays(list(names.keys()), library="np")

    # metadata: event number and mcp index within event
    i_events, i_mcps = [], []
    for i_event, branch in enumerate(data["MCParticle.PDG"]):
        n_mcps = len(branch)
        i_events.extend([i_event]*n_mcps)
        i_mcps.extend(list(range(n_mcps)))

    # count them up
    n_rows = len(i_events)
    if n_rows == 0:
        msg = f"No MCPs found"
        logger.error(msg)
        raise RuntimeError(msg)

    # numpify the data
    for key, value in names.items():
        data[value] = np.concatenate(data.pop(key))
    data["file"] = np.array([file_number]*n_rows)
    data["i_event"] = np.array(i_events)
    data["i_mcp"] = np.array(i_mcps)

    # convert to DataFrame
    mcps = pd.DataFrame(data)

    # post-process
    return mcps


def merge_mcp_info_into_hits(
    hits: pd.DataFrame,
    mcps: pd.DataFrame,
) -> pd.DataFrame:

    on_cols = ["file", "i_event", "i_mcp"]
    mcp_cols = [
        "mcp_pdg", "mcp_q",
        "mcp_px", "mcp_py", "mcp_pz",
        "mcp_vertex_x", "mcp_vertex_y", "mcp_vertex_z",
        "mcp_endpoint_x", "mcp_endpoint_y", "mcp_endpoint_z",
        ]
    sub_cols = on_cols + mcp_cols

    hits = hits.merge(mcps[sub_cols], on=on_cols)
    hits[mcp_cols] = hits[mcp_cols].fillna(0.0)

    return hits


def convert_one_file(
        slcio_file_path: str,
        file_number: int,
        load_geometry: bool,
        signal: bool,
        use_sim: bool,
        layers: dict[int, set[int]],
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

    # collections to process
    if use_sim:
        collections = [
            # Sim hits
            INNER_TRACKER_BARREL_COLLECTION,
            OUTER_TRACKER_BARREL_COLLECTION,
        ]
    else:
        if signal:
            collections = [
                # Relations, which contain sim hits and digi hits
                INNER_TRACKER_BARREL_RELATIONS,
                OUTER_TRACKER_BARREL_RELATIONS,
            ]
        else:
            collections = [
                # Digi hits
                INNER_TRACKER_BARREL_HITS,
                OUTER_TRACKER_BARREL_HITS,
            ]

    if len(collections) == 0:
        msg = f"No tracker collections selected for file {os.path.basename(slcio_file_path)}"
        logger.error(msg)
        raise ValueError(msg)

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

                if i_obj > 0 and i_obj % 1e6 == 0:
                    logger.info(f"Processing file {os.path.basename(slcio_file_path)} "
                                f"event {i_event} collection {collection} "
                                f"hit {i_obj}/{n_obj} ...")

                # consider a particular set of layers
                cellid0 = simhit.getCellID0() if use_sim else hit.getCellID0()
                system = np.right_shift(cellid0, 0) & 0b1_1111
                layer = np.right_shift(cellid0, 7) & 0b11_1111
                if system not in layers:
                    continue
                if layer not in layers[system]:
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
    df["mcp_pdg"] = df["mcp_pdg"].astype(np.int32)
    df["mcp_q"] = df["mcp_q"].astype(np.float32)

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
    df["simhit_glayer"] = df["simhit_layer"] + LAYER_OFFSET[df["simhit_system"]]
    if signal:
        df["simhit_R"] = np.sqrt(df["simhit_x"]**2 + df["simhit_y"]**2 + df["simhit_z"]**2)
        df["simhit_p"] = np.sqrt(df["simhit_px"]**2 + df["simhit_py"]**2 + df["simhit_pz"]**2)
        df["simhit_costheta"] = (df["simhit_x"] * df["simhit_px"] +
                                 df["simhit_y"] * df["simhit_py"] +
                                 df["simhit_z"] * df["simhit_pz"]) / (df["simhit_R"] * df["simhit_p"])
        df["simhit_from_fiducial_mcp"] = (
            (df["i_mcp"] != NO_MCP) &
            (np.abs(df["mcp_pdg"]) == MUON) &
            (df["mcp_q"] != 0) &
            (df["mcp_pt"] > ONE_POINT_FIVE_GEV) &
            (np.abs(df["mcp_eta"]) < BARREL_TRACKER_MAX_ETA) &
            (df["mcp_vertex_r"] < ZERO_POINT_ZERO_ONE_MM) &
            (np.abs(df["mcp_vertex_z"]) < ZERO_POINT_ZERO_ONE_MM)
        )
        df["simhit_first_exit"] = (
            (df["simhit_t_corrected"] < MAX_TIME) &
            (df["simhit_costheta"] > MIN_COSTHETA) &
            (df["simhit_p"] / df["mcp_p"] > MIN_SIMHIT_PT_FRACTION)
        )

    # remove unused columns
    drop_cols = [
        # "simhit_cellid0",
    ]
    if "simhit_t" in df.columns:
        drop_cols.append("simhit_t")
    if signal:
        drop_cols += [
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
    df["simhit_inside_bounds"] = df["simhit_inside_bounds"].astype(np.uint8)
    df["simhit_side"] = df["simhit_side"].astype(np.uint8)
    df["simhit_system"] = df["simhit_system"].astype(np.uint8)
    df["simhit_layer"] = df["simhit_layer"].astype(np.uint8)
    df["simhit_glayer"] = df["simhit_glayer"].astype(np.uint8)
    df["simhit_layer_div_2"] = df["simhit_layer_div_2"].astype(np.uint8)
    df["simhit_layer_mod_2"] = df["simhit_layer_mod_2"].astype(np.uint8)
    df["simhit_module"] = df["simhit_module"].astype(np.uint16)
    df["simhit_sensor"] = df["simhit_sensor"].astype(np.uint16)
    df["simhit_x"] = df["simhit_x"].astype(np.float32)
    df["simhit_y"] = df["simhit_y"].astype(np.float32)
    df["simhit_z"] = df["simhit_z"].astype(np.float32)
    df["simhit_r"] = df["simhit_r"].astype(np.float32)
    df["simhit_t_corrected"] = df["simhit_t_corrected"].astype(np.float32)
    if signal:
        df["simhit_p"] = df["simhit_p"].astype(np.float32)
        df["simhit_e"] = df["simhit_e"].astype(np.float32)
        df["simhit_pathlength"] = df["simhit_pathlength"].astype(np.float32)
        df["simhit_first_exit"] = df["simhit_first_exit"].astype(bool)
        df["simhit_from_fiducial_mcp"] = df["simhit_from_fiducial_mcp"].astype(bool)

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

            # make sure column is a bool
            mcps[column] = mcps[column].astype(bool)

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

