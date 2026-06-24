from glob import glob

signal_filepaths = {
    ("v01", "sim"): [
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_300.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_301.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_302.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_303.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_304.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_305.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_306.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_307.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_308.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_309.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_310.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_311.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_312.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_313.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_314.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_315.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_316.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_317.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_318.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_319.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_320.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_321.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_322.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_323.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_324.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_325.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_326.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_327.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_328.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_329.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_330.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_331.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_332.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_333.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_334.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_335.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_336.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_337.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_338.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_339.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_340.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_341.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_342.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_343.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_344.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_345.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_346.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_347.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_348.slcio",
        # "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_349.slcio",
    ],
    ("v05", "sim"): [
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v05/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_100.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v05/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_101.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v05/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_102.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v05/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_103.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v05/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_104.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v05/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_105.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v05/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_106.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v05/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_107.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v05/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_108.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v05/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_109.slcio",
    ],
    ("v01", "digi", "10um"): [
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/10um/muonGun_pT_2p0_2p1_digi_300.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/10um/muonGun_pT_2p0_2p1_digi_301.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/10um/muonGun_pT_2p0_2p1_digi_302.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/10um/muonGun_pT_2p0_2p1_digi_303.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/10um/muonGun_pT_2p0_2p1_digi_304.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/10um/muonGun_pT_2p0_2p1_digi_305.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/10um/muonGun_pT_2p0_2p1_digi_306.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/10um/muonGun_pT_2p0_2p1_digi_307.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/10um/muonGun_pT_2p0_2p1_digi_308.slcio",
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/10um/muonGun_pT_2p0_2p1_digi_309.slcio",
    ],
}

background100_filepaths = {
    ("v01", "sim"): [
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/neutrinoGun_n5_p15/neutrinoGun_digi_3.slcio",
    ],
    ("v01", "digi", "10um"): [
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/neutrinoGun/10um/neutrinoGun_digi_3.slcio",
    ],
}

background10_filepaths = {
    ("v01", "sim"): [
        # neutrinoGun 10% (-5, 15)
        "/ceph/users/atuna/work/maia/maia_noodling/samples/v01/neutrinoGun_n5_p15_0.10/neutrinoGun_digi_3.slcio",
    ]
}

def parse_filepaths(
    fnames: str | list[str],
) -> list[str]:
    names = []
    if isinstance(fnames, str):
        fnames = fnames.split(",")
    for fname in fnames:
        names.extend(sorted(glob(fname)))
    return names

def get_filepaths(
    geometry_version: str,
    signal: bool,
    background10: bool,
    background100: bool,
    sim: bool,
    digi: bool,
    smear: str,
) -> list[str]:
    if not sim and not digi:
        raise ValueError("Must specify either sim or digi")
    fpaths = []
    if sim:
        if signal:
            fpaths = signal_filepaths[(geometry_version, "sim")]
        elif background100:
            fpaths = background100_filepaths[(geometry_version, "sim")]
        elif background10:
            fpaths = background10_filepaths[(geometry_version, "sim")]
        else:
            raise ValueError("Must specify either signal or background filepaths.")
    elif digi:
        if signal:
            fpaths = signal_filepaths[(geometry_version, "digi", smear)]
        elif background100:
            fpaths = background100_filepaths[(geometry_version, "digi", smear)]
        elif background10:
            fpaths = background10_filepaths[(geometry_version, "digi", smear)]
        else:
            raise ValueError("Must specify either signal or background filepaths.")
    return parse_filepaths(fpaths)
