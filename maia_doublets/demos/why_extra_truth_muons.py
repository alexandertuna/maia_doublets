"""
I generated events with 10 muons per event and saved MCParticle info in a dataframe.
10 muons per event, 1000 events per file, 10 files gives 100000 muons.
But 100022 muon MCParticles are observed in the dataframe!
This explains whats happening.
"""

import numpy as np
import pandas as pd

FNAME = "/ceph/users/atuna/work/maia/maia_doublets/output/v06_signal_digi_10um/mcps.pkl"
MUON = 13
NEUTRINO = 14
PION = 211

def main():
    df = pd.read_pickle(FNAME)
    explain(df)
    df = df[np.abs(df["mcp_pdg"]) == MUON]


def explain(df):
    explain_pdg(df)
    explain_vertices(df)


def explain_pdg(df):
    for (particle, name) in [(False, "all"),
                             (MUON, "muon"),
                             (NEUTRINO, "muon neutrino"),
                             (PION, "pion"),
                             ]:
        count = (np.abs(df["mcp_pdg"]) == particle).sum() if particle else len(df)
        print(f"N(MCP): {count:>6} {name}")
    print()


def explain_vertices(df):
    muons = df[np.abs(df["mcp_pdg"]) == MUON]
    vertices = muons[["mcp_vertex_r", "mcp_vertex_z"]].drop_duplicates()
    print(f"N(unique muon vertices): {len(vertices)}")
    print(vertices)
    print()
    displaced = muons[(muons["mcp_vertex_r"] > 1)]
    cols = ["file", "i_event", "i_mcp", "mcp_vertex_r", "mcp_vertex_z"]
    print(f"N(displaced muon): {len(displaced)}")
    print(displaced[cols])
    print()

if __name__ == "__main__":
    main()
