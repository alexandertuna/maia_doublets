"""
slcio files can be converted to ROOT files like:
> lcio2edm4hep /ceph/users/atuna/work/maia/maia_noodling/samples/v01/neutrinoGun_n5_p15/neutrinoGun_digi_3.slcio neutrinoGun_digi_3.root
> lcio2edm4hep /ceph/users/atuna/work/maia/maia_noodling/samples/v01/muonGun_pT_2p0_2p1/muonGun_pT_2p0_2p1_sim_300.slcio muonGun_pT_2p0_2p1_sim_300.root
"""
import uproot

rf = uproot.open("neutrinoGun_digi_3.root")
# rf = uproot.open("muonGun_pT_2p0_2p1_sim_300.root")
ev = rf["events"]

for br in ev.branches:
    print(br)
print("*"*20)

for br in ev["MCParticle"]:
    print(br)
print("*"*20)

for br in ev["OuterTrackerBarrelCollection"]:
    print(br)
print("*"*20)

for br in ev["_OuterTrackerBarrelCollection_particle"]:
    print(br)
print("*"*20)

bnames = [
    "OuterTrackerBarrelCollection.position.x",
    "_OuterTrackerBarrelCollection_particle.index",
    "_OuterTrackerBarrelCollection_particle.collectionID",
    "MCParticle.PDG",
]
dic = ev.arrays(bnames, library="np")
pos_x = dic[bnames[0]]
print("len(pos_x)", len(pos_x))

pos_x_0 = pos_x[0]
print("len(pos_x_0)", len(pos_x_0))
print("*"*20)

p_index = dic[bnames[1]]
p_index_0 = p_index[0]
print("len(p_index_0)", len(p_index_0))
print("p_index_0", p_index_0)
print("*"*20)

col_id = dic[bnames[2]]
col_id_0 = col_id[0]
print("len(col_id_0)", len(col_id_0))
print("col_id_0", col_id_0)

pdg_id = dic[bnames[3]]
pdg_id_0 = pdg_id[0]
print("len(pdg_id_0)", len(pdg_id_0))
print("pdg_id_0", pdg_id_0)
