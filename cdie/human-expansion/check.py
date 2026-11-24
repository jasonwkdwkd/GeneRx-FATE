import math
from model import (Regulation, bottleneck, compare_coding_sequences, demo,
                   first_stop, nmd_heuristic, x_dosage)

d = demo()
family = d["paired_autosomal_and_sex_chromosome_alleles"]
assert family["autosome"]["a/a"] == .25
assert "x/Y" in family["XY_non_PAR"] and "x/X" in family["XX_non_PAR"]
assert x_dosage(1,0) == (1,0) and x_dosage(.5,1) == (1,1)
for h in (0,.37,1):
    dist = bottleneck(h,20)
    assert math.isclose(sum(p for _,p in dist),1)
    assert math.isclose(sum(x*p for x,p in dist),h)
assert d["mt_origins"]["heavy_start"] < d["mt_origins"]["light_start"] < d["mt_origins"]["complete"]
assert d["nuclear_fork_arrival"] == 6
assert first_stop("ATGTGATAA") == 6
assert first_stop("ATGTGATAA","vertebrate_mitochondrial") == 9
assert compare_coding_sequences("ATGCAATAA","ATGTAATAA",9)["stop_gain_without_frameshift"]
assert nmd_heuristic(10,100) and not nmd_heuristic(80,100) and not nmd_heuristic(80,None)
assert Regulation(cap_fraction=0).exportable_flux() == 0
assert Regulation(cleavage_poly_a_fraction=0).exportable_flux() == 0
assert Regulation(tata_present=False,tata_contribution=2).transcription_flux()==1
cases = d["nuclear_outputs"]
metric = lambda scenario,key: cases[scenario][0][key]
assert metric("early_stop_with_NMD","mRNA") < metric("reference","mRNA")
assert metric("late_stop_NMD_escape","mRNA") == metric("reference","mRNA")
assert metric("CFTR_like_gating_defect","delivered_protein") == metric("reference","delivered_protein")
assert metric("CFTR_like_gating_defect","functional_output") < metric("reference","functional_output")
assert metric("CFTR_like_trafficking_defect","delivered_protein") < metric("reference","delivered_protein")
mt = d["mt_translation_by_heteroplasmy"]
assert mt["0"]["usable_tRNA"] > mt["0.5"]["usable_tRNA"] > mt["1"]["usable_tRNA"]
same_mean = list(d["mt_tissue_same_mean_different_distribution"].values())
assert not math.isclose(*same_mean)
print("Passed synthetic inheritance, replication, RNA processing, translation and tissue aggregation checks.")
