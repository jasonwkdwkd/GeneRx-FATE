from math import isclose

from compartments import demo

d = demo()
assert d["nuclear_CCAAT_loss"]["NF_Y_CCAAT_occupancy"] == 0
assert d["nuclear_CCAAT_loss"]["pre_mRNA_flux"] < d["nuclear_reference"]["pre_mRNA_flux"]
assert d["mt_replication_low_SSBP1"]["relative_timing"]["complete"] > d["mt_replication_reference"]["relative_timing"]["complete"]
assert d["mt_high_heteroplasmy"]["HSP_primary_output"] == d["mt_reference"]["HSP_primary_output"]
assert d["mt_high_heteroplasmy"]["mt_encoded_subunit_flux"] < d["mt_reference"]["mt_encoded_subunit_flux"]
assert d["mt_low_RNase_P"]["processed_mRNA_proxy"] < d["mt_reference"]["processed_mRNA_proxy"]
assert d["mt_low_HSPA9"]["imported_nuclear_precursors"] < d["mt_reference"]["imported_nuclear_precursors"]
assert isclose(d["mt_low_HSPA9"]["mt_assembly_ready_subunits"],d["mt_reference"]["mt_assembly_ready_subunits"])
assert d["nuclear_transcription_interaction"]["contrast"] > 0
print("Passed compartment-specific synthetic checks.")
