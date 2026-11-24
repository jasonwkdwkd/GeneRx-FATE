"""Dimensionless human nuclear and mitochondrial CDIE example."""
import json
from dataclasses import dataclass
from math import exp

from model import Isoform, TransferRNA, mt_replication_times, probability


def occupancy(level, half_saturation=1.):
    if level < 0 or half_saturation <= 0:
        raise ValueError("Use a nonnegative level and positive half-saturation")
    return level / (half_saturation + level)


@dataclass(frozen=True)
class NuclearTranscription:
    dosage: float = 2.
    accessibility: float = .8
    ccaat_present: bool = True
    tata_present: bool = True
    nf_y: float = 1.
    tbp: float = 1.
    pol_ii: float = 2.
    enhancer_activity: float = .5
    enhancer_contact: float = .4
    repressor_drive: float = 0.
    cap_fraction: float = .95
    cleavage_poly_a_fraction: float = .9
    kozak_efficiency: float = .9

    def outputs(self, transcripts):
        if self.dosage < 0 or self.enhancer_activity < 0 or self.repressor_drive < 0:
            raise ValueError("Use nonnegative regulatory inputs")
        if abs(sum(tx.splice_fraction for tx in transcripts)-1.) > 1e-8:
            raise ValueError("Isoform fractions must sum to one")
        ccaat = occupancy(self.nf_y) if self.ccaat_present else 0.
        tata = occupancy(self.tbp) if self.tata_present else 0.
        upstream_drive = self.enhancer_activity * probability(self.enhancer_contact)
        precursor = (self.dosage * probability(self.accessibility) * occupancy(self.pol_ii)
                     * exp(ccaat + tata + upstream_drive - self.repressor_drive))
        exportable = precursor*probability(self.cap_fraction)*probability(self.cleavage_poly_a_fraction)
        return {"NF_Y_CCAAT_occupancy": ccaat, "TBP_TATA_occupancy": tata,
                "pre_mRNA_flux": precursor, "exportable_RNA_flux": exportable,
                "isoforms": [tx.output(exportable, probability(self.kozak_efficiency))
                             for tx in transcripts]}


@dataclass(frozen=True)
class MitochondrialReplication:
    polg: float = 1.
    twinkle: float = 1.
    ssbp1: float = 1.
    polrmt_priming: float = 1.
    origin_h: int = 0
    origin_l: int = 11_000
    genome_length: int = 16_569
    max_fork_speed: float = 100.

    def outputs(self):
        capacity = occupancy(self.polg)*occupancy(self.twinkle)*occupancy(self.ssbp1)
        priming = occupancy(self.polrmt_priming)
        speed = self.max_fork_speed*capacity
        timing = (mt_replication_times(self.genome_length,self.origin_h,self.origin_l,
                                       speed,speed,light_delay=1/priming)
                  if speed > 0 and priming > 0 else None)
        return {"fork_capacity":capacity,"priming_capacity":priming,"relative_timing":timing}


@dataclass(frozen=True)
class MitochondrialExpression:
    polrmt: float = 2.
    tfam: float = 2.
    tfb2m: float = 2.
    tefm: float = 2.
    lsp_strength: float = 1.
    hsp_strength: float = 1.
    rnase_p: float = 2.
    elac2: float = 2.
    nonpunctuated_processing: float = .8
    punctuated_fraction: float = .9
    mtpap: float = 2.
    incomplete_stop_fraction: float = .3
    free_28s: float = 1.
    free_39s: float = 1.
    half_tRNA_supply: float = 1.

    def outputs(self, copies, heteroplasmy, normal_tRNA, variant_tRNA):
        if min(copies,self.lsp_strength,self.hsp_strength) < 0 or self.half_tRNA_supply <= 0:
            raise ValueError("Invalid copy number, promoter strength or tRNA half-saturation")
        h = probability(heteroplasmy)
        initiation = copies*occupancy(self.polrmt)*occupancy(self.tfam)*occupancy(self.tfb2m)
        lsp, hsp = initiation*self.lsp_strength, initiation*self.hsp_strength
        elongation = occupancy(self.tefm)
        cleavage = occupancy(self.rnase_p)*occupancy(self.elac2)
        fraction = probability(self.punctuated_fraction)
        processing = fraction*cleavage+(1-fraction)*probability(self.nonpunctuated_processing)
        # The worked example treats MT-TL1 tRNA as heavy-strand encoded.
        wt_tRNA = hsp*elongation*cleavage*(1-h)*normal_tRNA.available()
        var_tRNA = hsp*elongation*cleavage*h*variant_tRNA.available()
        tRNA = wt_tRNA+var_tRNA
        partial_stop = probability(self.incomplete_stop_fraction)
        stop_completion = 1-partial_stop+partial_stop*occupancy(self.mtpap)
        rna = (lsp+hsp)*elongation*processing
        ribosomes = probability(self.free_28s)*probability(self.free_39s)
        translated = rna*stop_completion*ribosomes*tRNA/(self.half_tRNA_supply+tRNA)
        return {"LSP_primary_output":lsp,"HSP_primary_output":hsp,
                "punctuation_efficiency":cleavage,"processed_mRNA_proxy":rna,
                "stop_completion_fraction":stop_completion,
                "usable_WT_tRNA":wt_tRNA,"usable_variant_tRNA":var_tRNA,
                "mt_encoded_subunit_flux":translated}


@dataclass(frozen=True)
class MitochondrialProteinHandling:
    tom: float = 2.
    tim23: float = 2.
    hspa9: float = 2.
    hsp60: float = 2.
    hsp10: float = 2.
    membrane_potential_fraction: float = .9
    mt_insert_assembly_fraction: float = .8

    def outputs(self, nuclear_precursors, mt_subunits):
        if min(nuclear_precursors,mt_subunits) < 0:
            raise ValueError("Use nonnegative precursor amounts")
        imported = (nuclear_precursors*occupancy(self.tom)*occupancy(self.tim23)
                    *occupancy(self.hspa9)*probability(self.membrane_potential_fraction))
        matrix_folded = imported*occupancy(self.hsp60)*occupancy(self.hsp10)
        mt_assembly_ready = mt_subunits*probability(self.mt_insert_assembly_fraction)
        return {"imported_nuclear_precursors":imported,"folded_nuclear_subunit_proxy":matrix_folded,
                "mt_assembly_ready_subunits":mt_assembly_ready,
                "mixed_complex_capacity_proxy":min(matrix_folded,mt_assembly_ready)}


def interaction(run, a, b):
    """Apply the same outcome definition to all four factorial arms."""
    y0, ya, yb, yab = (run(x,y) for x,y in ((0,0),(a,0),(0,b),(a,b)))
    return {"control":y0,"A":ya,"B":yb,"AB":yab,"contrast":yab-ya-yb+y0}


def demo():
    tx = Isoform(1,1,0,False,.1,100,10,.9,.9,1,1)
    wt, var = TransferRNA(1,1,1,.9,1,.2), TransferRNA(.4,.9,.3,.6,.8,.1)
    mt_proteins = MitochondrialProteinHandling()

    def mt(h, rnase_p=2.):
        out = MitochondrialExpression(rnase_p=rnase_p).outputs(100,h,wt,var)
        return out | mt_proteins.outputs(200,out["mt_encoded_subunit_flux"])

    def contrast(a,b):
        return NuclearTranscription(nf_y=1+a,enhancer_activity=b).outputs([tx])["pre_mRNA_flux"]

    return {"nuclear_reference":NuclearTranscription().outputs([tx]),
            "nuclear_CCAAT_loss":NuclearTranscription(ccaat_present=False).outputs([tx]),
            "nuclear_transcription_interaction":interaction(contrast,1,.5),
            "mt_replication_reference":MitochondrialReplication().outputs(),
            "mt_replication_low_SSBP1":MitochondrialReplication(ssbp1=.2).outputs(),
            "mt_reference":mt(0), "mt_high_heteroplasmy":mt(.8),
            "mt_low_RNase_P":mt(0,rnase_p=.2),
            "mt_low_HSPA9":MitochondrialProteinHandling(hspa9=.2).outputs(
                200,mt(0)["mt_encoded_subunit_flux"])}


if __name__ == "__main__":
    print(json.dumps(demo(), indent=2))
