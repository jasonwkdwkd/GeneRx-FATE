"""Human CDIE extension. Runs a dimensionless, synthetic demonstration."""
import json
from dataclasses import asdict, dataclass
from itertools import product
from math import comb, exp, isclose


def probability(x):
    if not 0 <= x <= 1:
        raise ValueError("Expected a probability between zero and one")
    return x


def gamete_pair(maternal, paternal):
    """One non-recombining locus; probabilities of phased child allele pairs."""
    out = {}
    for pair in product(maternal, paternal):
        out[pair] = out.get(pair, 0) + 1 / (len(maternal) * len(paternal))
    return {"/".join(pair): chance for pair, chance in out.items()}


def inheritance(mother_autosome, father_autosome, mother_x, father_x, father_y):
    return {
        "autosome": gamete_pair(mother_autosome, father_autosome),
        "XX_non_PAR": gamete_pair(mother_x, (father_x,)),
        "XY_non_PAR": gamete_pair(mother_x, (father_y,)),
    }


def bottleneck(h_maternal, effective_copies):
    """Binomial approximation; each outcome gives offspring mtDNA fraction."""
    probability(h_maternal)
    if type(effective_copies) is not int or effective_copies <= 0:
        raise ValueError("Effective bottleneck size must be a positive integer")
    n = effective_copies
    return [(k/n, comb(n, k)*h_maternal**k*(1-h_maternal)**(n-k)) for k in range(n+1)]


def x_dosage(active_maternal, escape):
    """Tissue-specific allele weights for an XX locus."""
    p, q = probability(active_maternal), probability(escape)
    return p+q*(1-p), 1-p+q*p


def nuclear_arrival(locus, origins, fork_speed):
    """Earliest potential arrival from scheduled, bidirectional origins."""
    if not origins or fork_speed <= 0:
        raise ValueError("Supply origins and positive fork speed")
    return min(start + abs(locus-pos)/fork_speed for pos, start in origins)


def mt_replication_times(length, oh, ol, heavy_speed, light_speed, light_delay=0):
    """Circular, strand-displacement scenario with separate strand initiation."""
    if length <= 0 or not (0 <= oh < length and 0 <= ol < length) or oh == ol:
        raise ValueError("Supply distinct locations on the mtDNA circle")
    if min(heavy_speed, light_speed) <= 0 or light_delay < 0:
        raise ValueError("Invalid fork speed or initiation delay")
    light_start = (ol-oh) % length / heavy_speed + light_delay
    return {"heavy_start": 0, "light_start": light_start,
            "complete": max(length/heavy_speed, light_start+length/light_speed)}


@dataclass(frozen=True)
class Regulation:
    dosage: float = 1.
    promoter_strength: float = 1.
    accessibility: float = 1.
    cap_independent: bool = False
    poly_a_independent: bool = False
    tata_present: bool = True
    tata_contribution: float = 0.
    tf_drive: float = 0.
    enhancer_activity: float = 0.
    enhancer_contact: float = 0.
    cap_fraction: float = 1.
    cleavage_poly_a_fraction: float = 1.

    def transcription_flux(self):
        for x in (self.accessibility, self.enhancer_contact, self.cap_fraction,
                  self.cleavage_poly_a_fraction):
            probability(x)
        if min(self.dosage,self.promoter_strength,self.enhancer_activity) < 0:
            raise ValueError("Negative regulatory input")
        drive = self.tf_drive + self.enhancer_activity*self.enhancer_contact
        if self.tata_present:
            drive += self.tata_contribution
        return self.dosage*self.promoter_strength*self.accessibility*exp(drive)

    def exportable_flux(self):
        cap = 1 if self.cap_independent else self.cap_fraction
        tail = 1 if self.poly_a_independent else self.cleavage_poly_a_fraction
        return self.transcription_flux()*cap*tail


@dataclass(frozen=True)
class Isoform:
    splice_fraction: float
    decay: float
    nmd_extra_decay: float
    premature_stop: bool
    initiation_rate: float
    codons: int
    elongation_rate: float
    fold_fraction: float
    delivery_fraction: float
    specific_activity: float
    protein_decay: float

    def output(self, exportable_flux, ribosome_availability=1.):
        for x in (self.splice_fraction, self.fold_fraction, self.delivery_fraction):
            probability(x)
        if (min(self.decay, self.elongation_rate, self.protein_decay, self.codons) <= 0
                or min(self.initiation_rate,self.nmd_extra_decay,self.specific_activity) < 0):
            raise ValueError("Invalid kinetic parameter")
        decay = self.decay+self.nmd_extra_decay*int(self.premature_stop)
        rna = exportable_flux*self.splice_fraction/decay
        cycle_time = self.codons/self.elongation_rate
        initiation = self.initiation_rate*probability(ribosome_availability)
        flux = rna*initiation/(1+initiation*cycle_time)
        folded = flux*self.fold_fraction/self.protein_decay
        delivered = folded*self.delivery_fraction
        return {"mRNA": rna, "translation_flux": flux, "folded_protein": folded,
                "delivered_protein": delivered, "functional_output": delivered*self.specific_activity}


def coding_consequence(net_indel_bases, earlier_stop):
    frame_shift = net_indel_bases % 3 != 0
    return {"frameshift": frame_shift, "premature_stop": earlier_stop,
            "in_frame_indel": net_indel_bases != 0 and not frame_shift,
            "stop_gain_without_frameshift": earlier_stop and not frame_shift}


def first_stop(cds, code="nuclear"):
    """First in-frame stop, with zero-based exclusive end in supplied CDS."""
    stops = {"nuclear": {"TAA","TAG","TGA"},
             "vertebrate_mitochondrial": {"TAA","TAG","AGA","AGG"}}[code]
    cds = cds.upper().replace("U","T")
    if set(cds)-set("ACGT"):
        raise ValueError("Use an unambiguous coding sequence")
    return next((i+3 for i in range(0,len(cds)-2,3) if cds[i:i+3] in stops), None)


def compare_coding_sequences(reference, altered, normal_stop_altered_coordinate):
    """Interpret a preassembled CDS; coordinates must be supplied in altered CDS."""
    altered_stop = first_stop(altered)
    if first_stop(reference) is None or normal_stop_altered_coordinate <= 0:
        raise ValueError("Supply reference CDS with a stop and the mapped stop coordinate")
    earlier = altered_stop is not None and altered_stop < normal_stop_altered_coordinate
    return coding_consequence(len(altered)-len(reference),earlier) | {"stop_end":altered_stop}


def nmd_heuristic(stop_end, last_junction, cutoff=55):
    """Coordinates in the same spliced transcript. None denotes no junction."""
    return last_junction is not None and last_junction - stop_end > cutoff


def nuclear_outputs(regulation, isoforms, free_40s=1., free_60s=1.):
    if not isclose(sum(i.splice_fraction for i in isoforms),1):
        raise ValueError("Isoform fractions must sum to one")
    availability = probability(free_40s)*probability(free_60s)
    return [isoform.output(regulation.exportable_flux(),availability) for isoform in isoforms]


@dataclass(frozen=True)
class TransferRNA:
    amount_per_copy: float
    cca_fraction: float
    modified_fraction: float
    charged_given_cca: float
    modified_relative_decoding: float
    unmodified_relative_decoding: float

    def available(self):
        for x in (self.cca_fraction, self.modified_fraction, self.charged_given_cca,
                  self.modified_relative_decoding, self.unmodified_relative_decoding):
            probability(x)
        if self.amount_per_copy < 0:
            raise ValueError("Negative tRNA abundance")
        decoding = (self.modified_fraction*self.modified_relative_decoding
                    +(1-self.modified_fraction)*self.unmodified_relative_decoding)
        return self.amount_per_copy*self.cca_fraction*self.charged_given_cca*decoding


def mt_output(copies, heteroplasmy, normal_trna, variant_trna, half_supply, nuclear_capacity,
              free_28s=1., free_39s=1.):
    probability(heteroplasmy)
    if min(copies,nuclear_capacity) < 0 or half_supply <= 0:
        raise ValueError("Invalid mitochondrial capacity")
    h = heteroplasmy
    supply = copies*((1-h)*normal_trna.available()+h*variant_trna.available())
    translation = probability(free_28s)*probability(free_39s)*supply/(half_supply+supply)
    return {"usable_tRNA": supply, "mt_translation_capacity": translation,
            "respiratory_capacity_proxy": min(translation,nuclear_capacity)}


def tissue_average(group_weights, values):
    if (len(group_weights) != len(values) or min(group_weights, default=-1) < 0
            or not isclose(sum(group_weights),1)):
        raise ValueError("Tissue weights must have one entry per group and sum to one")
    return sum(w*v for w,v in zip(group_weights,values))


def demo():
    def isoform(**overrides):
        values = dict(splice_fraction=1, decay=1, nmd_extra_decay=0, premature_stop=False,
                      initiation_rate=.1, codons=100, elongation_rate=10, fold_fraction=.9,
                      delivery_fraction=.9, specific_activity=1, protein_decay=1)
        return Isoform(**(values | overrides))

    normal = TransferRNA(1,1,1,.9,1,.2)
    variant = TransferRNA(.4,.9,.3,.6,.8,.1)
    mt = lambda h: mt_output(100,h,normal,variant,40,.95)
    reg = Regulation(tf_drive=.3, enhancer_activity=.5, enhancer_contact=.4,
                     tata_contribution=.1,cap_fraction=.95,cleavage_poly_a_fraction=.9)
    ref = nuclear_outputs(reg,[isoform()])[0]["functional_output"]
    affected = nuclear_outputs(reg,[isoform(delivery_fraction=.1)])[0]["functional_output"]
    rows = {
        "paired_autosomal_and_sex_chromosome_alleles":
            inheritance(("A","a"),("A","a"),("X","x"),"X","Y"),
        "nuclear_fork_arrival": nuclear_arrival(40,[(0,1),(50,4),(100,0)],5),
        "mt_origins": mt_replication_times(16_569,0,11_000,50,50),
        "maternal_mtDNA_offspring_distribution": bottleneck(.5,10),
        "XX_allele_dosage": x_dosage(.8,.1),
        "coding_events": {
            "frameshift": coding_consequence(-1,True),
            "nonsense_scenario": coding_consequence(0,True),
            "in_frame_deletion": coding_consequence(-3,False),
        },
        "nuclear_regulation": asdict(reg) | {"transcription_flux": reg.transcription_flux(),
                                                "exportable_flux": reg.exportable_flux()},
        "synthetic_sequence_stop_gain": compare_coding_sequences("ATGCAATAA","ATGTAATAA",9),
        "nuclear_outputs": {
            "reference": nuclear_outputs(reg,[isoform()]),
            "early_stop_with_NMD": nuclear_outputs(reg,[isoform(
                premature_stop=True,nmd_extra_decay=3,codons=40,specific_activity=.1)]),
            "late_stop_NMD_escape": nuclear_outputs(reg,[isoform(
                premature_stop=True,codons=70,specific_activity=.2)]),
            "altered_splicing": nuclear_outputs(reg,[isoform(splice_fraction=.7),
                isoform(splice_fraction=.3,premature_stop=True,nmd_extra_decay=3,
                        codons=40,specific_activity=.1)]),
            "CFTR_like_trafficking_defect": nuclear_outputs(reg,[isoform(delivery_fraction=.1)]),
            "CFTR_like_gating_defect": nuclear_outputs(reg,[isoform(specific_activity=.1)]),
        },
        "tissue_nuclear_mosaicism": {"variant_cell_fraction":.3,
            "mixed_function":tissue_average([.7,.3],[ref,affected])},
        "mt_translation_by_heteroplasmy": {str(h):mt(h) for h in (0,.5,1)},
        "mt_tissue_same_mean_different_distribution": {
            "all_cells_h_0.5":mt(.5)["respiratory_capacity_proxy"],
            "half_WT_half_homoplasmic_variant":tissue_average(
                [.5,.5],[mt(h)["respiratory_capacity_proxy"] for h in (0,1)]),
        },
    }
    return rows


if __name__ == "__main__":
    print(json.dumps(demo(),indent=2))
