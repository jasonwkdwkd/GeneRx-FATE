# CDIE human nuclear and mitochondrial compartments

Requires Python 3.10 or newer for the two synthetic demonstrations and checks. Run in an extracted directory.

    python compartments.py
    python check_compartments.py
    python model.py
    python check.py

## What this extension adds

| Module | Inputs and outputs |
|---|---|
| Human nuclear DNA | Separate optional CCAAT/NF-Y and TATA/TBP occupancy terms. Enhancer contact, accessibility, repressors, and RNA polymerase II determine precursor flux. Cap and cleavage/polyadenylation fractions act on RNA. Kozak efficiency acts on translation initiation. Each supplied splice isoform retains its individual decay, premature-stop and protein parameters. |
| Mitochondrial replication | POLG, TWINKLE and mtSSB/SSBP1 scale an illustrative fork capacity. POLRMT priming affects the timing proxy, which uses separate heavy and light strand origins. |
| Mitochondrial gene expression | LSP and HSP initiation uses POLRMT, TFAM and TFB2M. TEFM acts on productive transcription. RNase P and ELAC2 act on tRNA-flanked precursor processing. Non-tRNA-flanked processing has a separate supplied efficiency. MTPAP acts on the fraction of mRNAs with incomplete stops. Normal and variant mtDNA copies contribute separately to mature usable tRNA and 28S/39S-limited translation. |
| Mitochondrial protein handling | A representative nuclear-encoded precursor uses TOM, TIM23 and mtHsp70/HSPA9 for the presequence pathway, then HSP60/HSP10 for folding of a chaperonin-dependent substrate. A separate insertion/assembly fraction represents mtDNA-encoded subunits. The limiting pool gives a normalized mixed-complex capacity proxy. |

Each occupancy is represented by L/(K+L). The demo uses dimensionless inputs and default K=1. A factorial contrast is evaluated for NF-Y and enhancer activity on precursor flux. In the mitochondrial example, the heavy-strand pool supplies the modeled MT-TL1 tRNA. The full transcript mixture is summarized as a single processed mRNA pool. Regulatory drive and serial efficiencies are chosen model forms rather than measured kinetic laws.

## Validation scope

The supplied checks demonstrate directional responses to CCAAT loss, reduced SSBP1, impaired RNA processing, heteroplasmy and reduced HSPA9. They check the intended behavior of these formulas on synthetic scenarios. Every rate, capacity, and output is illustrative. The mitochondrial replication calculation does not infer copy number. The mitochondrial expression calculation accepts the copy number separately. The complex output is neither measured respiration nor ATP production. The example does not infer promoter occupancy, splice outcomes, or variant-specific effects from genomic coordinates. Nuclear locus, ancestral inheritance, and tissue-mixture functions remain documented in Foundation.md. The original time-course script is included for reference, with a separate entry point.

## Biological references

- NF-Y binding to CCAAT: PMID 23332751
- Hillen et al. Structural Basis of Mitochondrial Transcription Initiation. doi:10.1016/j.cell.2017.10.036
- Ojala et al. tRNA punctuation model of RNA processing in human mitochondria. doi:10.1038/290470a0
- Holt and Reyes. Human mitochondrial DNA replication. doi:10.1101/cshperspect.a012971
- Human mtSSB: UniProt Q04837
- Human mitochondrial RNA processing and modifications: PMCID PMC8348895
- Importing mitochondrial proteins: machineries and mechanisms. PMID 19703392
