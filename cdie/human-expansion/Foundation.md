# Central Dogma Interaction Equation — human genetics extension
     
Run with Python 3.10 or newer. The extension is an independent, steady-state companion to the published time-dependent prototype bundled as original_cdie.py. The original and the extension have separate entry points.

    python model.py
    python check.py

## What is computed

- Inheritance enumerates the possible phased child allele pairs for one autosomal locus and for a non-pseudoautosomal sex-linked locus in an XX/XY cross. X-inactivation at an XX locus is modeled as the fraction of cells expressing each parental X, with an optional escape fraction. The mitochondrial module gives the binomial distribution of offspring heteroplasmy conditional on a maternal germline heteroplasmy and effective bottleneck size. The same inheritance calculation may be applied stepwise to subsequent generations. Transmission, tissue expression, and clinical penetrance are different quantities.
- Replication computes potential earliest fork arrival from scheduled nuclear origin zones. A separate illustrative mitochondrial strand-displacement model uses O_H and O_L, circular distance, fork speed, and light-strand initiation delay. Demonstration origin coordinates and rates are synthetic. This is a timing calculation, not a fitted copy-number or mutation-rate model. mtDNA copy number and heteroplasmy are supplied separately to the expression calculation.
- Transcription uses dosage, promoter strength, accessibility, optional TATA contribution, transcription-factor drive, and enhancer activity and contact. Cap and cleavage/polyadenylation are separate effective fractions with optional bypass switches for modeling cases outside the standard cap-dependent, polyadenylated transcript. This is a multiplicative approximation. The separate original time-course code includes histone modifications, CpG methylation, GR and GRE dynamics.
- RNA and protein use per-isoform splice proportions, RNA decay and NMD rates, initiation, ribosome availability, coding length, elongation, folding, delivery, activity, and protein turnover. Simple CDS comparisons distinguish net frame changes and stop gain for preassembled sequences. The NMD heuristic is separate from the rate, which must be measured or supplied. It is not an automatic consequence predictor for clinical variants.
- Mitochondrial tRNA supply distinguishes wild-type and mutant mtDNA copies, tRNA abundance, 3′-CCA completion, a simplified modification state, conditional aminoacylation, and decoding efficiency. Mitochondrial translation uses separate 28S/39S availability and a saturating tRNA term. The respiratory output uses a limiting-capacity proxy with nuclear-encoded capacity supplied as an input. It is not a fitted oxidative phosphorylation, oxygen flux, or ATP model.
- Tissue mixtures weight separately calculated functional outputs. Nuclear mosaic cell fraction and cell-level mtDNA heteroplasmy are distinct. The example shows two tissues with identical mean heteroplasmy and different functional proxies.

## Scope and validation

The numerical examples are synthetic. The checks establish expected behavior and reproducibility on those examples. The two scripts do not fit data or predict individual disease risk. The code does not infer allelic phase, variant effects, enhancer maps, actual splice proportions, germline mosaicism, or tissue composition from sequence. The inheritance examples assume a single non-recombining locus per chromosome class and typical XX/XY segregation. They do not simulate meiosis, recombination, aneuploidy, translocations or pseudoautosomal transmission. The mitochondrial stop lookup uses nominal vertebrate mitochondrial code stop assignments and does not model site-specific termination mechanisms or completion of partial stops by polyadenylation. The assays needed to fit individual rates include RNA junction counts, RNA decay, protein abundance, compartment-specific translation, and independent functional measurements.

## Selected references

- Petryk et al. Replication landscape of the human genome. doi:10.1038/ncomms10208
- Holt and Reyes. Human mitochondrial DNA replication. doi:10.1101/cshperspect.a012971
- Tukiainen et al. Landscape of X chromosome inactivation across human tissues. doi:10.1038/nature24265
- Jaganathan et al. Predicting splicing from primary sequence with deep learning. doi:10.1016/j.cell.2018.12.015
- Abou Tayoun et al. PVS1 variant interpretation, including nonsense-mediated decay. PMID:30192042
- Sosnay et al. Defining the disease liability of variants in CFTR. doi:10.1038/ng.2745
- Post-transcriptional methylation of mitochondrial-tRNA differentially contributes to mitochondrial pathology. doi:10.1038/s41467-024-53318-x
- Aryaman et al. Mitochondrial DNA density homeostasis and heteroplasmy. PMCID:PMC5705840
