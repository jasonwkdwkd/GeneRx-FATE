# Reproduction and evaluation protocol

This protocol separates code verification, parameter calibration, and testing on held-out measurements. The issued expansion currently completes only code verification. For the implemented formulas see [EQUATIONS.md](EQUATIONS.md), for biological mapping see [METHODS.md](METHODS.md) and [Foundation.md](Foundation.md), and for dataset access see [DATA_SOURCES.md](DATA_SOURCES.md).

## Reproduce the synthetic outputs

Use Python 3.10 or later. The expansion and its checks require only the standard library. Run these commands from this directory on a pinned Git revision.

```sh
git rev-parse HEAD
python --version
python check.py
python check_compartments.py
python audit_model.py
python compartments.py > reproduced_output.json
diff -u example_output.json reproduced_output.json
```

The checks exercise inheritance bounds and expected segregation, replication timing, cap/poly(A) dependence, stop position and NMD separation, compartment-specific RNA processing, tRNA availability, import and chaperone effects, and tissue aggregation. The audit checks a non-identifiability example and shows the effect of several isolated input changes. Differences from a Python build's floating point math should be assessed numerically before concluding that a model changed.

## Inputs, dimensions and model boundaries

| Inputs | Source in the released model | Measurement needed for a fitted model |
|---|---|---|
| Nuclear dosage, accessibility, CCAAT/TATA presence, NF-Y/TBP/Pol II drive, enhancer contact and repressor drive | Supplied synthetic settings | Genotype and copy number, chromatin accessibility, occupancy, enhancer contact, and nascent transcription in a defined tissue |
| Splice fractions, RNA turnover and NMD, cap/poly(A) processing | Supplied per-isoform fractions and dimensionless steady-state rate settings | Junction or isoform counts, independent RNA turnover and processing measurements; distinguish estimated and measured values |
| Ribosome availability, initiation/elongation, protein folding and delivery | Supplied synthetic efficiency or rate settings | Matched ribosome profiling and protein quantity, localization and functional measurements |
| mtDNA origins, fork rates, copy number and heteroplasmy | Origin timing and expression use separate synthetic inputs | Strand-specific replication evidence and time-resolved mtDNA copy number and allele fraction under matched conditions |
| mtRNA processing, tRNA abundance, modification, charging and decoding | Supplied rates or probabilities | Precursor/product abundance, modification-aware tRNA and charged/uncharged measurements plus translation readouts |
| Nuclear subunit import, chaperone folding and mixed complex assembly | Synthetic normalized supply terms | Subunit import, assembly and complex activity; measure respiration/ATP separately before modeling those outcomes |

No experimental measurement units were assigned to the synthetic coefficients. The locus, allele, isoform and tissue choices must be stated for every application. Replication time and mtDNA copy number are separate in the released code. The assembly output is a normalized proxy, and unobserved stages are not independently identifiable from a final phenotype alone.

## Prespecified protocol for real-data calibration and testing

1. **Define the endpoint and comparators.** For mitochondrial work, select an observed endpoint such as oxygen consumption and a simpler comparator such as heteroplasmy-only regression. For CFTR, define an assay-specific endpoint such as chloride current and compare against a model using protein abundance alone. State the prediction unit, assay units and direction in advance.
2. **Freeze input and splitting rules.** Version and hash downloaded inputs and record exclusions, replicate identities, batch, clone, donor, allele, and tissue. Keep related measurements together when splitting. Use donor-, clone-, variant- or residue-group holdouts as appropriate. Reserve the outer test set from model fitting and threshold selection.
3. **Fit only identifiable parameters.** Set assay-derived priors or constraints for directly observable components. Fit remaining parameters using inner validation and report intervals or profile likelihoods. Use separate rate/efficiency terms only when the available observations can distinguish them. If the data identify a product of efficiencies, report that product.
4. **Score matched baselines and ablations.** Compare identical inputs and held-out cases, report per-case errors and sample sizes, and remove stages one at a time. Report prediction intervals where justified. A failed direction or worse performance is a result to retain.
5. **Report scope.** Keep synthetic examples, standalone assay projections, and prospective or independent validation separately labeled. Publish software versions, split assignments, parameter provenance, and evaluation code alongside the outputs.

## Interaction interpretation

The factorial contrast is calculated on a prespecified endpoint and scale across all four conditions. An interaction on a measured molecular endpoint does not itself establish a causal interaction, a clinical diagnosis, or a personalized risk.
