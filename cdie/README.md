# Central Dogma Interaction Equation (CDIE)

Research prototype for DNA variation, gene regulation, expression, and molecular function. The executable forward model is a single human locus with two RNA branches and a synthetic genotype perturbation. The public-data scripts benchmark measured RNA and protein-function projections separately.

## Run

Python 3.11+ with NumPy, SciPy, pandas, and scikit-learn.

```sh
python cdie.py --variant variant_example.json --output simulation
python rna_benchmark.py GSE104714_DEX_EtOH_expression.tsv.gz
python brca1_benchmark.py BRCA1_HDR.csv BRCA1_CR.csv BRCA1_HDR_metadata.json BRCA1_CR_metadata.json
```

- `simulation/trajectories.csv` contains both genotype simulations and four input conditions.
- `simulation/interactions.csv` contains each output's A-by-B interaction and the difference of that interaction between synthetic genotype and reference.
- `A549_results.json` and `BRCA1_results.json` contain executable benchmark outputs and SHA-256 provenance.

## Public inputs

- A549 DEX and vehicle time course. https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104714/suppl/GSE104714_DEX_EtOH_expression.tsv.gz
- BRCA1 homology-directed repair score. https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00001222-b-2/scores
- BRCA1 cisplatin resistance score. https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00001222-a-2/scores
- Use the same two MaveDB URLs without `/scores` for JSON metadata with target protein sequences.

The simulator accepts parameter shifts, which a separately calibrated DNA-to-parameter model would have to predict for individual variants. The supplied shifts are synthetic. The code has no personal disease risk estimator. Public aggregate data and functional assay scores do not provide prospective individual outcomes. Neither clinical labels nor functional assay scores were used as surrogate incidence data.
