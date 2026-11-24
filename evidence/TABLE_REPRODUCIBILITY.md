# Manuscript table reproducibility

| Table | Public source |
|---|---|
| 1. Analytic cohort and evaluation designs | [Pair mapping](../S_pair_mapping.csv), [source integrity](../S_source_pairs_integrity.json), and the methods in the manuscript. |
| 2. Independent reconstruction | [Reconstruction error by pair](S_reconstruction_log1p.csv) and [reconstruction script](../reproduce.py). |
| 3. Repeated strict gene-holdout ablations | [Unique-pair archive](../S_frozen_observed_unique_pairs.csv), [appearance-level archive](../S_frozen_ablation_repeated_instances.csv), and [benchmark](../benchmark.py). |
| 4–5. Conditional feature-shuffling | [Archived null distribution](S_frozen_permutation_distribution.csv), [reproduced summary](permutation_reproduced.json), and [permutation script](../permutation.py). The manuscript states the seed and scope for each run. |
| 6. Spectral gene-disjoint evaluation | [Fold-level results](S_frozen_ablation_spectral_by_fold.csv) and [pooled results](S_frozen_ablation_spectral_summary.csv). These files are archived outputs. |
| 7. Matched K562 and THP-1 analyses | [K562 benchmark results](../benchmark.json) and [THP-1 replication archive](GeneRx_FATE_THP1_Replication.zip), including code, GEO input hashes, split counts, and per-pair errors. [Archive SHA-256](THP1_ARCHIVE_SHA256SUMS). |
| 8. Illustrative disease contexts | Literature citations in the manuscript. This table does not report a computational benchmark. |
| 9. Exploratory CDIE projections | See the [CDIE data-source register](../cdie/human-expansion/DATA_SOURCES.md), [A549 output](../cdie/A549_results.json) and [BRCA1 output](../cdie/BRCA1_results.json). These linked outputs reproduce selected components of the table. |

Source expression matrices are obtained from GEO as documented in the scripts and archive README. See the [repository README](../README.md) for software requirements and reconstruction instructions.
