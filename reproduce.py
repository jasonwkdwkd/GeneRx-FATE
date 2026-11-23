"""Rebuild GEO inputs and rerun the direct-ridge comparison."""
import argparse
import csv
import gzip
import json
import subprocess
import numpy as np
from collections import Counter
from pathlib import Path

from benchmark import evaluate, reconstruct

ROOT = Path(__file__).resolve().parent


def run(source, output):
    output.mkdir(parents=True, exist_ok=True)
    prefix = "GSE133344_filtered_"
    mapping = ROOT / "S_pair_mapping.csv"
    with mapping.open() as handle:
        labels = {r["condition"]: ("CONTROL" if r["category"] == "control" else r["target"])
                  for r in csv.DictReader(handle)}
    with gzip.open(source / f"{prefix}cell_identities.csv.gz", "rt") as handle:
        cells = {r["cell_barcode"]: labels[r["guide_identity"]] for r in csv.DictReader(handle)
                 if r["good_coverage"] == "True"}
    with gzip.open(source / f"{prefix}barcodes.tsv.gz", "rt") as handle:
        codes = [cells.get(line.strip()) for line in handle]
    counts = Counter(value for value in codes if value is not None)
    assert len(codes) == 111668 and sum(counts.values()) == 104507 and counts["CONTROL"] == 11183
    targets = sorted(counts)
    ids = {key: i for i, key in enumerate(targets)}
    labels_path = output / "groups.txt"
    labels_path.write_text("".join(f"{ids.get(value, -1)}\n" for value in codes))
    sums = output / "sums.bin"
    program = output / "aggregate"
    subprocess.run(["g++", "-O3", str(ROOT/"aggregate.cpp"), "-o", str(program)], check=True)
    with subprocess.Popen(["gzip", "-dc", str(source / f"{prefix}matrix.mtx.gz")],
                          stdout=subprocess.PIPE) as gzip_process:
        subprocess.run([str(program), str(len(targets)), str(labels_path), str(sums)],
                       stdin=gzip_process.stdout, check=True)
        if gzip_process.wait():
            raise ValueError("GEO matrix decompression failed")
    arrays = reconstruct(dict(targets=targets, counts=counts), sums, mapping)
    np.savez_compressed(output / "pairs_reconstructed.npz", **dict(zip(("a", "b", "y", "pairs"), arrays)))
    result = evaluate(*arrays, ROOT / "S_frozen_observed_unique_pairs.csv")
    assert max(result["max_absolute_per_pair_sse_difference_vs_archive"]) < 1e-4
    (output/"benchmark.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geo_directory", type=Path)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    run(args.geo_directory, args.output_directory)
