"""A549 held-out-time benchmark using GEO processed expression."""
import argparse
import gzip
import json
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.special import gammainc

def run(path):
    times = np.arange(1, 12, 2)
    with gzip.open(path, "rt") as f:
        header = f.readline().strip().split("\t")
        columns = [header.index(f"{arm}_{t}hr_Rep{r}")
                   for arm in ("DEX", "EtOH") for t in times for r in (1, 2)]
        counts = np.loadtxt(f, delimiter="\t", usecols=columns).reshape(-1, 2, 6, 2)
    keep = counts[:, :, :3, :].sum(axis=(1, 2, 3)) >= 120
    means = np.log2(counts[keep] / counts.sum(axis=0) * 1e6 + 0.5).mean(axis=-1)
    train, test = np.hsplit(means[:, 0] - means[:, 1], 2)
    prediction = dict(unchanged=0, mean_early=train.mean(axis=1, keepdims=True),
                      last_early=train[:, -1:])
    rates = np.geomspace(.025, 4, 101)
    def fit(order, count):
        shape = gammainc(order, rates[:, None] * times)
        early = shape[:, :count]
        projection = train[:, :count] @ early.T
        norm = np.sum(early ** 2, axis=1)
        best = np.argmax(np.sum(projection ** 2, axis=0) / norm)
        return projection[:, best, None] / norm[best] * shape[best]
    for order in (1, 2):
        prediction[f"stage_{order}"] = fit(order, 3)[:, 3:]
    # Choose blend strength using 5 h alone. This candidate was designed after
    # seeing first-round results and requires a fresh study for confirmation.
    empirical_5h = train[:, :2].mean(axis=1)
    delta_5h = fit(2, 2)[:, 2] - empirical_5h
    weight = np.clip(np.dot(delta_5h, train[:, 2] - empirical_5h) /
                     np.dot(delta_5h, delta_5h), 0, 1)
    prediction["exploratory_blend"] = (1-weight) * prediction["mean_early"] + weight * prediction["stage_2"]
    return dict(sha256=sha256(path.read_bytes()).hexdigest(), eligible=int(keep.sum()), exploratory_kinetic_weight=float(weight),
                rmse={k: float(np.sqrt(np.mean((v-test)**2))) for k, v in prediction.items()})

if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("path", type=Path)
    print(json.dumps(run(cli.parse_args().path), indent=2))
