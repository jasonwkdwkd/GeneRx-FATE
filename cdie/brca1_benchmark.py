"""Predict an independent BRCA1 assay readout with residue-block holdout."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.linear_model import Ridge

HGVS = r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$"
AA3 = "Ala Cys Asp Glu Phe Gly His Ile Lys Leu Met Asn Pro Gln Arg Ser Thr Val Trp Tyr".split()
AA = dict(zip(AA3, "ACDEFGHIKLMNPQRSTVWY"))
GROUPS = ["AVILM", "FWY", "STNQ", "KRH", "DE", "CGP"]
GROUP = {a: g for g, members in enumerate(GROUPS) for a in members}


def read_scores(path, metadata):
    target = json.loads(metadata.read_text())["targetGenes"][0]["targetSequence"]["sequence"]
    scores = pd.read_csv(path)
    variants = scores.hgvs_pro.str.extract(HGVS).set_axis(["ref", "position", "alt"], axis=1)
    variants["score"] = scores.score
    variants = variants.dropna()
    variants = variants[variants.ref.isin(AA) & variants.alt.isin(AA)].copy()
    variants["position"] = variants.position.astype(int)
    assert all(target[p - 1] == AA[a] for a, p in zip(variants.ref, variants.position))
    return target, variants.groupby(["position", "ref", "alt"], as_index=False).score.mean()


def run(hdr, cr, hdr_meta, cr_meta):
    sequence, repair = read_scores(hdr, hdr_meta)
    fragment, growth = read_scores(cr, cr_meta)
    assert sequence.count(fragment) == 1
    offset = sequence.index(fragment)
    growth["position"] += offset
    data = repair.merge(growth, on=["position", "ref", "alt"], suffixes=("_HDR", "_CR"))
    positions = np.sort(data.position.unique())
    folds = np.array_split(positions, 5)
    repair_score = data.score_HDR.to_numpy()
    y = data.score_CR.to_numpy()
    ref, alt = data.ref.map(AA), data.alt.map(AA)
    same = (ref.map(GROUP) == alt.map(GROUP)).astype(float).to_numpy()
    seq_features = np.column_stack([same, *(alt.eq(a) for a in AA.values()),
                                    *(ref.eq(a) for a in AA.values())])
    predictions = {}
    for name, x in {"mean": None, "substitution_only": seq_features,
                    "repair_assay_only": repair_score[:, None],
                    "repair_assay_and_substitution": np.c_[repair_score, seq_features]}.items():
        pred = np.empty_like(y)
        for fold in folds:
            test = data.position.isin(fold).to_numpy()
            pred[test] = y[~test].mean() if x is None else Ridge(alpha=10).fit(
                x[~test], y[~test]).predict(x[test])
        predictions[name] = {"rmse": float(np.sqrt(np.mean((pred-y)**2))),
                             "pearson": None if x is None else float(pearsonr(pred,y).statistic)}
    return {"variant_count": len(data), "residue_count": len(positions),
            "coordinate_offset": offset, "split": "five contiguous blocks of observed residues",
            "endpoint": "cisplatin-resistance growth assay score",
            "input_assay": "homology-directed DNA repair score, available for test variants",
            "scope": "two assays from the same publication, not clinical risk or DNA-only prediction",
            "sha256": {name: sha256(path.read_bytes()).hexdigest() for name, path in
                       zip(("HDR_csv", "CR_csv", "HDR_metadata", "CR_metadata"),
                           (hdr, cr, hdr_meta, cr_meta))}, "models": predictions}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", type=Path, nargs=4)
    print(json.dumps(run(*parser.parse_args().paths), indent=2))
