"""Compare six-feature direct ridge and rank-one regression on frozen holdout rules."""
import argparse
import csv
import json
from ast import literal_eval
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

import generx_fate as gf


def reconstruct(groups, raw_sums, mapping):
    index = {key:i for i,key in enumerate(groups["targets"])}
    size = len(index)
    sums = np.fromfile(raw_sums, dtype=np.float64).reshape(size, 33694)
    means = sums / np.array([groups["counts"][key] for key in index])[:,None]
    control = means[index["CONTROL"]]
    with mapping.open() as f:
        pair_labels = sorted({row["target"] for row in csv.DictReader(f)
                              if row["in_frozen_pairs"] == "True"})
    pairs = np.array([literal_eval(x) for x in pair_labels])
    a = np.stack([means[index[u]] - control for u in pairs[:,0]])
    b = np.stack([means[index[v]] - control for v in pairs[:,1]])
    y = np.stack([means[index[label]] - control for label in pair_labels])
    assert pairs.shape == (131, 2) and len(np.unique(pairs)) == 73
    return a,b,y,pairs


def evaluate(a,b,y,pairs, archive):
    add = a+b
    resid = y-add
    x = gf.feature_matrix(a,b)
    mask = gf.panels(a,b)
    genes = np.unique(pairs)
    output = {}
    with archive.open() as f:
        frozen = list(csv.DictReader(f))
    archived = {(r["pair_a"],r["pair_b"]):r for r in frozen}
    for seed in range(100,120):
        held = np.random.default_rng(seed).choice(genes, size=18, replace=False)
        indicator = np.isin(pairs,held)
        tr,te = np.flatnonzero(~indicator.any(1)), np.flatnonzero(indicator.all(1))
        base = gf.Baseline("B3").fit(add[tr],resid[tr],mask[tr])
        base_pred = add[te]+base.predict(add[te])
        targets = np.empty_like(resid[tr])
        for st,sv in KFold(n_splits=5, shuffle=True, random_state=seed).split(tr):
            train,valid = tr[st],tr[sv]
            sub = gf.Baseline("B3").fit(add[train],resid[train],mask[train])
            targets[sv] = resid[valid]-sub.predict(add[valid])
        scaler = StandardScaler().fit(x[tr])
        xr, xt = scaler.transform(x[tr]),scaler.transform(x[te])
        pca=PCA(n_components=1,svd_solver="full").fit(targets)
        predictions = [add[te],base_pred,base_pred+targets.mean(axis=0),
                       base_pred+pca.inverse_transform(Ridge(alpha=100).fit(xr,pca.transform(targets)).predict(xt).reshape(-1,1)),
                       base_pred+Ridge(alpha=100).fit(xr,targets).predict(xt)]
        sses = np.array([gf._sse(y,pred,mask,te) for pred in predictions]).T
        for idx,sse in zip(te,sses):
            output.setdefault(tuple(pairs[idx]),[]).append(sse)
    assert set(output)==set(archived)
    assert sum(map(len,output.values()))==172 and len(output)==92
    means=np.stack([np.mean(output[k],axis=0) for k in sorted(output)])
    expected=np.array([[float(archived[k][col]) for col in
                        ("sse_add","sse_b3","sse_mean","sse_full_observed")] for k in sorted(output)])
    deltas=np.max(np.abs(means[:,:4]-expected),axis=0)
    total=means.sum(axis=0)
    return {"unique_pairs":len(output),"test_appearances":172,"sum_sse":total.tolist(),
            "additive_relative_r2":(1-total/total[0]).tolist(),
            "labels":["additive","B3","mean_only","rank_one_PCA_ridge","direct_ridge"],
            "max_absolute_per_pair_sse_difference_vs_archive":deltas.tolist(),
            "full_vs_direct_sse_reduction":float(1-total[3]/total[4])}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("groups",type=Path)
    parser.add_argument("sums",type=Path)
    parser.add_argument("mapping",type=Path)
    parser.add_argument("archive",type=Path)
    arg=parser.parse_args()
    print(json.dumps(evaluate(*reconstruct(json.loads(arg.groups.read_text()), arg.sums, arg.mapping),arg.archive),indent=2))
