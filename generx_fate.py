"""GeneRx-FATE v0.3 — research scaffold for nonadditive CRISPRa response prediction.

NOT a validated therapeutic predictor. No biological or clinical claim follows
from running this code.

Changes from v0.2 (all driven by the 2026-11-10 Colab run on GSE133344)
-----------------------------------------------------------------------
1. PCA mean handling refactored: the PCA object now restores the residual mean
   itself, instead of the caller subtracting and re-adding it. Verified
   numerically identical to v0.2 (max abs difference ~1e-15), so this changes
   no result — it only removes the hazard of a double mean addition.
2. Rank-1 coefficient reshape made robust across sklearn versions.
3. residual_mode: 'crossfit' (default, as in v0.2) or 'infold'. The Colab run
   used in-fold residuals as the PCA/ridge target. Both are now selectable so
   the choice can be reported as a robustness check rather than silently made.
4. FIXED, v0.2 bug: permutation_control(block='all') permuted a and b, which
   also permutes the additive baseline and the evaluation panels, inflating the
   denominator of R-squared. The null was therefore not comparable to the
   observed value. Permutation now shuffles rows of the feature matrix only,
   leaving a, b, y, the panels and the baselines untouched.
5. Permutation control re-runs inner model selection inside every permutation
   by default, and defaults to 200 permutations.
6. gene_disjoint_folds(method='spectral'): partitions genes to keep more pairs
   inside a block, raising the fraction of pairs that can be tested under the
   strict no-shared-gene constraint. method='random' reproduces v0.2 behaviour.
7. repeat_cv() over seeds and bootstrap_ci() on pooled R-squared.

RNA-fate features ("extra") must be measured independently of the held-out
double-perturbation outcomes. Anything derived from y is leakage.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

EPS = 1e-9
PANEL_SIZE = 200
N_BASE_FEATURES = 6


# ----------------------------------------------------------------- features


def single_features(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Six order-invariant features of a gene pair's two single responses."""
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    s = a + b
    ns2 = float(s @ s)
    return np.array([
        np.log1p(min(na, nb)),
        np.log1p(max(na, nb)),
        np.log1p(np.sqrt(ns2)),
        float(a @ b) / (na * nb + EPS),
        np.log1p(np.linalg.norm(a - b)),
        float(a @ b) / (ns2 + EPS),
    ])


def feature_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.vstack([single_features(ai, bi) for ai, bi in zip(a, b)])


def panels(a: np.ndarray, b: np.ndarray, size: int = PANEL_SIZE) -> np.ndarray:
    """Boolean mask (pairs, genes): top-`size` genes by |dA| + |dB| per pair.

    Depends only on single-perturbation responses, so it is computable for a
    held-out pair without seeing its double-perturbation outcome.
    """
    score = np.abs(a) + np.abs(b)
    size = min(size, a.shape[1])
    idx = np.argpartition(-score, size - 1, axis=1)[:, :size]
    m = np.zeros(a.shape, dtype=bool)
    np.put_along_axis(m, idx, True, axis=1)
    return m


# ----------------------------------------------------------------- baselines


@dataclass
class Baseline:
    """Gamma-hat = mean_resid * use_mean + coef . additive-terms.

    kind: 'B0' zero, 'B1' training mean, 'B2a' scaled additive (all genes),
    'B2b' scaled additive (panel-fitted), 'B3' mean + scaled additive,
    'B3b' B3 with scaling varying in log additive magnitude.
    """

    kind: str = "B3"
    mean_resid: np.ndarray | None = None
    coef: np.ndarray = field(default_factory=lambda: np.zeros(0))

    def fit(self, add: np.ndarray, resid: np.ndarray, mask: np.ndarray) -> "Baseline":
        k = self.kind
        self.mean_resid = resid.mean(axis=0) if k in ("B1", "B3", "B3b") else np.zeros(add.shape[1])
        if k in ("B0", "B1"):
            self.coef = np.zeros(0)
            return self
        w = mask if k in ("B2b", "B3", "B3b") else np.ones_like(mask, dtype=bool)
        target = (resid - self.mean_resid)[w]
        cols = [add[w]]
        if k == "B3b":
            mag = np.log1p(np.abs(add)).mean(axis=1, keepdims=True) * np.ones_like(add)
            cols.append((add * mag)[w])
        X = np.column_stack(cols)
        self.coef = np.linalg.lstsq(X, target, rcond=None)[0]
        return self

    def predict(self, add: np.ndarray) -> np.ndarray:
        out = np.broadcast_to(self.mean_resid, add.shape).copy()
        if self.coef.size == 0:
            return out
        out = out + self.coef[0] * add
        if self.kind == "B3b":
            mag = np.log1p(np.abs(add)).mean(axis=1, keepdims=True)
            out = out + self.coef[1] * add * mag
        return out


# ----------------------------------------------------------------- splitting


def pair_disjoint_folds(n: int, k: int, seed: int):
    return list(KFold(n_splits=k, shuffle=True, random_state=seed).split(np.arange(n)))


def _blocks_random(genes, k, rng):
    return {g: i % k for i, g in enumerate(rng.permutation(genes))}


def _blocks_spectral(genes, gene_pairs, k, rng):
    """Group genes so that many pairs fall inside one block, with balanced blocks.

    Pairs bridging two blocks can never be tested under the strict no-shared-gene
    rule, so a partition that cuts few edges leaves more pairs testable. Block
    sizes are capacity-constrained, because an unbalanced partition produces a
    fold that tests most of the cohort while training on very few pairs.

    The partition uses only the gene-pair graph, i.e. which combinations the
    screen happened to contain. It never consults y, so it cannot leak outcome
    information. Falls back to the random partition if anything is unavailable.
    """
    try:
        from sklearn.cluster import KMeans
        from sklearn.manifold import SpectralEmbedding
    except Exception:
        return _blocks_random(genes, k, rng)
    idx = {g: i for i, g in enumerate(genes)}
    n = len(genes)
    adj = np.zeros((n, n))
    for x, y in gene_pairs:
        if x in idx and y in idx:
            adj[idx[x], idx[y]] += 1.0
            adj[idx[y], idx[x]] += 1.0
    adj += 1e-6  # keep the affinity graph connected
    seed = int(rng.integers(1 << 31))
    try:
        emb = SpectralEmbedding(n_components=max(1, min(k - 1, n - 2)),
                                affinity="precomputed",
                                random_state=seed).fit_transform(adj)
        centers = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(emb).cluster_centers_
    except Exception:
        return _blocks_random(genes, k, rng)
    # capacity-constrained greedy assignment: balanced blocks, nearest centre first
    dist = np.linalg.norm(emb[:, None, :] - centers[None, :, :], axis=2)
    cap = int(np.ceil(n / k))
    order = np.dstack(np.unravel_index(np.argsort(dist, axis=None), dist.shape))[0]
    block, counts = {}, {j: 0 for j in range(k)}
    for gi, bj in order:
        g = genes[int(gi)]
        if g in block or counts[int(bj)] >= cap:
            continue
        block[g] = int(bj)
        counts[int(bj)] += 1
    for g in genes:  # safety net, should not trigger
        block.setdefault(g, int(np.argmin([counts[j] for j in range(k)])))
    return block


def gene_disjoint_folds(gene_pairs, k: int = 5, seed: int = 42,
                        method: str = "spectral", min_train: int = 6,
                        min_train_frac: float = 0.30):
    """Folds where no gene of a test pair occurs in any training pair.

    Genes are partitioned into k blocks; a pair is a test pair of block j only
    if BOTH its genes fall in block j, and the training set of that fold keeps
    only pairs with NEITHER gene in block j. Pairs bridging blocks are therefore
    never tested, so coverage is partial by construction. Report the covered
    fraction alongside any gene-disjoint result.

    A fold is DROPPED when its training set is smaller than max(min_train,
    min_train_frac * n_pairs). Without that guard an unbalanced partition can
    produce a fold that tests most of the cohort while training on a handful of
    pairs: total coverage then looks like 100% while the result is worthless.
    Counting "pairs that lie inside some block" is NOT coverage for this reason
    — only pairs in a fold that survives the guard are actually evaluated.
    """
    rng = np.random.default_rng(seed)
    gene_pairs = [tuple(p) for p in gene_pairs]
    genes = sorted({g for p in gene_pairs for g in p})
    block = (_blocks_spectral(genes, gene_pairs, k, rng) if method == "spectral"
             else _blocks_random(genes, k, rng))
    floor = max(int(min_train), int(np.ceil(min_train_frac * len(gene_pairs))))
    folds, dropped = [], []
    for j in range(k):
        te = [i for i, (x, y) in enumerate(gene_pairs) if block[x] == j and block[y] == j]
        tr = [i for i, (x, y) in enumerate(gene_pairs) if block[x] != j and block[y] != j]
        if not te:
            continue
        if len(tr) < floor:
            dropped.append({"block": j, "n_train": len(tr), "n_test": len(te)})
            continue
        folds.append((np.array(tr), np.array(te)))
    if dropped:
        import warnings as _w
        _w.warn(f"gene_disjoint_folds dropped {len(dropped)} degenerate fold(s) with "
                f"fewer than {floor} training pairs: {dropped}. On a graph with one "
                f"giant component this is expected; prefer repeated_gene_holdout().",
                stacklevel=2)
    return folds


def graph_summary(gene_pairs):
    """Connected-component structure of the gene-pair graph.

    Gene-disjoint splitting is bounded by this structure: a partition can only
    avoid cutting edges where the graph is already disconnected. If one giant
    component holds most pairs, splitting by component alone is not an option
    and blocks must be cut inside that component, which is why coverage stays
    well below 100%.
    """
    gene_pairs = [tuple(p) for p in gene_pairs]
    genes = sorted({g for p in gene_pairs for g in p})
    parent = {g: g for g in genes}

    def find(g):
        while parent[g] != g:
            parent[g] = parent[parent[g]]
            g = parent[g]
        return g

    for x, y in gene_pairs:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry
    comps = {}
    for g in genes:
        comps.setdefault(find(g), []).append(g)
    rows = []
    for members in comps.values():
        s = set(members)
        np_ = sum(1 for x, y in gene_pairs if x in s and y in s)
        rows.append({"genes": len(members), "pairs": np_,
                     "pair_fraction": np_ / len(gene_pairs) if gene_pairs else 0.0})
    rows.sort(key=lambda r: -r["pairs"])
    deg = {}
    for x, y in gene_pairs:
        deg[x] = deg.get(x, 0) + 1
        deg[y] = deg.get(y, 0) + 1
    return {"n_genes": len(genes), "n_pairs": len(gene_pairs),
            "n_components": len(rows), "components": rows[:10],
            "largest_component_pair_fraction": rows[0]["pair_fraction"] if rows else 0.0,
            "mean_degree": float(np.mean(list(deg.values()))) if deg else 0.0,
            "max_degree": int(max(deg.values())) if deg else 0}


def repeated_gene_holdout(a, b, y, gene_pairs, extra=None, seeds=range(100, 120),
                          holdout_frac=0.25, config=("B3", 1, 100.0),
                          panel_size=PANEL_SIZE, residual_mode="crossfit",
                          feature_perm=None, grid=None, inner=3):
    """Repeated random gene holdout — the preferred gene-disjoint design.

    Each run holds out a random `holdout_frac` of GENES. Test pairs are those
    with BOTH genes held out; training pairs are those with NEITHER. No gene is
    shared between train and test within a run, so every run is strictly
    gene-disjoint. Unlike a single k-block split, repeated runs accumulate
    coverage across the cohort.

    A pair can be tested in several runs, so two summaries are returned:
      - observation-weighted: pools every test instance, over-weighting pairs
        that happen to be drawn more often;
      - pair-averaged: averages SSE within a unique pair first, one vote each.
    Report the pair-averaged figure as primary.

    config fixes one (baseline, rank, alpha). Pass grid=[...] instead to run
    inner model selection inside every run, which is slower but keeps selection
    inside the resampling loop.
    """
    a, b, y = (np.asarray(v, float) for v in (a, b, y))
    gene_pairs = [tuple(map(str, p)) for p in gene_pairs]
    n = a.shape[0]
    genes = sorted({g for p in gene_pairs for g in p})
    x = build_features(a, b, extra)
    if feature_perm is not None:
        rng0 = np.random.default_rng(12345)
        idx = rng0.permutation(n)
        if feature_perm == "all":
            x = x[idx]
        elif feature_perm == "extra":
            x = np.column_stack([x[:, :N_BASE_FEATURES], x[idx, N_BASE_FEATURES:]])
        else:
            raise ValueError("feature_perm must be None, 'all' or 'extra'")
    mask = panels(a, b, panel_size)
    add = a + b

    inst = []  # (pair_index, model_sse, additive_sse, baseline_sse)
    run_rows = []
    for s in seeds:
        rng = np.random.default_rng(int(s))
        k = max(1, int(round(holdout_frac * len(genes))))
        held = set(rng.permutation(genes)[:k])
        te = np.array([i for i, (u, v) in enumerate(gene_pairs) if u in held and v in held])
        tr = np.array([i for i, (u, v) in enumerate(gene_pairs)
                       if u not in held and v not in held])
        if te.size == 0 or tr.size < 10:
            continue
        if grid is None:
            pick = tuple(config)
        else:
            scores = {}
            folds_in = list(KFold(n_splits=int(min(inner, len(tr))), shuffle=True,
                                  random_state=int(s)).split(tr))
            for cfg in grid:
                per = [float(np.mean(_sse(y, _fit_predict(a, b, y, x, mask, tr[i1], tr[i2],
                                                          *cfg, seed=int(s),
                                                          residual_mode=residual_mode),
                                          mask, tr[i2]))) for i1, i2 in folds_in]
                scores[cfg] = float(np.mean(per))
            pick = min(scores, key=scores.get)
        p = _fit_predict(a, b, y, x, mask, tr, te, *pick, seed=int(s),
                         residual_mode=residual_mode)
        pb = _fit_predict(a, b, y, x, mask, tr, te, pick[0], 0, pick[2], seed=int(s),
                          residual_mode=residual_mode)
        m_s, a_s, b_s = (_sse(y, p, mask, te), _sse(y, add[te], mask, te),
                         _sse(y, pb, mask, te))
        for j, i in enumerate(te):
            inst.append((int(i), float(m_s[j]), float(a_s[j]), float(b_s[j])))
        run_rows.append({"seed": int(s), "n_train": int(tr.size), "n_test": int(te.size),
                         "config": pick,
                         "r2_vs_additive": float(1 - m_s.sum() / (a_s.sum() + EPS)),
                         "r2_vs_baseline": float(1 - m_s.sum() / (b_s.sum() + EPS))})

    if not inst:
        raise ValueError("no run produced test pairs; lower holdout_frac or add seeds")
    arr = np.array([(i, m, ad, bs) for i, m, ad, bs in inst], float)
    idx_all, m_all, a_all, b_all = arr[:, 0].astype(int), arr[:, 1], arr[:, 2], arr[:, 3]

    uniq = np.unique(idx_all)
    pm = np.array([m_all[idx_all == u].mean() for u in uniq])
    pa_ = np.array([a_all[idx_all == u].mean() for u in uniq])
    pb_ = np.array([b_all[idx_all == u].mean() for u in uniq])

    return {
        "n_runs": len(run_rows), "n_instances": int(len(idx_all)),
        "n_unique_pairs": int(len(uniq)), "coverage": float(len(uniq) / n),
        "holdout_frac": holdout_frac,
        "observation_weighted": {
            "r2_vs_additive": float(1 - m_all.sum() / (a_all.sum() + EPS)),
            "baseline_only_r2_vs_additive": float(1 - b_all.sum() / (a_all.sum() + EPS)),
            "r2_vs_baseline": float(1 - m_all.sum() / (b_all.sum() + EPS)),
            "win_rate_vs_additive": float(np.mean(m_all < a_all))},
        "pair_averaged": {
            "r2_vs_additive": float(1 - pm.sum() / (pa_.sum() + EPS)),
            "baseline_only_r2_vs_additive": float(1 - pb_.sum() / (pa_.sum() + EPS)),
            "r2_vs_baseline": float(1 - pm.sum() / (pb_.sum() + EPS)),
            "win_rate_vs_additive": float(np.mean(pm < pa_))},
        "runs_positive_vs_additive": int(sum(r["r2_vs_additive"] > 0 for r in run_rows)),
        "runs_positive_vs_baseline": int(sum(r["r2_vs_baseline"] > 0 for r in run_rows)),
        "runs": run_rows,
        "per_pair_model_sse": pm, "per_pair_additive_sse": pa_,
        "per_pair_baseline_sse": pb_, "unique_pair_indices": uniq,
        "pooled_r2_vs_additive": float(1 - pm.sum() / (pa_.sum() + EPS)),
        "pooled_r2_vs_baseline": float(1 - pm.sum() / (pb_.sum() + EPS)),
    }


def holdout_coverage_plan(gene_pairs, holdout_frac=0.25, n_runs=20, trials=200, seed=0):
    """Expected unique-pair coverage for a repeated gene-holdout design.

    Use this to pick holdout_frac and the number of seeds BEFORE running, rather
    than reporting whatever coverage happened to come out.
    """
    rng = np.random.default_rng(seed)
    gene_pairs = [tuple(map(str, p)) for p in gene_pairs]
    genes = sorted({g for p in gene_pairs for g in p})
    k = max(1, int(round(holdout_frac * len(genes))))
    covs = []
    for _ in range(trials):
        seen = set()
        for _ in range(n_runs):
            held = set(rng.permutation(genes)[:k])
            seen |= {i for i, (u, v) in enumerate(gene_pairs) if u in held and v in held}
        covs.append(len(seen) / len(gene_pairs))
    return {"holdout_frac": holdout_frac, "n_runs": n_runs,
            "mean_coverage": float(np.mean(covs)),
            "p05": float(np.quantile(covs, 0.05)),
            "p95": float(np.quantile(covs, 0.95))}


def fold_coverage(folds, n):
    """Coverage and balance of a fold set. Report BOTH with any gene-disjoint result:
    a high pooled R-squared from one fold that trains on very few pairs is not a
    cohort-level claim."""
    tested = sorted({int(i) for _, te in folds for i in te})
    tr_sizes = [int(len(tr)) for tr, _ in folds]
    return {"n_tested": len(tested), "n_total": int(n),
            "coverage": len(tested) / n if n else 0.0,
            "train_sizes": tr_sizes,
            "test_sizes": [int(len(te)) for _, te in folds],
            "min_train_size": min(tr_sizes) if tr_sizes else 0}


# ----------------------------------------------------------------- core model


def _fit_predict(a, b, y, x, mask, tr, te, kind, rank, alpha,
                 inner_folds=5, seed=0, residual_mode="crossfit"):
    """Fit baseline + PCA-ridge residual correction on tr, predict te.

    residual_mode='crossfit': each training pair's PCA/ridge target comes from a
    baseline fitted without it. 'infold': targets come from the baseline fitted
    on the whole training split. Neither mode lets any held-out pair influence
    the fit; 'crossfit' additionally removes the baseline's in-sample shrinkage
    from the target.
    """
    add = a + b
    resid = y - add
    g = y.shape[1]

    base = Baseline(kind).fit(add[tr], resid[tr], mask[tr])
    pred_te = add[te] + base.predict(add[te])

    if residual_mode == "crossfit":
        target = np.zeros_like(resid[tr])
        nf = int(min(inner_folds, len(tr)))
        for sub_tr, sub_te in KFold(n_splits=nf, shuffle=True, random_state=seed).split(tr):
            bl = Baseline(kind).fit(add[tr][sub_tr], resid[tr][sub_tr], mask[tr][sub_tr])
            target[sub_te] = resid[tr][sub_te] - bl.predict(add[tr][sub_te])
    elif residual_mode == "infold":
        target = resid[tr] - base.predict(add[tr])
    else:
        raise ValueError("residual_mode must be 'crossfit' or 'infold'")

    if rank and rank > 0:
        k = int(min(rank, len(tr) - 1, g))
        if k >= 1:
            # PCA restores the residual mean on inverse_transform; do not add it again.
            pca = PCA(n_components=k, svd_solver="full").fit(target)
            sc = StandardScaler().fit(x[tr])
            reg = Ridge(alpha=alpha).fit(sc.transform(x[tr]), pca.transform(target))
            coef_te = np.asarray(reg.predict(sc.transform(x[te]))).reshape(len(te), k)
            pred_te = pred_te + pca.inverse_transform(coef_te)
    return pred_te


def _sse(y, pred, mask, rows):
    d = (y[rows] - pred) * mask[rows]
    return np.sum(d ** 2, axis=1)


def build_features(a, b, extra=None):
    x = feature_matrix(a, b)
    if extra is not None:
        extra = np.asarray(extra, float)
        if extra.ndim != 2 or extra.shape[0] != a.shape[0] or not np.isfinite(extra).all():
            raise ValueError("extra must be a finite (pairs, features) array")
        x = np.column_stack([x, extra])
    return x


DEFAULT_GRID = [(k, r, al) for k in ("B3", "B3b")
                for r in (1, 2, 4, 8, 16) for al in (0.1, 1, 10, 100, 1000)]


def nested_cv(a, b, y, extra=None, grid=None, outer=5, inner=3, seed=42,
              folds=None, panel_size=PANEL_SIZE, residual_mode="crossfit",
              feature_perm=None):
    """Pair-disjoint nested CV with a one-standard-error selection rule.

    feature_perm: None, 'all' or 'extra'. Shuffles rows of the feature matrix
    (the whole matrix, or only the RNA-fate block) while leaving a, b, y, the
    evaluation panels and the baselines in their true order. This keeps the
    R-squared denominator fixed, so the result is a valid null for the
    feature-to-pair correspondence.
    """
    a, b, y = (np.asarray(v, float) for v in (a, b, y))
    if not (a.shape == b.shape == y.shape) or a.ndim != 2:
        raise ValueError("a, b, y must share one 2D shape (pairs, genes)")
    n = a.shape[0]
    x = build_features(a, b, extra)

    if feature_perm is not None:
        rng = np.random.default_rng(seed)
        idx = rng.permutation(n)
        if feature_perm == "all":
            x = x[idx]
        elif feature_perm == "extra":
            if x.shape[1] <= N_BASE_FEATURES:
                raise ValueError("feature_perm='extra' needs an extra feature block")
            x = np.column_stack([x[:, :N_BASE_FEATURES], x[idx, N_BASE_FEATURES:]])
        else:
            raise ValueError("feature_perm must be None, 'all' or 'extra'")

    mask = panels(a, b, panel_size)
    add = a + b
    grid = DEFAULT_GRID if grid is None else list(grid)
    if folds is None:
        folds = pair_disjoint_folds(n, outer, seed)

    pred = np.full_like(y, np.nan)
    base_pred = np.full_like(y, np.nan)
    tested = np.zeros(n, bool)
    chosen, fold_r2 = [], []

    for tr, te in folds:
        tr = np.asarray(tr); te = np.asarray(te)
        scores = {}
        inner_folds = list(KFold(n_splits=int(min(inner, len(tr))), shuffle=True,
                                 random_state=seed).split(tr))
        for cfg in grid:
            per_fold = []
            for i_tr, i_te in inner_folds:
                p = _fit_predict(a, b, y, x, mask, tr[i_tr], tr[i_te], *cfg,
                                 seed=seed, residual_mode=residual_mode)
                per_fold.append(float(np.mean(_sse(y, p, mask, tr[i_te]))))
            sem = (float(np.std(per_fold, ddof=1) / np.sqrt(len(per_fold)))
                   if len(per_fold) > 1 else 0.0)
            scores[cfg] = (float(np.mean(per_fold)), sem)
        best = min(scores, key=lambda c: scores[c][0])
        thr = scores[best][0] + scores[best][1]
        within = [c for c in grid if scores[c][0] <= thr]
        pick = sorted(within, key=lambda c: (c[0] != "B3", c[1], -c[2]))[0]
        chosen.append({"selected": pick, "lowest_loss": best,
                       "n_within_1se": len(within), "n_grid": len(grid)})

        p = _fit_predict(a, b, y, x, mask, tr, te, *pick,
                         seed=seed, residual_mode=residual_mode)
        pred[te] = p
        # same baseline, PCA/ridge term switched off: isolates the feature-driven part
        base_pred[te] = _fit_predict(a, b, y, x, mask, tr, te, pick[0], 0, pick[2],
                                     seed=seed, residual_mode=residual_mode)
        tested[te] = True
        fold_r2.append(float(1 - _sse(y, p, mask, te).sum() /
                             (_sse(y, add[te], mask, te).sum() + EPS)))

    rows = np.where(tested)[0]
    m_sse = _sse(y, pred[rows], mask, rows)
    a_sse = _sse(y, add[rows], mask, rows)
    b_sse = _sse(y, base_pred[rows], mask, rows)
    return {
        "pooled_r2_vs_additive": float(1 - m_sse.sum() / (a_sse.sum() + EPS)),
        # increment of the PCA/ridge term over its own baseline. Use THIS as the
        # statistic for feature-permutation tests: the baseline contribution is
        # unaffected by shuffling features and would otherwise swamp the signal.
        "pooled_r2_vs_baseline": float(1 - m_sse.sum() / (b_sse.sum() + EPS)),
        "baseline_only_r2_vs_additive": float(1 - b_sse.sum() / (a_sse.sum() + EPS)),
        "win_rate_vs_additive": float(np.mean(m_sse < a_sse)),
        "win_rate_vs_baseline": float(np.mean(m_sse < b_sse)),
        "n_pairs_evaluated": int(len(rows)),
        "coverage": float(len(rows) / n),
        "outer_fold_r2": fold_r2,
        "selected_configs": chosen,
        "residual_mode": residual_mode,
        "tested_indices": rows,
        "per_pair_model_sse": m_sse,
        "per_pair_additive_sse": a_sse,
        "per_pair_baseline_sse": b_sse,
        "prediction": pred,
        "baseline_prediction": base_pred,
        "panel_mask": mask,
    }


def baseline_table(a, b, y, outer=5, seed=42, panel_size=PANEL_SIZE,
                   residual_mode="crossfit", folds=None):
    """Out-of-sample pooled R2 and win rate for each baseline, no PCA term."""
    a, b, y = (np.asarray(v, float) for v in (a, b, y))
    x = feature_matrix(a, b)
    mask = panels(a, b, panel_size)
    add = a + b
    folds = pair_disjoint_folds(a.shape[0], outer, seed) if folds is None else folds
    out = {}
    for kind in ("B0", "B1", "B2a", "B2b", "B3", "B3b"):
        pred = np.full_like(y, np.nan)
        tested = np.zeros(a.shape[0], bool)
        for tr, te in folds:
            pred[te] = _fit_predict(a, b, y, x, mask, np.asarray(tr), np.asarray(te),
                                    kind, 0, 1.0, seed=seed, residual_mode=residual_mode)
            tested[te] = True
        rows = np.where(tested)[0]
        m, ad = _sse(y, pred[rows], mask, rows), _sse(y, add[rows], mask, rows)
        out[kind] = {"pooled_r2_vs_additive": float(1 - m.sum() / (ad.sum() + EPS)),
                     "win_rate_vs_additive": float(np.mean(m < ad))}
    return out


# ----------------------------------------------------------------- controls


def oracle_projection(a, b, y, ranks=(1, 2, 4, 8, 16, 32), outer=5, seed=42,
                      panel_size=PANEL_SIZE, folds=None):
    """Representational ceiling: PCA coefficients fitted to the OBSERVED held-out
    residual. This is not predictive performance and must never be reported as such."""
    a, b, y = (np.asarray(v, float) for v in (a, b, y))
    mask = panels(a, b, panel_size)
    add = a + b
    resid = y - add
    folds = pair_disjoint_folds(a.shape[0], outer, seed) if folds is None else folds
    res = {}
    for r in ranks:
        pred = np.full_like(y, np.nan)
        tested = np.zeros(a.shape[0], bool)
        for tr, te in folds:
            tr = np.asarray(tr); te = np.asarray(te)
            k = int(min(r, len(tr) - 1, y.shape[1]))
            pca = PCA(n_components=k, svd_solver="full").fit(resid[tr])
            pred[te] = add[te] + pca.inverse_transform(pca.transform(resid[te]))
            tested[te] = True
        rows = np.where(tested)[0]
        m, ad = _sse(y, pred[rows], mask, rows), _sse(y, add[rows], mask, rows)
        res[int(r)] = float(1 - m.sum() / (ad.sum() + EPS))
    return res


def permutation_control(a, b, y, extra=None, block="all", n_perm=200, seed=42,
                        nested_each=True, grid=None, stat="pooled_r2_vs_baseline", **kw):
    """Null distribution of pooled R2 under shuffled feature-to-pair assignment.

    Only the feature matrix rows are shuffled; a, b, y, the panels and the
    baselines keep their true order, so the R-squared denominator is unchanged
    and the null is comparable to the observed value.

    nested_each=True repeats inner model selection inside every permutation,
    which is the honest null when selection is part of the procedure. Set it to
    False, with a single-configuration `grid`, only to save time, and say so.
    """
    g = grid if (nested_each or grid is None) else list(grid)[:1]
    obs = nested_cv(a, b, y, extra=extra, seed=seed, grid=g, **kw)[stat]
    null = []
    for i in range(n_perm):
        r = nested_cv(a, b, y, extra=extra, seed=seed + 1 + i, grid=g,
                      feature_perm=block, **kw)
        null.append(r[stat])
    null = np.asarray(null)
    return {"statistic": stat, "observed": float(obs), "n_perm": int(n_perm), "block": block,
            "nested_selection_each_permutation": bool(nested_each),
            "null_mean": float(null.mean()),
            "null_sd": float(null.std(ddof=1)) if n_perm > 1 else 0.0,
            "null_max": float(null.max()),
            "empirical_p": float((np.sum(null >= obs) + 1) / (n_perm + 1)),
            "null": null}


def repeat_cv(a, b, y, extra=None, seeds=(1, 2, 3, 4, 5), **kw):
    """Repeat the whole nested CV across seeds; reports spread of pooled R2."""
    vals, wins = [], []
    for s in seeds:
        r = nested_cv(a, b, y, extra=extra, seed=s, **kw)
        vals.append(r["pooled_r2_vs_additive"])
        wins.append(r["win_rate_vs_additive"])
    v = np.asarray(vals)
    return {"seeds": list(seeds), "pooled_r2": vals,
            "mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
            "min": float(v.min()), "max": float(v.max()),
            "win_rate_mean": float(np.mean(wins))}


def bootstrap_ci(result, n_boot=2000, seed=42, alpha=0.05):
    """Percentile CI for pooled R2, resampling the evaluated PAIRS.

    This propagates between-pair variability only. It does not cover cell
    sampling noise, the model-selection step, or the choice of fold partition.
    """
    rng = np.random.default_rng(seed)
    m = np.asarray(result["per_pair_model_sse"], float)
    ad = np.asarray(result["per_pair_additive_sse"], float)
    n = len(m)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        k = rng.integers(0, n, n)
        draws[i] = 1 - m[k].sum() / (ad[k].sum() + EPS)
    lo, hi = np.quantile(draws, [alpha / 2, 1 - alpha / 2])
    return {"point": float(result["pooled_r2_vs_additive"]),
            "ci_low": float(lo), "ci_high": float(hi),
            "n_boot": int(n_boot), "n_pairs": int(n),
            "note": "pair-level resampling only; excludes cell noise, selection and fold choice"}


# ----------------------------------------------------------------- demo / cli


def _clean(r):
    drop = {"prediction", "baseline_prediction", "panel_mask", "per_pair_model_sse",
            "per_pair_additive_sse", "per_pair_baseline_sse", "tested_indices", "null"}
    return {k: v for k, v in r.items() if k not in drop}


def run_demo(seed=7):
    rng = np.random.default_rng(seed)
    n, g = 131, 400
    a = rng.normal(size=(n, g)) * 0.25
    b = rng.normal(size=(n, g)) * 0.25
    x = feature_matrix(a, b)
    extra = rng.normal(size=(n, 4))
    latent = rng.normal(size=g)
    latent /= np.linalg.norm(latent)
    resid = (2.5 * extra[:, 0] + 1.5 * x[:, 0])[:, None] * latent[None, :] \
        + rng.normal(size=(n, g)) * 0.10
    y = a + b + resid

    print("SYNTHETIC DATA. Not Norman et al., not patient data, not a result.\n")
    print("baselines:", json.dumps(baseline_table(a, b, y), indent=2))
    for mode in ("crossfit", "infold"):
        r = nested_cv(a, b, y, residual_mode=mode)
        print(f"single_signature [{mode}]:", json.dumps(_clean(r), indent=2, default=str))
    ext = nested_cv(a, b, y, extra=extra)
    print("with_rna_fate:", json.dumps(_clean(ext), indent=2, default=str))
    print("bootstrap CI:", json.dumps(bootstrap_ci(ext, n_boot=500), indent=2))
    print("oracle ceiling (NOT predictive):", json.dumps(oracle_projection(a, b, y), indent=2))


def main():
    p = argparse.ArgumentParser(description="GeneRx-FATE v0.3 research scaffold")
    p.add_argument("--npz", help="NPZ with arrays a, b, y and optional extra, pairs")
    p.add_argument("--demo", action="store_true")
    p.add_argument("--gene-disjoint", action="store_true")
    p.add_argument("--split-method", default="spectral", choices=("spectral", "random"))
    p.add_argument("--blocks", type=int, default=5,
                   help="gene blocks for gene-disjoint CV; fewer blocks raise coverage "
                        "but shrink each training set")
    p.add_argument("--residual-mode", default="crossfit", choices=("crossfit", "infold"))
    p.add_argument("--permute", type=int, default=0)
    p.add_argument("--seeds", type=int, default=1, help="repeat the whole CV over N seeds")
    args = p.parse_args()
    if args.demo or not args.npz:
        return run_demo()

    d = np.load(args.npz, allow_pickle=True)
    a, b, y = d["a"], d["b"], d["y"]
    extra = d["extra"] if "extra" in d else None
    folds = None
    if args.gene_disjoint:
        if "pairs" not in d:
            raise SystemExit("--gene-disjoint needs a 'pairs' array of (geneA, geneB)")
        pairs = [tuple(map(str, q)) for q in d["pairs"]]
        print("graph:", json.dumps(graph_summary(pairs), indent=2))
        folds = gene_disjoint_folds(pairs, args.blocks, 42, method=args.split_method)
        print("fold coverage:", json.dumps(fold_coverage(folds, len(a)), indent=2))

    kw = dict(folds=folds, residual_mode=args.residual_mode)
    print("baselines:", json.dumps(baseline_table(a, b, y, folds=folds,
                                                  residual_mode=args.residual_mode), indent=2))
    r = nested_cv(a, b, y, extra=extra, **kw)
    print("model:", json.dumps(_clean(r), indent=2, default=str))
    print("bootstrap CI:", json.dumps(bootstrap_ci(r), indent=2))
    if args.seeds > 1:
        print("seed repeats:", json.dumps(repeat_cv(a, b, y, extra=extra,
                                                    seeds=tuple(range(1, args.seeds + 1)),
                                                    **kw), indent=2))
    if args.permute:
        pc = permutation_control(a, b, y, extra=extra, n_perm=args.permute, **kw)
        print("permutation control:", json.dumps(_clean(pc), indent=2))


if __name__ == "__main__":
    main()
