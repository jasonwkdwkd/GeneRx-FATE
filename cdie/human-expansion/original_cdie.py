"""Central Dogma Interaction Equation. Illustrative, dimensionless amounts; time in hours."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

FIELDS = ("ligand", "GR", "H3K27ac", "H3K27me3", "CpG_methylation", "accessibility",
          "GRE_occupancy", "promoter_on", "pre_mRNA", "productive_mRNA", "retained_intron_mRNA",
          "charged_tRNA", "initiation_40S", "elongating_80S", "unfolded_protein",
          "folded_protein", "misfolded_protein", "mature_GFP", "free_40S", "free_60S",
          "functional_protein", "total_mRNA", "transcription_flux", "translation_flux")
DEFAULTS = dict(gre_log_affinity=0., promoter_log_bias=0., splice_logit_shift=0.,
                fold_log_rate=0., aggregate_log_rate=0., function_log_activity=0.)
ARMS = {"control": (0., 0.), "A": (1., 0.), "B": (0., .85), "AB": (1., .85)}


def model(state, ligand, writer, variant):
    L, G, ac, me, dna, acc, gre, q, u, mrna, ri, trna, c, e, un, native, bad, gfp = state
    free40 = max(65 - c - e, 0)
    free60 = max(65 - e, 0)
    productive = expit(1.386294 + variant["splice_logit_shift"] - .4 * gre - .8 * writer)
    transcribe = 18 * q
    recruit = .002 * mrna * free40
    join = .04 * c * free60
    finish = 7 * e * max(trna, 0) / (10 + max(trna, 0))
    fold = 1.1 * np.exp(variant["fold_log_rate"]) * un
    misfold = .15 * np.exp(variant["aggregate_log_rate"]) * un
    rates = np.array([
        2 * ligand - 1.7 * L, 2.3 * (L / (.3 + L) - G),
        .8 * (expit(-1.098612 + 1.8 * G + writer) - ac),
        .22 * (expit(-.619039 - .65 * G - .7 * writer) - me),
        .06 * (expit(-1.386294 - .4 * writer) - dna),
        1.15 * (expit(-.9 + 4 * ac - 2 * me - 1.2 * dna) - acc),
        2 * (acc * G / (.22 * np.exp(-variant["gre_log_affinity"]) + G) - gre),
        .9 * np.exp(2.6 * gre + 1.2 * (acc - expit(-.84)) + variant["promoter_log_bias"])
        * (1 - q) - .6 * q,
        transcribe - 1.3 * u, u * productive - .25 * mrna, u * (1 - productive) - .25 * ri,
        .9 * np.exp(-.7 * G * writer) * (60 - trna) - .05 * trna - 12 * finish,
        recruit - join - .2 * c, join - finish,
        finish - fold - misfold - .12 * un, fold - .12 * native, misfold - .12 * bad,
        .75 * (native - gfp) - .12 * gfp,
    ])
    output = np.r_[state, free40, free60, native * np.exp(variant["function_log_activity"]),
                   mrna + ri, transcribe, finish]
    return rates, output


def solve(initial, times, signals, variant, tol):
    result = solve_ivp(lambda _, y: model(y, *signals, variant)[0], (0, times[-1]),
                       initial, t_eval=times, method="LSODA", rtol=tol, atol=tol / 100)
    if not result.success:
        raise RuntimeError(result.message)
    return result.y.T


def simulate(variant, times, tol):
    initial = solve(np.zeros(18), np.array([500.]), ARMS["control"], variant, tol)[-1]
    np.testing.assert_allclose(model(initial, *ARMS["control"], variant)[0], 0, atol=1e-7)
    return np.array([
        [model(y, *signals, variant)[1] for y in solve(initial, times, signals, variant, tol)]
        for signals in ARMS.values()
    ])


def run(config, output):
    shifts = json.loads(config.read_text())
    if shifts.keys() - DEFAULTS.keys() or not all(np.isfinite(x) and abs(x) <= 5 for x in shifts.values()):
        raise ValueError("Use recognized finite synthetic shifts between -5 and 5")
    variant = {**DEFAULTS, **shifts}
    times = np.linspace(0, 24, 97)
    reference, altered = [simulate(params, times, 1e-9) for params in (DEFAULTS, variant)]
    for params, values in ((DEFAULTS, reference), (variant, altered)):
        np.testing.assert_allclose(values, simulate(params, times, 1e-11), rtol=2e-5, atol=2e-7)
        assert np.all(np.isfinite(values)) and np.min(values) >= -1e-8
        assert np.max(values[:, :, 1:8]) <= 1 + 1e-8
        assert np.max(values[:, :, 11]) <= 60
        assert np.max(values[:, :, 12] + values[:, :, 13]) <= 65
        assert np.all(values[:, :, 17] <= values[:, :, 15] + 1e-8)
        protein = values[:, :, 14:17].sum(axis=2)
        # Check total protein balance directly against the modeled completion flux.
        assert np.all(protein >= 0)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "trajectories.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["genotype", "condition", "hour", *FIELDS])
        for name, values in (("reference", reference), ("variant", altered)):
            for arm, rows in zip(ARMS, values):
                writer.writerows([name, arm, t, *row] for t, row in zip(times, rows))
    with (output / "interactions.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["hour", "endpoint", "gamma_reference", "gamma_variant", "genotype_by_A_by_B"])
        contrast = np.array([1, -1, -1, 1])
        r, v = [np.einsum("c,ctf->tf", contrast, x) for x in (reference, altered)]
        for i, t in enumerate(times):
            writer.writerows([t, name, r[i, j], v[i, j], v[i, j] - r[i, j]]
                             for j, name in enumerate(FIELDS))
    print("Verified: equilibrium, solver convergence, bounded states, resource constraints.")
    i = FIELDS.index("functional_protein")
    print("24-h AB functional protein, reference vs synthetic variant:", reference[-1, -1, i], altered[-1, -1, i])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", type=Path, default=Path(__file__).with_name("variant_example.json"))
    parser.add_argument("--output", type=Path, default=Path("simulation"))
    args = parser.parse_args()
    run(args.variant, args.output)
