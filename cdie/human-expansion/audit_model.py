"""Deterministic checks and one-at-a-time sensitivity examples."""
from dataclasses import replace
from math import isclose

from compartments import MitochondrialExpression, NuclearTranscription
from model import Isoform, TransferRNA


def main():
    transcript = Isoform(1, 1, 0, False, .1, 100, 10, .9, .9, 1, 1)
    default = NuclearTranscription()
    first = replace(default, cap_fraction=.8, cleavage_poly_a_fraction=.5)
    second = replace(default, cap_fraction=.5, cleavage_poly_a_fraction=.8)
    y1 = first.outputs([transcript])
    y2 = second.outputs([transcript])
    assert isclose(y1["exportable_RNA_flux"], y2["exportable_RNA_flux"])
    assert y1["isoforms"] == y2["isoforms"]

    wt = TransferRNA(1, 1, 1, .9, 1, .2)
    variant = TransferRNA(.4, .9, .3, .6, .8, .1)
    baseline = MitochondrialExpression()
    endpoint = "mt_encoded_subunit_flux"
    y0 = baseline.outputs(100, .5, wt, variant)[endpoint]
    for parameter in ("rnase_p", "elac2", "mtpap", "free_28s"):
        lower = replace(baseline, **{parameter: getattr(baseline, parameter)*.9})
        y = lower.outputs(100, .5, wt, variant)[endpoint]
        assert 0 <= y <= y0
        print(f"{parameter}: -10% input, {(y/y0-1)*100:.2f}% output (synthetic)")
    print("Cap and poly(A) efficiencies are not individually identifiable from this endpoint.")


if __name__ == "__main__":
    main()
