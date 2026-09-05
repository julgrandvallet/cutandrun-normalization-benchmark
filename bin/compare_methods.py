#!/usr/bin/env python3
"""Score each normalization method against the dataset's known ground truth.

Why this dataset can be scored at all
-------------------------------------
GSE331454 contains two contrasts whose true direction is known before any
analysis is run:

  ESR1-KO vs ESR1-WT   CRISPR disruption of ESR1 removes the protein being
                       profiled, so ER binding must collapse GENOME-WIDE.
  E2 vs vehicle        Estradiol activates ER, so binding must increase
                       GENOME-WIDE.

Both are global shifts, and that is the discriminating case. Library-size
normalization forces the bulk of regions to a log fold-change of zero by
construction: it cannot represent "everything went down". A method that
reports a symmetric, zero-centred result for the ESR1 knockout is not being
conservative, it is reporting an artifact of its own assumption.

Scored per method
-----------------
median_lfc          Bulk shift. Should be strongly negative (KO) / positive (E2).
directional_purity  Of significant regions, the fraction moving the expected way.
                    ~0.5 means the method resolved a global shift into noise.
n_sig               Detection count. Read WITH purity, never alone.
engine_jaccard      edgeR/DESeq2 agreement on the significant set. Low values
                    mean the result is driven by engine internals, not signal.
recovers_truth      median_lfc has the expected sign AND purity >= 0.8.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PURITY_PASS = 0.80


def score_one(df: pd.DataFrame, expected: str, fdr: float) -> dict:
    """df: one method x one contrast, columns edger_/deseq_ lfc and fdr."""
    sign = 1.0 if expected == "up" else -1.0

    med = float(np.nanmedian(df["edger_lfc"]))
    sig = df[df["edger_fdr"] < fdr]
    n_sig = int(len(sig))
    purity = float((np.sign(sig["edger_lfc"]) == sign).mean()) if n_sig else float("nan")

    e = set(df.index[df["edger_fdr"] < fdr])
    d = set(df.index[df["deseq_fdr"] < fdr])
    jac = len(e & d) / len(e | d) if (e | d) else float("nan")

    recovers = bool(np.sign(med) == sign and n_sig > 0 and purity >= PURITY_PASS)
    return {"median_lfc": med, "n_sig": n_sig, "directional_purity": purity,
            "engine_jaccard": jac, "recovers_truth": recovers}


def score(results: pd.DataFrame, contrasts: pd.DataFrame, fdr: float) -> pd.DataFrame:
    exp = dict(zip(contrasts["numerator"] + "_vs_" + contrasts["denominator"],
                   contrasts["expected_direction"]))
    rows = []
    for (method, contrast), g in results.groupby(["method", "contrast"], sort=False):
        if contrast not in exp:
            print(f"  no expected direction for {contrast}, skipping", file=sys.stderr)
            continue
        rows.append({"method": method, "contrast": contrast,
                     "expected": exp[contrast],
                     **score_one(g.reset_index(drop=True), exp[contrast], fdr)})
    return pd.DataFrame(rows).sort_values(["contrast", "method"])


def plot(results: pd.DataFrame, summary: pd.DataFrame, outdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    for contrast, sub in summary.groupby("contrast"):
        methods = list(sub["method"])
        fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))

        r = results[results["contrast"] == contrast]
        ax[0].boxplot([r.loc[r["method"] == m, "edger_lfc"].dropna() for m in methods],
                      labels=methods, showfliers=False)
        ax[0].axhline(0, color="0.6", lw=1, ls="--")
        ax[0].set_ylabel("log2 fold change (all regions)")
        ax[0].set_title(f"Bulk shift\n{contrast}")

        ax[1].bar(methods, sub["directional_purity"])
        ax[1].axhline(PURITY_PASS, color="crimson", lw=1, ls="--")
        ax[1].set_ylim(0, 1)
        ax[1].set_ylabel("fraction of significant regions\nin the expected direction")
        ax[1].set_title("Directional purity")

        ax[2].bar(methods, sub["n_sig"])
        ax[2].set_ylabel(f"significant regions")
        ax[2].set_title("Detections")

        for a in ax:
            a.tick_params(axis="x", rotation=45)
            for lab in a.get_xticklabels():
                lab.set_ha("right")
        fig.tight_layout()
        fig.savefig(outdir / f"benchmark_{contrast}.svg")
        plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", nargs="+", type=Path, required=True)
    p.add_argument("--contrasts", type=Path, required=True)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--outdir", type=Path, default=Path("."))
    a = p.parse_args()

    results = pd.concat([pd.read_csv(f, sep="\t") for f in a.results], ignore_index=True)
    contrasts = pd.read_csv(a.contrasts)
    summary = score(results, contrasts, a.fdr)

    a.outdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(a.outdir / "benchmark_summary.tsv", sep="\t", index=False)

    print("=" * 100)
    print("NORMALIZATION BENCHMARK vs KNOWN GROUND TRUTH")
    print("=" * 100)
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    for contrast, sub in summary.groupby("contrast"):
        good = list(sub.loc[sub["recovers_truth"], "method"])
        bad = list(sub.loc[~sub["recovers_truth"], "method"])
        print(f"\n{contrast}  (expected {sub['expected'].iloc[0]})")
        print(f"  recovers ground truth : {', '.join(good) or 'NONE'}")
        print(f"  fails                 : {', '.join(bad) or 'none'}")

    plot(results, summary, a.outdir / "figures")
    print(f"\nwrote {a.outdir}/benchmark_summary.tsv and {a.outdir}/figures/")


def demo() -> None:
    """Self-check: a fabricated global-loss contrast must expose library-size norm."""
    rng = np.random.default_rng(0)
    n = 2000
    frames = []
    # Truth: every region lost signal. Spike-in preserves the shift; library-size
    # normalization re-centres it on zero, which is exactly the failure mode.
    for method, shift in [("spikein_dm6", -2.0), ("library_size", 0.0)]:
        lfc = rng.normal(shift, 0.5, n)
        fdr = np.where(np.abs(lfc - shift) < 5, 0.001, 0.9)
        frames.append(pd.DataFrame({"method": method, "contrast": "ESR1_KO_vs_ESR1_WT",
                                    "edger_lfc": lfc, "edger_fdr": fdr,
                                    "deseq_lfc": lfc, "deseq_fdr": fdr}))
    contrasts = pd.DataFrame({"numerator": ["ESR1_KO"], "denominator": ["ESR1_WT"],
                              "expected_direction": ["down"]})
    s = score(pd.concat(frames, ignore_index=True), contrasts, 0.05).set_index("method")

    assert bool(s.loc["spikein_dm6", "recovers_truth"]) is True
    assert bool(s.loc["library_size", "recovers_truth"]) is False
    # The failure is specifically loss of direction, not loss of detections.
    assert s.loc["library_size", "directional_purity"] < 0.6
    assert s.loc["library_size", "n_sig"] > 0
    print("demo OK")


if __name__ == "__main__":
    demo() if "--demo" in sys.argv else main()
