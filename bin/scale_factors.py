#!/usr/bin/env python3
"""Compute every scale factor under test from bowtie2 alignment logs.

One row per sample, one column per normalization method. Each method is a
single well-defined function of read counts, so the downstream differential
analysis differs ONLY in which column it consumes.

Methods
-------
none            1.0 for every sample. The null comparator.
library_size    target_pairs scaled to the dataset mean. What DESeq2/edgeR
                default normalization approximates in the absence of a
                composition shift.
spikein_dm6     C / dm6_pairs        (Drosophila carrier chromatin)
spikein_sacCer3 C / sacCer3_pairs    (yeast carrier DNA)
igg             target_pairs of the IgG control, scaled to the dataset mean.

Read counts are PAIRS, taken as (paired_aligned_one + paired_aligned_multi)
from the bowtie2 log. Multi-mapping pairs are included deliberately: spike-in
genomes carry repetitive sequence, and dropping multimappers removes a
condition-independent constant from the denominator, which biases the ratio
between samples rather than cancelling out of it.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

SPIKEIN_CONSTANT = 10_000  # arbitrary; only between-sample ratios matter

_PAIR_ONE = re.compile(r"^\s*(\d+) \([\d.]+%\) aligned concordantly exactly 1 time", re.M)
_PAIR_MULTI = re.compile(r"^\s*(\d+) \([\d.]+%\) aligned concordantly >1 times", re.M)


def aligned_pairs(log_text: str) -> int:
    """Concordantly aligned pairs (unique + multi) from a bowtie2 stderr log."""
    one = _PAIR_ONE.search(log_text)
    multi = _PAIR_MULTI.search(log_text)
    if one is None:
        raise ValueError("no 'aligned concordantly exactly 1 time' line in log")
    return int(one.group(1)) + (int(multi.group(1)) if multi else 0)


def parse_logs(log_dir: Path) -> pd.DataFrame:
    """log_dir/<genome>/<sample>.<genome>.bowtie2.log -> tidy counts."""
    rows = []
    for log in sorted(log_dir.rglob("*.bowtie2.log")):
        sample, genome, _, _ = log.name.split(".")
        rows.append({"sample": sample, "genome": genome,
                     "pairs": aligned_pairs(log.read_text())})
    if not rows:
        raise SystemExit(f"no bowtie2 logs under {log_dir}")
    return pd.DataFrame(rows).pivot(index="sample", columns="genome", values="pairs")


def scale_factors(counts: pd.DataFrame, meta: pd.DataFrame, target: str,
                  spikeins: list[str], igg_target: str = "IgG") -> pd.DataFrame:
    """counts: sample x genome pairs. meta: sample, target, condition."""
    m = meta.set_index("sample").loc[counts.index]
    sf = pd.DataFrame(index=counts.index)
    sf["target_pairs"] = counts[target]

    sf["none"] = 1.0
    sf["library_size"] = counts[target].mean() / counts[target]

    for sp in spikeins:
        if sp not in counts.columns:
            print(f"  spike-in genome '{sp}' absent from logs, skipping", file=sys.stderr)
            continue
        sf[f"spikein_{sp}"] = SPIKEIN_CONSTANT / counts[sp]
        sf[f"{sp}_pairs"] = counts[sp]
        # Fraction of the library that is carrier. If this varies wildly the
        # spike-in was not added proportionally and the method is unsafe.
        sf[f"{sp}_frac"] = counts[sp] / (counts[sp] + counts[target])

    igg = m.index[m["target"] == igg_target]
    if len(igg):
        # One IgG control for the dataset: every sample is referenced to it, so
        # this column is constant. It becomes informative only with per-condition
        # IgG, which this design does not have. Reported for completeness and
        # explicitly flagged as non-discriminating.
        sf["igg"] = counts.loc[igg, target].mean() / counts[target]
    return sf


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--log-dir", type=Path, required=True)
    p.add_argument("--samplesheet", type=Path, required=True)
    p.add_argument("--target-genome", default="hg38")
    p.add_argument("--spikein-genomes", nargs="+", default=["dm6", "sacCer3"])
    p.add_argument("--out", type=Path, default=Path("scale_factors.tsv"))
    a = p.parse_args()

    counts = parse_logs(a.log_dir)
    meta = pd.read_csv(a.samplesheet)
    sf = scale_factors(counts, meta, a.target_genome, a.spikein_genomes)
    sf.to_csv(a.out, sep="\t")

    print("=" * 88)
    print("SCALE FACTORS (relative between-sample scaling is what matters)")
    print("=" * 88)
    print(sf.to_string(float_format=lambda v: f"{v:.5g}"))

    for sp in a.spikein_genomes:
        col = f"{sp}_frac"
        if col in sf:
            lo, hi = sf[col].min(), sf[col].max()
            verdict = "OK" if hi / max(lo, 1e-12) < 5 else "UNSTABLE"
            print(f"\n{sp} carrier fraction {lo:.4%}-{hi:.4%}  "
                  f"max/min={hi / max(lo, 1e-12):.2f}x  [{verdict}]")
            if verdict == "UNSTABLE":
                print(f"  !! {sp} carrier fraction varies >5x across samples. The "
                      f"spike-in was not incorporated proportionally; scale factors "
                      f"derived from it are not trustworthy for this dataset.")
    print(f"\nwrote {a.out}")


def demo() -> None:
    """Self-check on hand-computable inputs."""
    log = ("100 reads; of these:\n"
           "  80 (80.00%) aligned concordantly exactly 1 time\n"
           "  15 (15.00%) aligned concordantly >1 times\n")
    assert aligned_pairs(log) == 95

    counts = pd.DataFrame({"hg38": [1000, 3000], "dm6": [100, 100]},
                          index=["a", "b"])
    meta = pd.DataFrame({"sample": ["a", "b"], "target": ["ER", "ER"],
                         "condition": ["x", "y"]})
    sf = scale_factors(counts, meta, "hg38", ["dm6"])
    # library_size: mean 2000, so a -> 2.0, b -> 0.667; ratio a/b == 3.
    assert abs(sf.loc["a", "library_size"] / sf.loc["b", "library_size"] - 3.0) < 1e-9
    # equal spike-in counts => equal spike-in scale factors, regardless of depth.
    assert abs(sf.loc["a", "spikein_dm6"] - sf.loc["b", "spikein_dm6"]) < 1e-12
    # and that is the whole point: the two methods disagree by 3x here.
    print("demo OK")


if __name__ == "__main__":
    demo() if "--demo" in sys.argv else main()
