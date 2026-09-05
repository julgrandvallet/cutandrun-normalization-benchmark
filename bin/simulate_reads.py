#!/usr/bin/env python3
"""Simulate a CUT&RUN experiment whose answer is known before it is analysed.

The benchmark scores normalization methods against a ground truth. That scoring
is itself code, and code that is never wrong on a case with a known answer has
not been tested. This generates such a case.

The physical model
------------------
Spike-in carrier is added per cell, so its ABSOLUTE amount does not depend on
the target. Sequencing then fills a fixed number of reads from whatever is in
the library. So when target signal collapses, fewer target fragments compete
for that fixed output and everything else, carrier included, is sampled more
deeply. That is the mechanism spike-in normalization relies on.

Pre-sequencing material (arbitrary units), for a true 4x loss at peaks:

    condition   peak    bg     spikein   total
    WT          30000   5000   10000     45000
    KO           7500   5000   10000     22500

Both libraries are then sequenced to the same 45000 reads, so KO is sampled at
2x. That is what actually lands in the FASTQ:

    condition   peak    bg     spikein   total   target   library_size
                                                          correction
    WT          30000   5000   10000     45000   35000    1.00
    KO          15000  10000   20000     45000   25000    1.40x up

  spikein      compares 15000/20000 against 30000/10000  ->  log2 = -2.00, exact.
  library_size compares 15000/25000 against 30000/35000  ->  log2 = -0.51.

Library-size normalization absorbs three quarters of a real fourfold loss into
its own depth correction, because it cannot tell a global biological change
from a depth difference. That is the failure this benchmark exists to measure.

If a future change breaks the scoring, or the scale-factor units, or the
engine parameterization, this case stops reproducing and CI fails.
"""
from __future__ import annotations

import argparse
import gzip
import random
import sys
from pathlib import Path

READ_LEN = 50
FRAG_MIN, FRAG_MAX = 120, 220
PEAK_WIDTH = 400
N_PEAKS = 60

# reads per sample: (peak, background, spikein). Totals are equal by design, so
# the two methods cannot be told apart by depth alone.
DESIGN = {
    "ESR1_WT": (30_000, 5_000, 10_000),
    "ESR1_KO": (15_000, 10_000, 20_000),
}
TRUE_LOG2FC = -2.0  # spike-in-normalized: (15000/20000) / (30000/10000)


def read_fasta(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">")).upper()


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def emit_pairs(seq: str, starts, rng: random.Random, r1, r2, tag: str, n0: int) -> int:
    """Write one read pair per fragment start. Returns the next read index."""
    n = n0
    for s in starts:
        frag = rng.randint(FRAG_MIN, FRAG_MAX)
        s = max(0, min(s, len(seq) - frag - 1))
        f = seq[s:s + frag]
        if len(f) < READ_LEN or f.count("N") > len(f) // 4:
            continue
        a, b = f[:READ_LEN], revcomp(f[-READ_LEN:])
        q = "I" * READ_LEN
        r1.write(f"@{tag}_{n}/1\n{a}\n+\n{q}\n")
        r2.write(f"@{tag}_{n}/2\n{b}\n+\n{q}\n")
        n += 1
    return n


def simulate(target: str, spikein: str, peaks, counts, rng: random.Random,
             out1: Path, out2: Path, tag: str) -> None:
    n_peak, n_bg, n_spike = counts
    with gzip.open(out1, "wt") as r1, gzip.open(out2, "wt") as r2:
        # Enriched fragments, concentrated in the peak regions.
        starts = [rng.choice(peaks) + rng.randint(0, PEAK_WIDTH - FRAG_MAX)
                  for _ in range(n_peak)]
        n = emit_pairs(target, starts, rng, r1, r2, f"{tag}_peak", 0)
        # Uniform background across the target genome.
        starts = [rng.randrange(0, len(target) - FRAG_MAX) for _ in range(n_bg)]
        n = emit_pairs(target, starts, rng, r1, r2, f"{tag}_bg", n)
        # Carrier genome.
        starts = [rng.randrange(0, len(spikein) - FRAG_MAX) for _ in range(n_spike)]
        emit_pairs(spikein, starts, rng, r1, r2, f"{tag}_spike", 0)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--target-fasta", type=Path, required=True)
    p.add_argument("--spikein-fasta", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--replicates", type=int, default=2)
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()

    rng = random.Random(a.seed)
    target, spikein = read_fasta(a.target_fasta), read_fasta(a.spikein_fasta)
    fq = a.outdir / "fastq"
    fq.mkdir(parents=True, exist_ok=True)

    # Shared peak positions: both conditions bind the SAME sites, only the
    # amount differs. A benchmark of normalization must not confound a
    # magnitude change with a location change.
    peaks = sorted(rng.randrange(1000, len(target) - PEAK_WIDTH - 1000)
                   for _ in range(N_PEAKS))

    rows = ["sample,fastq_1,fastq_2,target,condition,replicate"]
    for cond, counts in DESIGN.items():
        for rep in range(1, a.replicates + 1):
            s = f"{cond}_R{rep}"
            f1, f2 = fq / f"{s}_1.fastq.gz", fq / f"{s}_2.fastq.gz"
            simulate(target, spikein, peaks, counts,
                     random.Random(a.seed + hash(s) % 10_000), f1, f2, s)
            rows.append(f"{s},{f1.resolve()},{f2.resolve()},ER,{cond},{rep}")
            print(f"  {s}: peak={counts[0]} bg={counts[1]} spikein={counts[2]}")

    (a.outdir / "samplesheet.csv").write_text("\n".join(rows) + "\n")
    (a.outdir / "contrasts.csv").write_text(
        "id,numerator,denominator,expected_direction,rationale\n"
        "sim_ko_vs_wt,ESR1_KO,ESR1_WT,down,"
        f"simulated {TRUE_LOG2FC} log2 loss at every peak with constant "
        "absolute spike-in\n")
    print(f"\n{N_PEAKS} shared peaks; true log2FC at peaks = {TRUE_LOG2FC}")
    print(f"wrote {a.outdir}/samplesheet.csv and contrasts.csv")


def demo() -> None:
    """Self-check: the design must actually separate the two methods.

    This asserts the arithmetic of the simulated design, before any reads are
    written. If DESIGN is edited into something that no longer distinguishes
    the methods, or no longer encodes TRUE_LOG2FC, this fails immediately
    rather than producing a benchmark that silently proves nothing.
    """
    import math

    wt_peak, wt_bg, wt_spike = DESIGN["ESR1_WT"]
    ko_peak, ko_bg, ko_spike = DESIGN["ESR1_KO"]
    wt_target, ko_target = wt_peak + wt_bg, ko_peak + ko_bg

    assert sum(DESIGN["ESR1_WT"]) == sum(DESIGN["ESR1_KO"]), \
        "equal sequencing totals, or depth alone distinguishes the methods"

    # Spike-in: correct by carrier, which the target cannot influence.
    spikein_lfc = math.log2((ko_peak / ko_spike) / (wt_peak / wt_spike))
    # Library size: correct by target depth, which the biology just changed.
    libsize_lfc = math.log2((ko_peak / ko_target) / (wt_peak / wt_target))

    assert abs(spikein_lfc - TRUE_LOG2FC) < 1e-9, \
        f"spike-in must recover the truth exactly, got {spikein_lfc:.3f}"
    assert abs(libsize_lfc - TRUE_LOG2FC) > 1.0, \
        f"library-size error must be large enough to detect, got {libsize_lfc:.3f}"
    assert libsize_lfc > spikein_lfc, \
        "library size must UNDERSTATE the loss; that is the failure mode"

    print(f"demo OK  truth={TRUE_LOG2FC:.2f}  spikein={spikein_lfc:.2f}  "
          f"library_size={libsize_lfc:.2f}  "
          f"(library size hides {100 * (1 - libsize_lfc / TRUE_LOG2FC):.0f}% of the loss)")


if __name__ == "__main__":
    demo() if "--demo" in sys.argv else main()
