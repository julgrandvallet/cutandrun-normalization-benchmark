# CUT&RUN normalization benchmark

**Which normalization method recovers a genome-wide change in CUT&RUN signal, and which one hides it?**

[![ci](https://github.com/julgrandvallet/cutandrun-normalization-benchmark/actions/workflows/ci.yml/badge.svg)](https://github.com/julgrandvallet/cutandrun-normalization-benchmark/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.10-0DC09D)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

> **Status:** validated end to end on simulated data with a planted ground truth — real bowtie2,
> samtools, bedtools, edgeR and DESeq2, asserted in CI on every commit (see
> [Validation](#validation-the-benchmark-tests-itself)). The run on the real dataset, GSE331454,
> is in progress; [Results](#results-gse331454) will be filled in from that run and not before.
> Every number on this page was produced by this pipeline.

---

## The problem

Differential binding analysis assumes most regions do not change. Library-size normalization
is built on that assumption: it scales samples so the bulk of regions sits at a log fold-change
of zero.

When the assumption holds, this is fine. When a treatment changes signal **genome-wide**, it is
not conservative, it is wrong in a specific and dangerous way: the method cannot represent
"everything went down", so it reports the global shift as noise and returns a symmetric,
zero-centred result that looks perfectly reasonable. Spike-in normalization exists to break
that assumption by referencing an exogenous genome that the treatment cannot affect.

The open question is not whether spike-in normalization is theoretically better. It is how much
it matters in practice, on real CUT&RUN data, and whether the answer survives changing the
statistical engine.

## Why this dataset can be scored

Most normalization comparisons have no ground truth, so they report concordance between methods
and stop. [GSE331454](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE331454) has one built in.

ER CUT&RUN in MCF7 cells, with **two independent spike-ins** (*Drosophila* dm6 and yeast sacCer3),
across four conditions:

| Contrast | Manipulation | Known truth | Why it discriminates |
|---|---|---|---|
| `ESR1_KO` vs `ESR1_WT` | CRISPR disruption of *ESR1* | ER binding collapses | The profiled protein is gone. A global loss. |
| `estradiol` vs `charcoal_stripped` | E2 activates ER | ER binding increases | A global gain, opposite direction. |

Both are global shifts, in opposite directions. That is exactly the regime where library-size
normalization is expected to fail, and having both directions rules out a method that simply
biases everything one way.

Two spike-in species in the same libraries also allow the spike-ins to be checked against
*each other*, which is a control on the control.

## What is measured

Every method consumes the **same consensus regions** and the **same count matrix**. The only thing
that varies is the scale factor. Any difference in the result is therefore attributable to
normalization and nothing else.

| Method | Scale factor | Assumption |
|---|---|---|
| `none` | 1.0 | none; the null comparator |
| `library_size` | mean target pairs / target pairs | most regions do not change |
| `spikein_dm6` | C / dm6 pairs | carrier chromatin is added proportionally |
| `spikein_sacCer3` | C / sacCer3 pairs | as above, different species |

Each is scored against the known truth:

- **`median_lfc`** — bulk shift. Must carry the expected sign.
- **`directional_purity`** — of significant regions, the fraction moving the expected way.
  ~0.5 means the method dissolved a global shift into noise.
- **`n_sig`** — detections. Meaningless without purity; a method can be made arbitrarily
  sensitive and arbitrarily wrong at the same time.
- **`engine_jaccard`** — edgeR/DESeq2 agreement. Low values mean the result is driven by
  engine internals rather than by signal.
- **`recovers_truth`** — correct sign **and** purity ≥ 0.8.

## The parameterization trap

The two statistical engines take an external scale factor in **different units**, and the
difference is silent. Pass the same numbers to both and you get results that disagree by tens of
percent, with no warning, no error, and entirely plausible-looking output.

```
edgeR    effective library size = lib.size * norm.factors
         -> norm.factors is a CORRECTION to depth: (1/s) / lib.size, centred so the product is 1

DESeq2   normalized counts = counts / sizeFactors
         -> sizeFactors is the COMPLETE divisor: 1/s, depth already included
```

Multiplying a DESeq2 size factor by library size, as if it were an edgeR norm factor,
double-counts sequencing depth. In our hands this manufactured a **~40% apparent discrepancy
between two engines that in fact agreed** — a discrepancy that looks exactly like a real
biological finding until you check the units.

[`bin/differential_binding.R`](bin/differential_binding.R) hands both engines the same `s`, then
asserts that their effective per-sample scaling agrees to within 2% after removing the global
constant. If the parameterization ever drifts back, the pipeline stops rather than reporting
a wrong number.

## Quick start

```bash
# Full run: downloads 9 SRA runs (~9.6 GB) and builds three bowtie2 indices
nextflow run julgrandvallet/cutandrun-normalization-benchmark -profile docker

# On a SLURM cluster with Singularity
NXF_SLURM_QUEUE=your_queue \
nextflow run julgrandvallet/cutandrun-normalization-benchmark -profile slurm

# Validate the whole DAG in ~30 seconds, no data required
nextflow run julgrandvallet/cutandrun-normalization-benchmark -profile test -stub

# Prove the benchmark actually discriminates, on simulated data (~2 min, needs the toolchain)
tests/simulation_test.sh
```

Your own data: point `--samplesheet` and `--contrasts` at your own CSVs
(see [`assets/`](assets/)). `expected_direction` in the contrasts file is what makes a contrast
scoreable; omit it and the contrast is analysed but not graded.

## Validation: the benchmark tests itself

A benchmark that scores methods is itself code, and code that has never been run on a case with
a known answer has not been tested. [`bin/simulate_reads.py`](bin/simulate_reads.py) generates
that case.

Spike-in carrier is added per cell, so its absolute amount does not depend on the target.
Sequencing then fills a fixed number of reads from whatever is in the library. When target
signal collapses, fewer target fragments compete for that fixed output, so the carrier is
sampled more deeply. Simulating a true **fourfold loss** at 60 shared peak regions:

```
              pre-sequencing              after sequencing to equal depth
condition   peak    bg    spike      peak    bg     spike    target
WT         30000  5000    10000     30000  5000     10000     35000
KO          7500  5000    10000     15000 10000     20000     25000
```

`tests/simulation_test.sh` then runs the real toolchain over those reads. Result:

| method | median log2FC | n_sig | directional purity | recovers truth |
|---|---|---|---|---|
| `spikein_dm6` | −0.10 | 104 | **0.98** | **yes** |
| `library_size` | **+1.38** | 480 | **0.16** | no |
| `none` | +0.89 | 405 | 0.21 | no |

**Library-size normalization reports a fourfold loss as a significant gain at 480 regions.**
Not a loss of power, not a conservative result: the wrong sign, with high confidence. It cannot
distinguish "this sample has less signal" from "this sample was sequenced less deeply", so it
corrects away the biology and then, because the unchanged background is now scaled up, reports
that background as newly enriched.

Spike-in normalization recovers the planted truth at 98% directional purity.

CI asserts all of this on every commit. If the benchmark ever stops discriminating between the
two methods, the build fails.

## Results: GSE331454

Pending the full run on the real dataset. This section will report the scored table from
`results/benchmark/benchmark_summary.tsv` and the figures in `results/benchmark/figures/`.

The simulation above establishes that the machinery detects the effect when it is present. It
does not establish the size of the effect in real data, which is the point of the full run.

## Layout

```
main.nf                        workflow
modules/local/                 one process per file
bin/scale_factors.py           bowtie2 logs -> every scale factor under test
bin/differential_binding.R     edgeR + DESeq2 from one shared count matrix
bin/compare_methods.py         scoring against ground truth, figures
bin/simulate_reads.py          reads with a planted ground truth, for validation
tests/simulation_test.sh       end-to-end integration test, asserted in CI
assets/                        samplesheet, contrasts, tiny test genomes
docs/decisions.md              analysis choices and why
```

All three Python scripts carry `--demo` self-checks that run in seconds and are executed by CI and
by the Docker build, so a broken assumption fails at build time rather than six hours into a run.
`simulate_reads.py --demo` asserts the simulated design still separates the methods before a
single read is written.

## Design decisions

Recorded in [`docs/decisions.md`](docs/decisions.md). The ones most likely to affect a result:
duplicates are **kept** by default, multi-mapping pairs are **counted** in spike-in totals, and
region filtering happens on **raw** counts before any normalization is applied.

## References

- Benchmarking normalisation methods for differential binding analysis in CUT&RUN.
  [GSE331453](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE331453) (BRG1),
  [GSE331454](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE331454) (ER).
- "The Wild West of Spike-in Normalization": benchmarking ChIP-seq spike-in normalization
  methods. [GSE273915](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273915).
- Meers MP, Bryson TD, Henikoff JG, Henikoff S. Improved CUT&RUN chromatin profiling tools.
  *eLife* 2019;8:e46314.
- Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis
  of RNA-seq data. *Genome Biology* 2010;11:R25.

## Author

Julian Grandvallet-Contreras — bioinformatics analyst, cancer genomics and epigenomics,
University of Colorado Anschutz Medical Campus.
[GitHub](https://github.com/julgrandvallet) ·
[ORCID](https://orcid.org/0000-0001-8021-070X) ·
[Scholar](https://scholar.google.com/citations?user=WKPGteEAAAAJ)

This repository uses only public data. It reimplements, on a public dataset, a normalization
question that arose in unpublished institutional work; no institutional data, results, or
sample metadata are included.
