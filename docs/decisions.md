# Analysis decisions

Choices that change results, and the reasoning. Recorded because a normalization benchmark is
only interpretable if the things around normalization are held fixed and stated.

## Duplicates are kept by default

`params.remove_duplicates = false`.

CUT&RUN targets are cleaved at defined positions, so identical fragments arise for biological
reasons rather than only from PCR. For a low-abundance target, deduplication removes real signal,
and it removes it *unevenly*: the sample with more true binding loses more reads. That converts a
biological difference into an apparent depth difference, which is precisely the quantity
normalization is trying to estimate.

The flag is exposed rather than hidden so the choice can be reversed and the effect measured.

## Multi-mapping pairs are counted in spike-in totals

`aligned_pairs = paired_aligned_one + paired_aligned_multi`.

Spike-in genomes carry repetitive sequence. Discarding multimappers removes a roughly
condition-independent fraction of the spike-in, which does not cancel out of the between-sample
ratio; it shifts it. Since only relative scaling matters, including multimappers is both simpler
and less biased.

## Spike-in reads are aligned with the same aligner and settings as target reads

Each sample is aligned separately against the target genome and against each spike-in genome
using identical `bowtie2` parameters. A scale factor is a ratio of two read counts; if the two
counts come from different alignment regimes, the ratio carries the difference between the
regimes as well as the difference between the samples.

Reads that map to both the target and a spike-in genome are counted in both. The genomes here are
distant enough that this is a small effect, but it is an approximation, not an exact partition.

## Dovetailed fragments are retained

`--no-overlap` and `--no-dovetail` are deliberately not set. CUT&RUN fragments are short and
frequently dovetail. Excluding them truncates the fragment-size distribution, and for spike-ins it
biases the scale factor itself.

## Region filtering happens on raw counts, before normalization

`rowSums(counts) >= 10 * n_samples`, applied once, before any scale factor is used.

Filtering after normalization would let each method under test select its own region set, and the
methods would then no longer be compared on the same data. This is the single most important
control in the benchmark.

## One consensus region set for every method

Peaks are called per sample, merged once, and counted once. The resulting matrix is shared by
every method and both engines. Differences in the output are therefore attributable to
normalization alone.

## The IgG library is a control, not a sample

The single IgG library is used as the MACS2 control and is excluded from differential testing.
It is reported in the scale-factor table for completeness, but with one IgG library for the whole
dataset the `igg` column is constant across samples and therefore cannot discriminate between
conditions. It becomes a usable method only with per-condition IgG, which this design does not
have. It is reported and flagged rather than silently dropped.

## Carrier fraction is checked, not assumed

`bin/scale_factors.py` reports the spike-in fraction of each library and flags the method as
`UNSTABLE` if it varies more than fivefold across samples. Spike-in normalization assumes the
carrier was added in proportion to cell number; when that fails, the resulting scale factors are
confidently wrong, and nothing downstream can detect it.
