#!/usr/bin/env bash
# End-to-end integration test on simulated data with a planted ground truth.
#
# The stub run (-profile test -stub) proves the DAG is wired correctly. It does
# not execute a single aligner or statistical model, so it cannot tell whether
# the benchmark actually measures anything. This does: real bowtie2, real
# samtools, real bedtools, real edgeR and DESeq2, on reads whose answer is
# known before they are generated.
#
# It asserts the benchmark's central claim. If spike-in normalization stops
# recovering the planted truth, or library-size normalization stops failing to,
# the benchmark has broken and this exits non-zero.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:-$(mktemp -d)}"
THREADS="${THREADS:-4}"
cd "$WORK"
echo "workdir: $WORK"

echo "== 1. simulate =="
"$ROOT/bin/simulate_reads.py" \
    --target-fasta  "$ROOT/assets/test/hg38.fa.gz" \
    --spikein-fasta "$ROOT/assets/test/dm6.fa.gz" \
    --outdir .

echo "== 2. index and align =="
mkdir -p idx bam logs/hg38 logs/dm6
for g in hg38 dm6; do
    gunzip -c "$ROOT/assets/test/$g.fa.gz" > "idx/$g.fa"
    bowtie2-build -q --threads "$THREADS" "idx/$g.fa" "idx/$g"
done
SAMPLES=(ESR1_KO_R1 ESR1_KO_R2 ESR1_WT_R1 ESR1_WT_R2)
for s in "${SAMPLES[@]}"; do
    for g in hg38 dm6; do
        bowtie2 --very-sensitive-local --no-unal --no-mixed --no-discordant \
                -I 10 -X 700 -p "$THREADS" -x "idx/$g" \
                -1 "fastq/${s}_1.fastq.gz" -2 "fastq/${s}_2.fastq.gz" \
                2> "logs/$g/$s.$g.bowtie2.log" \
          | samtools sort -@ "$THREADS" -o "bam/$s.$g.bam" -
        samtools index "bam/$s.$g.bam"
    done
done

echo "== 3. scale factors from the real bowtie2 logs =="
"$ROOT/bin/scale_factors.py" --log-dir logs --samplesheet samplesheet.csv \
    --target-genome hg38 --spikein-genomes dm6 --out scale_factors.tsv

echo "== 4. count over fixed windows =="
# Fixed windows rather than called peaks: this test is about normalization, and
# peak calling would add a second thing that could fail. Windows also keep the
# unchanged background in the matrix, which is where library-size normalization
# produces its false positives.
samtools faidx idx/hg38.fa
cut -f1,2 idx/hg38.fa.fai > hg38.genome
bedtools makewindows -g hg38.genome -w 400 > regions.bed
{ printf 'chr\tstart\tend'; printf '\t%s' "${SAMPLES[@]}"; printf '\n'; } > counts.tsv
bedtools multicov -bams "${SAMPLES[@]/#/bam/}" -bed regions.bed \
    2>/dev/null >> counts.tsv || {
    BAMS=(); for s in "${SAMPLES[@]}"; do BAMS+=("bam/$s.hg38.bam"); done
    bedtools multicov -bams "${BAMS[@]}" -bed regions.bed >> counts.tsv
}

echo "== 5. differential binding under each method =="
RESULTS=()
for m in none library_size spikein_dm6; do
    Rscript "$ROOT/bin/differential_binding.R" \
        --counts counts.tsv --samplesheet samplesheet.csv \
        --scale-factors scale_factors.tsv --method "$m" \
        --numerator ESR1_KO --denominator ESR1_WT --fdr 0.05 \
        --out "db.$m.tsv"
    RESULTS+=("db.$m.tsv")
done

echo "== 6. score =="
"$ROOT/bin/compare_methods.py" --results "${RESULTS[@]}" \
    --contrasts contrasts.csv --fdr 0.05 --outdir bench

echo "== 7. assert the benchmark still discriminates =="
python3 - bench/benchmark_summary.tsv <<'PY'
import sys, csv

rows = {r["method"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
spike, lib = rows["spikein_dm6"], rows["library_size"]

def check(ok, msg):
    print(("  PASS  " if ok else "  FAIL  ") + msg)
    return ok

sp_pur, lb_pur = float(spike["directional_purity"]), float(lib["directional_purity"])
ok = all([
    check(spike["recovers_truth"] == "True",
          f"spike-in recovers the planted truth (purity {sp_pur:.3f})"),
    check(lib["recovers_truth"] == "False",
          f"library size fails to (purity {lb_pur:.3f})"),
    check(sp_pur > 0.90, f"spike-in directional purity {sp_pur:.3f} > 0.90"),
    check(lb_pur < 0.50, f"library-size directional purity {lb_pur:.3f} < 0.50"),
    check(float(spike["median_lfc"]) < 0,
          f"spike-in bulk shift {float(spike['median_lfc']):+.3f} is negative"),
    check(float(lib["median_lfc"]) > float(spike["median_lfc"]),
          "library size understates the loss relative to spike-in"),
])
print("\nSIMULATION TEST " + ("PASSED" if ok else "FAILED"))
sys.exit(0 if ok else 1)
PY
