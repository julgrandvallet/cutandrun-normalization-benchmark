// One consensus region set, counted identically for every sample. Every
// normalization method under test then operates on the SAME count matrix, so
// differences in the results are attributable to normalization alone.
process CONSENSUS_COUNTS {
    label 'process_medium'
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path peaks
    path bams
    path bais
    val  sample_order

    output:
    path "consensus_peaks.bed", emit: bed
    path "counts.tsv",          emit: counts

    script:
    """
    cat ${peaks} | cut -f1-3 | sort -k1,1 -k2,2n \\
      | bedtools merge -i - -d ${params.consensus_merge_distance} > consensus_peaks.bed

    printf 'chr\\tstart\\tend\\t%s\\n' "\$(echo '${sample_order.join(" ")}' | tr ' ' '\\t')" > counts.tsv
    bedtools multicov -bams ${bams} -bed consensus_peaks.bed >> counts.tsv
    """

    stub:
    """
    touch consensus_peaks.bed counts.tsv
    """
}
