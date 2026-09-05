// Duplicate handling is a normalization decision, not a QC detail. CUT&RUN
// libraries from low-abundance targets carry genuine duplicate fragments, so
// params.remove_duplicates defaults to false and the flag is exposed rather
// than buried. See docs/decisions.md.
process DEDUP {
    tag   "${meta.sample}"
    label 'process_medium'
    publishDir "${params.outdir}/align/target", mode: 'copy', pattern: '*.metrics.txt'

    input:
    tuple val(meta), val(genome_name), path(bam), path(bai)

    output:
    tuple val(meta), path("*.final.bam"), path("*.final.bam.bai"), emit: bam
    path "*.metrics.txt",                                          emit: metrics

    script:
    def prefix = "${meta.sample}"
    def action = params.remove_duplicates ? '-r' : ''
    """
    samtools collate -@ ${task.cpus} -O ${bam} \\
      | samtools fixmate -m -@ ${task.cpus} - - \\
      | samtools sort -@ ${task.cpus} - \\
      | samtools markdup ${action} -@ ${task.cpus} -f ${prefix}.markdup.metrics.txt - ${prefix}.final.bam
    samtools index ${prefix}.final.bam
    """

    stub:
    """
    touch ${meta.sample}.final.bam ${meta.sample}.final.bam.bai ${meta.sample}.markdup.metrics.txt
    """
}
