process MACS2_CALLPEAK {
    tag   "${meta.sample}"
    label 'process_medium'
    publishDir "${params.outdir}/peaks", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path  control_bam

    output:
    tuple val(meta), path("*.narrowPeak"), emit: peaks

    script:
    def ctrl = control_bam.name != 'NO_CONTROL' ? "-c ${control_bam}" : ''
    """
    macs2 callpeak -t ${bam} ${ctrl} -f BAMPE -g ${params.macs_gsize} \\
        -n ${meta.sample} -q ${params.peak_qvalue} --keep-dup all --outdir .
    """

    stub:
    """
    touch ${meta.sample}_peaks.narrowPeak
    """
}
