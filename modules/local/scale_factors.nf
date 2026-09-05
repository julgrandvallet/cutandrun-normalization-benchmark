process SCALE_FACTORS {
    label 'process_single'
    publishDir "${params.outdir}/normalization", mode: 'copy'

    input:
    path logs, stageAs: 'logs/*/*'
    path samplesheet

    output:
    path "scale_factors.tsv", emit: tsv
    path "scale_factors.log", emit: log

    script:
    """
    scale_factors.py \\
        --log-dir logs \\
        --samplesheet ${samplesheet} \\
        --target-genome ${params.target_genome} \\
        --spikein-genomes ${params.spikein_genomes.join(' ')} \\
        --out scale_factors.tsv | tee scale_factors.log
    """

    stub:
    """
    touch scale_factors.tsv scale_factors.log
    """
}
