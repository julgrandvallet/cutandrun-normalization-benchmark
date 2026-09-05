process DIFFERENTIAL_BINDING {
    tag   "${method}|${contrast.id}"
    label 'process_medium'
    publishDir "${params.outdir}/differential", mode: 'copy'

    input:
    tuple val(method), val(contrast)
    path counts
    path samplesheet
    path scale_factors

    output:
    path "*.tsv", emit: results

    script:
    """
    differential_binding.R \\
        --counts ${counts} \\
        --samplesheet ${samplesheet} \\
        --scale-factors ${scale_factors} \\
        --method ${method} \\
        --numerator ${contrast.numerator} \\
        --denominator ${contrast.denominator} \\
        --fdr ${params.fdr} \\
        --out db.${method}.${contrast.id}.tsv
    """

    stub:
    """
    touch db.${method}.${contrast.id}.tsv
    """
}

process COMPARE_METHODS {
    label 'process_single'
    publishDir "${params.outdir}/benchmark", mode: 'copy'

    input:
    path results
    path contrasts

    output:
    path "benchmark_summary.tsv", emit: summary
    path "figures/*",             emit: figures
    path "benchmark.log",         emit: log

    script:
    """
    compare_methods.py --results ${results} --contrasts ${contrasts} \\
        --fdr ${params.fdr} --outdir . | tee benchmark.log
    """

    stub:
    """
    mkdir -p figures && touch benchmark_summary.tsv benchmark.log figures/placeholder.svg
    """
}
