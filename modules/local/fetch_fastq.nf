process FETCH_FASTQ {
    tag   "${meta.sample}"
    label 'process_medium'
    storeDir "${params.outdir}/fastq"

    input:
    val meta

    output:
    tuple val(meta), path("${meta.sample}_{1,2}.fastq.gz"), emit: reads

    script:
    """
    prefetch --max-size 100G ${meta.srr}
    fasterq-dump --split-files --threads ${task.cpus} --outdir . ${meta.srr}
    pigz -p ${task.cpus} -c ${meta.srr}_1.fastq > ${meta.sample}_1.fastq.gz
    pigz -p ${task.cpus} -c ${meta.srr}_2.fastq > ${meta.sample}_2.fastq.gz
    rm -f ${meta.srr}_*.fastq
    """

    stub:
    """
    printf '' | gzip > ${meta.sample}_1.fastq.gz
    printf '' | gzip > ${meta.sample}_2.fastq.gz
    """
}
