process BOWTIE2_BUILD {
    tag   "${genome_name}"
    label 'process_high'
    storeDir { "${params.genome_cache}/${genome_name}" }

    input:
    tuple val(genome_name), val(fasta_url)

    output:
    tuple val(genome_name), path("index"), emit: index

    script:
    """
    mkdir -p index
    curl -fsSL ${fasta_url} | gunzip -c > index/${genome_name}.fa
    bowtie2-build --threads ${task.cpus} index/${genome_name}.fa index/${genome_name}
    """

    stub:
    """
    mkdir -p index && touch index/${genome_name}.rev.1.bt2
    """
}
