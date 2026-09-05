// Aligns one sample against one genome. Called once per genome (target +
// each spike-in) so the spike-in read counts come from the same aligner and
// the same settings as the target counts.
process BOWTIE2_ALIGN {
    tag   "${meta.sample}|${genome_name}"
    label 'process_high'
    publishDir { "${params.outdir}/align/${genome_name}" }, mode: 'copy', pattern: '*.log'

    input:
    tuple val(meta), path(reads), val(genome_name), path(index_dir)

    output:
    tuple val(meta), val(genome_name), path("*.bam"), path("*.bam.bai"), emit: bam
    tuple val(meta), val(genome_name), path("*.log"),                    emit: log

    script:
    def prefix = "${meta.sample}.${genome_name}"
    // --no-overlap/--no-dovetail are deliberately OFF: CUT&RUN fragments are
    // short and frequently dovetail. Excluding them biases the fragment-size
    // distribution and, for spike-ins, the scale factor itself.
    """
    IDX=\$(find -L ${index_dir} -name '*.rev.1.bt2' | head -n1 | sed 's/\\.rev\\.1\\.bt2\$//')
    bowtie2 --very-sensitive-local --no-unal --no-mixed --no-discordant \\
            -I 10 -X 700 -p ${task.cpus} -x \$IDX \\
            -1 ${reads[0]} -2 ${reads[1]} 2> ${prefix}.bowtie2.log \\
      | samtools sort -@ ${task.cpus} -o ${prefix}.bam -
    samtools index ${prefix}.bam
    """

    stub:
    """
    touch ${meta.sample}.${genome_name}.bam ${meta.sample}.${genome_name}.bam.bai
    echo "0 reads; of these:" > ${meta.sample}.${genome_name}.bowtie2.log
    """
}
