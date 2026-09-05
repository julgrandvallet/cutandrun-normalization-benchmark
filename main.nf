#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * cutandrun-normalization-benchmark
 *
 * Which normalization method recovers a known genome-wide change in CUT&RUN
 * signal? Every method under test consumes the SAME consensus regions and the
 * SAME count matrix, so any difference in the differential-binding result is
 * attributable to normalization and nothing else.
 */

include { BOWTIE2_BUILD        } from './modules/local/bowtie2_build'
include { FETCH_FASTQ          } from './modules/local/fetch_fastq'
include { BOWTIE2_ALIGN        } from './modules/local/bowtie2_align'
include { DEDUP                } from './modules/local/dedup'
include { MACS2_CALLPEAK       } from './modules/local/macs2'
include { CONSENSUS_COUNTS     } from './modules/local/consensus_counts'
include { SCALE_FACTORS        } from './modules/local/scale_factors'
include { DIFFERENTIAL_BINDING } from './modules/local/differential'
include { COMPARE_METHODS      } from './modules/local/differential'

workflow {

    samplesheet = file(params.samplesheet, checkIfExists: true)
    contrasts   = file(params.contrasts,   checkIfExists: true)

    ch_samples = Channel.fromPath(samplesheet)
        .splitCsv(header: true)
        .map { [sample: it.sample, srr: it.srr, target: it.target,
                condition: it.condition, replicate: it.replicate] }

    // One index per genome: the target plus every spike-in.
    ch_genomes = Channel.fromList(
        params.genomes.collect { name, url -> tuple(name, url) })

    BOWTIE2_BUILD(ch_genomes)
    FETCH_FASTQ(ch_samples)

    // Cross every sample with every genome. Spike-in counts therefore come from
    // the same aligner and the same parameters as the target counts, which is
    // what makes the ratio between them meaningful.
    BOWTIE2_ALIGN(FETCH_FASTQ.out.reads.combine(BOWTIE2_BUILD.out.index))

    ch_target_bam = BOWTIE2_ALIGN.out.bam
        .filter { meta, genome, bam, bai -> genome == params.target_genome }

    DEDUP(ch_target_bam)

    // The IgG library is the MACS2 control, not a sample to be tested.
    ch_igg = DEDUP.out.bam
        .filter { meta, bam, bai -> meta.target == params.igg_label }
        .map { meta, bam, bai -> bam }
        .ifEmpty { file("${projectDir}/assets/NO_CONTROL") }

    ch_signal = DEDUP.out.bam.filter { meta, bam, bai -> meta.target != params.igg_label }

    MACS2_CALLPEAK(ch_signal, ch_igg.collect())

    // Sort by sample name so the count-matrix columns and the samplesheet rows
    // stay in the same order. bedtools multicov emits columns in -bams order.
    ch_counting = ch_signal.toSortedList { a, b -> a[0].sample <=> b[0].sample }

    CONSENSUS_COUNTS(
        MACS2_CALLPEAK.out.peaks.map { meta, peak -> peak }.collect(),
        ch_counting.map { rows -> rows.collect { it[1] } },
        ch_counting.map { rows -> rows.collect { it[2] } },
        ch_counting.map { rows -> rows.collect { it[0].sample } }
    )

    SCALE_FACTORS(BOWTIE2_ALIGN.out.log.map { meta, genome, log -> log }.collect(),
                  samplesheet)

    ch_contrasts = Channel.fromPath(contrasts)
        .splitCsv(header: true)
        .map { [id: it.id, numerator: it.numerator, denominator: it.denominator] }

    DIFFERENTIAL_BINDING(
        Channel.fromList(params.methods).combine(ch_contrasts),
        CONSENSUS_COUNTS.out.counts,
        samplesheet,
        SCALE_FACTORS.out.tsv
    )

    COMPARE_METHODS(DIFFERENTIAL_BINDING.out.results.collect(), contrasts)
}

