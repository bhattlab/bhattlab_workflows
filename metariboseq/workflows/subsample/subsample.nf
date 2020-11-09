#!/usr/bin/env nextflow

metagenomic = Channel.from(params.sampleSpecs)
    .map { it.metagenomic }

fastq = Channel.from(params.sampleSpecs)
    .map { it.metariboseq }
    .concat(metagenomic)

process subsample {
    echo true
    publishDir "${params.outputPrefix}", mode: "copy"

    input:
    val(name) from fastq

    output:
    path("${name}.fq.gz")

    shell:
    """
    module load seqtk
    seqtk sample !{params.inputPrefix}/!{name}.fq.gz !{params.subSamplePercent} > !{name}.fq
    gzip !{name}.fq
    """
}
