#!/usr/bin/env nextflow

// Create two channels, one each for the metagenomic and metariboseq
// samples, in order to label each sample file appropriately.
metagenomicTrim = Channel.from(params.sampleSpecs)
    .map { [name:it.name, type:"metagenomic", file:it.metagenomic] }

metariboseqTrim = Channel.from(params.sampleSpecs)
    .map { [name:it.name, type:"metariboseq", file:it.metariboseq] }

// Concatenate the labeled samples into a single channel for trimming.
// For n samples this contain 2*n entries, one for metagenomic and one
// for metariboseq data.
toTrim = metagenomicTrim.concat(metariboseqTrim)

// Trim adapters from all samples and save the trim report to the results
// directory.
process adapterTrim  {
    publishDir "${params.resultsPrefix}", pattern: "*trimming_report.txt", mode: "copy"

    input:
    val sample from toTrim

    output:
        tuple val("${sample.name}"), sample, path("${sample.file}_trimmed.fq.gz") into trimmedAll
        path("${sample.file}.fq.gz_trimming_report.txt")

    shell:
    '''
    trim_galore  !{params.inputPrefix}/!{sample.file}.fq.gz
    '''
}

// Create a channel with all samples that are to be aligned and a separate
// channel with the samples to be assembled. Only the metagenomic samples
// are used to generate assemblies.
trimmedAll.into{ alignmentCandidates; assemblyCandidates}
toAssemble = assemblyCandidates.filter { it[1].type == "metagenomic" }

// Create an assembly from the metagenomic reads that will be used for
// both metagenomic and metariboseq alignments.
process metagenomicAssembly {
    memory "${params.assemblyMemory}"
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
    tuple name, sample, trimmed from toAssemble

    output:
    tuple name, sample, path("${sample.name}-assembly") into assembled
    path("${sample.name}-assembly/contigs.fasta")

    shell:
    '''
    spades.py !{params.spadesOptions} -o !{sample.name}-assembly -s !{trimmed}
    '''
}

process bowtieIndex {
    publishDir "${params.resultsPrefix}/${sample.name}-assembly/", mode: "copy"

    input:
    tuple name, sample, assembly from assembled

    output:
    tuple name, sample, assembly, path("contigs.*") into assembledAndIndexed
    path("contigs.*")

    shell:
    '''
    bowtie-build !{params.bowtieIndexOptions} !{assembly}/contigs.fasta contigs
    '''
}

// Combine the assembly information with every sample to be aligned using the
// name of the sample as the key. The resulting channel will contain 2*n records
// for n samples, one record for every metagenomic and metariboseq data file.
// Each of those records is paired with the corresponding assembly by the
// comebine operator.
toAlign = alignmentCandidates.combine(assembledAndIndexed, by: [0]).map {
    // it[1] is the sample object.
    // it[2] the name of the trimmed output file.
    // it[3] is the sample object again, which is dropped.
    // it[4] is the name of the directory containing the assembly.
    // it[5] the bowtie index for the assembly.
    [ it[1], it[2], it[4], it[5] ]
}

process alignment {
    memory "${params.alignmentMemory}"
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
    tuple sample, trimmed, assembly, path(index) from toAlign

    output:
    path("${sample.file}.bam")

    shell:
    '''
    bowtie !{params.bowtieAlignmentOptions} --sam -l 20 contigs <(zcat !{trimmed}) !{sample.file}.sam
    samtools view -b -F 4 !{sample.file}.sam | samtools sort - > !{sample.file}.bam
    '''
}
