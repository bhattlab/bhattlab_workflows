#!/usr/bin/env nextflow

workflow {
    // Dedup inputs since they can be repeated for multiple samples and the
    // metagenomic assembly is by far the most expensive step so it makes
    // sense to avoid reprocessing the same data.
    metariboseq_files = Channel.from(params.sampleSpecs)
        .map { [sample_type: 'metariboseq', path: it.metariboseq] }
        .unique()

    metagenomic_files = Channel.from(params.sampleSpecs)
        .map { [sample_type: 'metagenomic', path: it.metagenomic] }
        .unique()

    metariboseq_samples = Channel.from(params.sampleSpecs)
        .map { [it.metariboseq, it.metagenomic] }

    main:
        adapterTrim(metariboseq_files.concat(metagenomic_files))
        trimmed = adapterTrim.out.trimmed.branch{
            metagenomic: it[0] == 'metagenomic'
            metariboseq: it[0] == 'metariboseq'
        }

        assemblies = metagenomicAssembly(trimmed.metagenomic
            .map{[it[1], it[2]]}) // drop sample_type
        indexed = bowtieIndex(assemblies)

        metariboseq = trimmed.metariboseq
            .map{[it[1], it[2]]} // drop sample_type

        // combine the metagenomic assemblies with the metariboseq samples
        // that are to be aligned against them.
        combined = metariboseq_samples.combine(metariboseq, by: 0)
            .map{ [it[1], it[0], it[2]]} // move the sample name to index 0
            .combine(indexed, by: 0)
            .map{ [it[1], it[2], it[3]]} // sample name, trimmed riboseq, mag

        combined.view()

        alignmtMetariboseqAgainstAssemblies(combined).view()
}

// Trim adapters from all samples and save the trim report to the results
// directory.
process adapterTrim {
    cpus "${params.trimGaloreCPUs}"
    publishDir "${params.resultsPrefix}", pattern: "*trimming_report.txt", mode: "copy"

    input:
    tuple val(type), val(path)

    output:
    tuple val(type), val(path), path("${path}_trimmed.fq.gz"), emit: trimmed
    // The report is not needed after this process has run, simply publish it.
    path("${path}.fq.gz_trimming_report.txt"), emit: report 

    shell:
    '''
    trim_galore !{params.trimGaloreOptions} !{params.inputPrefix}/!{path}.fq.gz
    '''

    stub:
    """
    touch "${path}_trimmed.fq.gz" "${path}.fq.gz_trimming_report.txt"
    """
}

// Create an assembly from the metagenomic reads that will be used for
// both metagenomic and metariboseq alignments. Build a bowtie index for
// the MAG.
process metagenomicAssembly {
    cpus "${params.assemblyCPUs}" 
    memory "${params.assemblyMemory}"
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
    tuple val(name), path(trimmed)

    output:
    tuple val(name), path("${name}-assembly")

    shell:
    '''
    spades.py !{params.spadesOptions} -o !{name}-assembly -s !{trimmed}
    '''

    stub:
    """
    mkdir "${name}-assembly"
    touch "${name}-assembly/contigs.fasta"
    """

}

// Build a bowtie index for every assembly.
process bowtieIndex {
    cpus "${params.alignmentCPUs}"
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
    tuple val(name), path(assembly)

    output:
    tuple val(name), path(assembly)
  
    shell:
    '''
    bowtie-build !{params.bowtieIndexOptions} !{assembly}/contigs.fasta !{assembly}/contigs
    '''

    stub:
    """
    touch "${assembly}/contigs.1.ebwt"
    """
}


// Alignment every input file, metagenomic and metariboseq against the
// assembly.
process alignmtMetariboseqAgainstAssemblies {
    memory "${params.alignmentMemory}"
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
    tuple val(name), path(trimmed), path(assembly)

    output:
    path("${name}.bam")

    shell:
    '''
    bowtie !{params.bowtieAlignmentOptions} --sam -l 20 !{assembly}/contigs <(zcat !{trimmed}) !{name}.sam
    samtools view -b -F 4 !{name}.sam | samtools sort - > !{name}.bam
    '''

    stub:
    """
    touch "${name}.bam"
    """

}
