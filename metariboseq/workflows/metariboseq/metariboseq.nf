#!/usr/bin/env nextflow

workflow metagenomic_trim {
    take:
        toTrim
    main:
        adapterTrim(toTrim)
    emit:
        adapterTrim.out.trimmed
}

workflow metariboseq_trim {
    take:
        toTrim
    main:
        adapterTrim(toTrim)
    emit:
        adapterTrim.out.trimmed
}

workflow {
    metariboseqTrim = Channel.from(params.sampleSpecs)
        .map { [key: it.name, path: it.metariboseq, metagenomic_key: it.metagenomic] }

    // Dedup metagenomic inputs since they can be repeated for multiple samples
    // and the metagenomic assembly is by far the most expensive step so it makes
    // sense to avoid reprocessing the same data.
    // Note that the key field is the same as the file for the metagenomic
    // data at this stage.
    metagenomicTrim = Channel.from(params.sampleSpecs)
        .map { [key: it.metagenomic, path: it.metagenomic, metagenomic_key: it.metagenomic] }
        .unique()

    main:
        metagenomic = metagenomic_trim(metagenomicTrim)
        metariboseq = metariboseq_trim(metariboseqTrim)
        .map { [it[0], it[1], it[2]]} // key, trimmed, mag

        assemblies = metagenomicAssembly(metagenomic)
        indexed = bowtieIndex(assemblies)
             .map { [it[0], it[1], it[2]] } // key, assembly, mag

        // combine the metagenomic assemblies with the metariboseq samples
        // that are to aligned against them.
        combined = metariboseq.combine(indexed, by: 2).map{
            {[key: it[1], trimmed: it[2], contigs: it[4], metagenomic: it[0]]}
        }
        combined.view()
        alignmtMetariboseqAgainstAssemblies(combined)
}

// Trim adapters from all samples and save the trim report to the results
// directory.
process adapterTrim {
   publishDir "${params.resultsPrefix}", pattern: "*trimming_report.txt", mode: "copy"

    input:
    tuple val(key), val(trimpath), val(metagenomic)

    output:
    tuple val(key), path("${trimpath}_trimmed.fq.gz"), val(metagenomic), emit: trimmed
    // The report is not needed after this process has run, simply publish it.
    path("${trimpath}.fq.gz_trimming_report.txt"), emit: report 

    shell:
    '''
    trim_galore !{params.inputPrefix}/!{trimpath}.fq.gz
    '''
}

// Create an assembly from the metagenomic reads that will be used for
// both metagenomic and metariboseq alignments. Build a bowtie index for
// the MAG.
process metagenomicAssembly {
    memory "${params.assemblyMemory}"
    cpus "${params.assemblyCPUs}"
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
    tuple val(key), path(trimmed), val(metagenomic)

    output:
    tuple val(key), path("${key}-assembly"), val(metagenomic)

    shell:
    '''
    spades.py !{params.spadesOptions} -o !{key}-assembly -s !{trimmed}
    '''
}

// Build a bowtie index for every assembly.
process bowtieIndex {
    cpus "${params.assemblyCPUs}"
    publishDir "${params.resultsPrefix}/${key}-assembly/", mode: "copy"

    input:
    tuple val(key), path(assembly), val(metagenomic)

    output:
    tuple val(key), path(assembly), val(metagenomic)
  
    shell:
    '''
    bowtie-build !{params.bowtieIndexOptions} !{assembly}/contigs.fasta !{assembly}/contigs
    '''
}


// Alignment every input file, metagenomic and metariboseq against the
// assembly.
process alignmtMetariboseqAgainstAssemblies {
    memory "${params.alignmentMemory}"
    cpus "${params.alignmentCPUs}"
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
    tuple val(key), path(trimmed), path(assembly), val(metagenomic)

    output:
    path("${key}.bam")

    shell:
    '''
    bowtie !{params.bowtieAlignmentOptions} --sam -l 20 !{assembly}/contigs <(zcat !{trimmed}) !{key}.sam
    samtools view -b -F 4 !{key}.sam | samtools sort - > !{key}.bam
    '''
}
