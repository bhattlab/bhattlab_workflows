#!/usr/bin/env nextflow

toAnalyze = Channel.from(params.sampleSpecs)

// Plots are generated using the metagenomic assembly and the metariboseq
// alignments.
process generatePlots {
    memory "${params.assemblyMemory}"

    input:
        val(sample) from toAnalyze
        path script from "${params.analysisScript}"

    publishDir "${params.resultsPrefix}/${sample.name}-plots/", mode: "copy"

    output:
        path("*.pdf")

    shell:
    """
    Rscript !{script} !{params.resultsPrefix}/!{sample.name}-assembly/contigs.fasta !{params.resultsPrefix}/!{sample.metariboseq}.bam
    """
}
