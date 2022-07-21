#!/usr/bin/env nextflow

workflow {
    main:
        prepare_fastaCDA = Channel.from(params.sampleSpecs)
            .map{ [ name: it.metagenomic,
                    metagenomic: "${params.resultsPrefix}/" + it.metagenomic + "-assembly/contigs.fasta",
                    script: "${params.analysisScriptDir}/prepare-fastaCDS.R"] }
            .unique()


        prepare_riboDat = Channel.from(params.sampleSpecs)
            .map{ [ name: it.name,
                    metariboseq: "${params.resultsPrefix}/" + it.metariboseq + ".bam",
                    script: "${params.analysisScriptDir}/prepare-ribodat.R"] }

        samples = Channel.from(params.sampleSpecs)
            .map { [it.name, it.metagenomic, it.metariboseq] }

        fastaCDS = prepareFastaCDA(prepare_fastaCDA)
        riboDat = prepareRiboDat(prepare_riboDat)

        combined = samples.combine(riboDat, by: 0)
            .map{ [it[1], it[0], it[2], it[3]]} // move the metagenome to the front
            .combine(fastaCDS, by: 0)
            .map{ [it[1], it[4], it[3]]} // keep the sample name, fastaCDS and ribodat

        prepare_counts = combined
            .map{ [it[0], it[1], it[2], "${params.analysisScriptDir}/prepare-counts.R"]}

        with_counts = prepareCounts(prepare_counts)

        plots = with_counts
            .map{ [it[0], it[1], it[2], it[3], "${params.analysisScriptDir}/plotVisualsRibo.R"] }

        generatePlots(plots)
}

process prepareFastaCDA {
    memory "${params.assemblyMemory}"

    input:
        tuple val(name), path(metagenomic), path(script)

    output:
        tuple val(name), path("${name}-fastaCDS.RData")

    shell:
    """
    Rscript !{script} !{metagenomic} !{name}-fastaCDS.RData
    """

    stub:
    """
    touch "${name}-fastaCDS.RData"
    """
}

process prepareRiboDat {
    memory "${params.assemblyMemory}"

    input:
        tuple val(name), path(metariboseq), path(script)

    output:
        tuple val(name), path("${name}-ribodata.RData")

    shell:
    """
    Rscript !{script} !{metariboseq} !{name}-ribodata.RData
    """

    stub:
    """
    touch "${name}-ribodata.RData"
    """
}

process prepareCounts {
    memory "${params.assemblyMemory}"

    input:
        tuple val(name), path(fastaCDS), path(ribodat), path(script)

    output:
        tuple val(name), path(fastaCDS), path(ribodat), path("${name}-counts.RData")

    shell:
    """
    Rscript !{script} !{fastaCDS} !{ribodat} !{name}-counts.RData
    """

    stub:
    """
    touch "${name}-counts.RData"
    """
}

// Plots are generated using the metagenomic assembly and the metariboseq
// alignments.
process generatePlots {
    memory "${params.assemblyMemory}"
    publishDir "${params.resultsPrefix}/${name}-plots/", mode: "copy"

    input:
        tuple val(name), path(metagenomic),path(metariboseq), path(counts), path(script)

    output:
        path("*.pdf")

    shell:
    """
    Rscript !{script} !{metagenomic} !{metariboseq} !{counts}
    """

    stub:
    """
    echo Rscript ${script} ${metagenomic} ${metariboseq} ${counts} > dummy-stub.pdf
    """
}
