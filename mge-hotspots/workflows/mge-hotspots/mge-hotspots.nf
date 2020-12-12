#!/usr/bin/env nextflow

process preprocess {
    output:
        path("outputs") into preprocessed

    shell:
    '''
    python "!{projectDir}/02.make_mobile_insertion_bed_files.py" \
        --species=!{params.speciesInput} \
        --genotypes=!{params.genotypesInput} \
        --clusters=!{params.clustersInput} \
        --output-prefix=./outputs
    python "!{projectDir}/03.get_genome_lengths.py" \
        --species=!{params.speciesInput} \
        --genomes=!{params.referenceDir} \
        --output-prefix=./outputs
    python "!{projectDir}/04.make_sliding_windows.py" \
        --species=!{params.speciesInput} \
        --output-prefix=./outputs
    python "!{projectDir}/05.get_window_coverage.py" \
        --species=!{params.speciesInput} \
        --output-prefix=./outputs
    '''
}

process analyze_window_coverage {
    input:
        path("outputs") from preprocessed
    output:
        path("outputs") into analyzed

    shell:
    '''
    python "!{projectDir}/06.analyze_window_coverage.py" \
        --species=!{params.speciesInput} \
        --output-prefix=./outputs
    '''
}

process plot_window_coverage {
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
        path("outputs") from analyzed

    output:
        path("outputs/window_significance_plots/*")
        path("outputs/window_significance_results/*")


    shell:
    '''
    python "!{projectDir}/07.plot_window_coverage.py" \
        --species=!{params.speciesInput} \
        --output-prefix=./outputs
    '''
}
