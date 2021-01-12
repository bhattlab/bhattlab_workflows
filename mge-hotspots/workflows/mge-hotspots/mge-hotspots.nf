#!/usr/bin/env nextflow

// preprocess performs simple data manipulation and preparation for the
// hotspot identification computation carried out by analyze_window_coverage.
process preprocess {
    input:
        path species from "${params.speciesInput}"
        path genotypes from "${params.genotypesInput}"
        path clusters from "${params.clustersInput}"
        path reference from "${params.referenceDir}"
        path make_insertions_py from "${projectDir}/02.make_mobile_insertion_bed_files.py"
        path genome_lengths_py from "${projectDir}/03.get_genome_lengths.py"
        path sliding_windows_py from "${projectDir}/04.make_sliding_windows.py" 
        path window_coverage_py from "${projectDir}/05.get_window_coverage.py"
    output:
        path("outputs/insert_beds/*") into insert_beds
        path("outputs/window_counts/*") into window_counts
        path("outputs/sliding_windows/*") into sliding_windows
        path("outputs/genome_lengths/*") into genome_lengths

    shell:
    '''
    python !{make_insertions_py} --species=!{species} \
        --genotypes=!{genotypes} --clusters=!{clusters} --output-prefix=./outputs
    python !{genome_lengths_py} --species=!{species} \
        --genomes=!{reference} --output-prefix=./outputs
    python !{sliding_windows_py} --species=!{species} --output-prefix=./outputs
    python !{window_coverage_py} --species=!{species} --output-prefix=./outputs
    '''
}

insert_beds.into{ analyze_insertions; analyze_genomewide_insertions}

// prokka_annotations generates the prokka annotations for the species
// configured in the species file. It takes about 20 minutes to run on
// a MacBookPro16,1 i19.
process prokka_annotations {
    input:
        path species from "${params.speciesInput}"
        path reference from "${params.referenceDir}"
        path prokka_py from "${projectDir}/01.create_prokka_annotations.py"
    output:
        path("outputs/prokka_genomes") into prokka_genomes

    shell:
    '''
    python "!{prokka_py}" --species=!{species} --genomes=!{reference} \
        --output-prefix=./outputs
    '''
}

// analyze_window_coverage identifies insertion hotspots based on
// the window coverage data.
process analyze_window_coverage {
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
        path species from "${params.speciesInput}"
        path analyze_coverage_py from "${projectDir}/06.analyze_window_coverage.py"
        path("outputs/insert_beds/*") from analyze_insertions
        path("outputs/window_counts/*") from window_counts
    output:
        path("outputs/window_significance_results/*") into analyzed

    shell:
    '''
    python "!{analyze_coverage_py}" --species=!{species} --output-prefix=./outputs
    '''
}

analyzed.into{ to_summarize; to_plot_hotspots}

process summarize_window_coverage {
    input:
        path species from "${params.speciesInput}"
        path summarize_coverage_py from "${projectDir}/08.summarize_window_coverage.py"
        path("outputs/window_significance_results/*") from to_summarize

    output:
        path("outputs/all_top_windows_detailed.tsv") into plot_genomewide_hotspots
        path("outputs/all_top_windows_summarized.bed") into merge_ranges

    shell:
    '''
    python "!{summarize_coverage_py}" --species=!{species} --output-prefix=./outputs
    '''
}

process merge_ranges {
    input:
        path("outputs/all_top_windows_summarized.bed") from merge_ranges

    output:
        path("outputs/all_top_windows_summarized.merged.bed") into merged_ranges

    shell:
    '''
    bedtools merge -i outputs/all_top_windows_summarized.bed  -c 4 -o min > outputs/all_top_windows_summarized.merged.bed
    '''
}

process get_prokka_coding_sequences {
    input:
        path species from "${params.speciesInput}"
        path get_coding_sequences_py from "${projectDir}/11.get_prokka_coding_sequences.py"
        path("outputs/prokka_genomes") from prokka_genomes

    output:
        path("outputs/prokka_annot/*") into prokka_annotations
        path("outputs/all_prokka_genomes_name_conversion.tsv") into prokka_name_conversions

    shell:
    '''
    python "!{get_coding_sequences_py}" --species=!{species} --output-prefix=./outputs
    '''
}

process genomewide_insertion_hotspots {
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
        path species from "${params.speciesInput}"
        path find_closest_genes_py from "${projectDir}/12.find_closest_genes.py"
        path("outputs/all_top_windows_summarized.merged.bed") from merged_ranges
        path("outputs/insert_beds/*") from analyze_genomewide_insertions
        path("outputs/genome_lengths/*") from genome_lengths
        path("outputs/prokka_annot/*") from prokka_annotations
        path("outputs/sliding_windows/*") from sliding_windows

    output:
        path("outputs/closest_genes/*")

    shell:
    '''
    python !{find_closest_genes_py} --species=!{species} --output-prefix=./outputs
    '''
}

process plot_hotspots {
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
        path species from "${params.speciesInput}"
        path coverage_plots_py from "${projectDir}/07.plot_window_coverage.py"
        path("outputs/window_significance_results/*") from to_plot_hotspots

    output:
        path("outputs/window_significance_plots/*")

    shell:
    '''
    python !{coverage_plots_py} --species=!{species} --output-prefix=./outputs
    '''
}

process plot_genomewide_hotspots {
    echo true
    publishDir "${params.resultsPrefix}", mode: "copy"

    input:
        path species from "${params.speciesInput}"
        path genome_wide_plots_py from "${projectDir}/09.plot_genomewide_insertions.py"
        path("outputs/all_top_windows_detailed.tsv") from plot_genomewide_hotspots

    output:
        path("outputs/genomewide_insertion_hotspots.pdf")

    shell:
    '''
    python !{genome_wide_plots_py} --species=!{species} --output-prefix=./outputs
    '''
}
