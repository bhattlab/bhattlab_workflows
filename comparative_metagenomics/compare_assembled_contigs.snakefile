# compare assembled contigs pipeline
# largely follows compare aligned assembled reads
# start with sample table: name to contigs
# bwa mem to reference geonome
# then follow with standard comparisons

import re
import sys
from os.path import join, abspath
from os import listdir
from itertools import combinations
import string
import random

def random_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

localrules: collect_reports_contigs, nucmer_heatmaps_contigs, filter_contigs

# need at least one reference genome. But can also use a multifasta for alignment
reference_genome_single = config['reference_genome_single']
if 'reference_genome_multi' in config:
    if config['reference_genome_multi'] is not None:
        reference_genome_multi = config['reference_genome_multi']
    else:
        reference_genome_multi = reference_genome_single
else:
    reference_genome_multi = reference_genome_single
    
outdir = config['outdir']
include_reference = config['include_reference_heatmap']
reference_name = config['reference_name']

# convert outdir to absolute path
if outdir[0] == '~':
    outdir = expanduser(outdir)
outdir = abspath(outdir)

def get_sample_contigs(sample_file):
    sample_contigs = {}
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            sample = s[0]
            contigs=s[1]
            if sample in sample_contigs:
                raise ValueError("Non-unique sample encountered!")
            sample_contigs[sample] = contigs
    return sample_contigs


# read in sample info and reads from the sample_contig_list
sample_contigs = get_sample_contigs(config['sample_contig_list'])
sample_names =  list(sample_contigs.keys())
# sort so combinations are consistent
sample_names.sort()

print('Sample names: ')
print(sample_names)
# print(sample_contigs)
# print(reference_genome_single)
# print(reference_genome_multi)

# make pairwise combinations
pairwise_combination_list = list(combinations(sample_names, 2))
pairwise_combination_names = [a + '__' + b for (a,b) in pairwise_combination_list]
pairwise_combination_contigs = {a + '__' + b:(
    join(outdir, "01_align/", a, "contigs_filtered_aligned.fasta"),
    join(outdir, "01_align/", b, "contigs_filtered_aligned.fasta")) for (a,b) in pairwise_combination_list}
print('Made all pairwise comparisons expected')

# if including the reference genome, add that to pairwise combinations
if include_reference:
    pairwise_combination_reference = [(a, reference_name) for a in sample_names]
    pairwise_combination_reference_names = [a + '__' + b for (a,b) in pairwise_combination_reference]
    pairwise_combination_reference_contigs = {a + '__' + b:(
        join(outdir, "01_align/", a, "contigs_filtered_aligned.fasta"),
        reference_genome_single) for (a,b) in pairwise_combination_reference}
    # append to exisitng
    pairwise_combination_names = pairwise_combination_names + pairwise_combination_reference_names
    pairwise_combination_contigs.update(pairwise_combination_reference_contigs)
    # print(pairwise_combination_names)
    # print(pairwise_combination_contigs)

rule all:
    input:
        expand(join(outdir, '01_align/{name}/contigs_filtered_aligned.fasta'), name=sample_names),
        join(outdir, "01_align/quast_report_merged.tsv"),
        # expand(join(outdir, "02_nucmer_pairwise/contig_reports/{name}.report"), name=pairwise_combination_names),
        # expand(join(outdir, "02_nucmer_pairwise/contig_dotplots/{name}_dotplot.png"), name=pairwise_combination_names),
        # expand(join(outdir, "02_nucmer_pairwise/contig_dotplots_snps/{name}_dotplot_snps.png"), name=pairwise_combination_names),
        # join(outdir, "02_nucmer_pairwise/contig_report_stats.txt"),
        join(outdir, "02_nucmer_pairwise/contig_plots/nucmer_identity_heatmap_complete.pdf"),


rule filter_contigs:
    input:
        lambda wildcards: sample_contigs[wildcards.name],
    output:
        join(outdir, "01_align/{name}/contigs_filtered.fasta")
    params:
        min_length = 500
    shell: """
        if [[ $(wc -l <{input}) -ge 1 ]]; then
            samtools faidx {input}
            cat {input}.fai | awk '{{if ($2 > {params.min_length}) print $1}}' | xargs samtools faidx {input} > {output}
        else
            touch {output}
        fi
    """

rule align_to_reference:
    input:
        contigs = join(outdir, "01_align/{name}/contigs_filtered.fasta"),
        reference = reference_genome_multi,
    output:
        bam = join(outdir, "01_align/{name}/contigs_filtered_aligned.bam"),
    threads: 4
    params:
        mem = 64,
        time = 6,
    shell: """
        # if an index needs to be built, use bwa index ref.fa
        if [[ ! -f {input.reference}.bwt ]]; then
            bwa index {input.reference}
        fi
        # align contigs to multi reference
        # Maybe should be using different params, but start with default for now
        bwa mem -t {threads} {input.reference} {input.contigs} | samtools view -hu -F 4 | samtools sort -@{threads} -o {output.bam} -n -
    """

rule extract_fasta:
    input:
        bam = join(outdir, "01_align/{name}/contigs_filtered_aligned.bam"),
    output:
        fasta = join(outdir, "01_align/{name}/contigs_filtered_aligned.fasta"),
    shell: """
        samtools fasta {input.bam} > {output.fasta}
    """

rule quast_contigs:
    input: 
        contigs = join(outdir, "01_align/{name}/contigs_filtered_aligned.fasta"),
        ref = reference_genome_single,
    output:
        join(outdir, "01_align/{name}/quast/report.tsv"),
    params:
        outdir = join(outdir, "01_align/{name}/quast/")
    threads: 1
    shell: """
        if [[ $(wc -l <{input.contigs}) -ge 1 ]]; then
            quast {input.contigs} -o {params.outdir} -t {threads} -r {input.ref} -s 
        else
            touch {output}
        fi 
    """

rule combine_quast_reports_R:
    input:
        expand(join(outdir, "01_align/{name}/quast/report.tsv"), name=sample_names),
    output:
        join(outdir, "01_align/quast_report_merged.tsv")
    params:
        sample_names = sample_names,
        assembly_dir = join(outdir, "01_align/")
    script: "scripts/combine_quast_reports.R"



########################################################################
##### NUCMER PIPELINE (CONTIGS) ########################################
########################################################################
rule nucmer_contigs:
    input:
        file1 = lambda wildcards: pairwise_combination_contigs[wildcards.name][0],
        file2 = lambda wildcards: pairwise_combination_contigs[wildcards.name][1],
    output:
        report = join(outdir, "02_nucmer_pairwise/contig_reports/{name}.report"),
        dotplot = join(outdir, "02_nucmer_pairwise/contig_dotplots/{name}_dotplot.png"),
        dotplot_snps = join(outdir, "02_nucmer_pairwise/contig_dotplots_snps/{name}_dotplot_snps.png"),
    params:
        outdir_reports = join(outdir, '02_nucmer_pairwise/contig_reports'),
        outdir_dotplots = join(outdir, '02_nucmer_pairwise/contig_dotplots'),
        outdir_dotplots_snps = join(outdir, '02_nucmer_pairwise/contig_dotplots_snps'),
        outdir_original = outdir,
        tmpdir = lambda wildcards: join(outdir, "tmp_" + wildcards.name)
        # tmpdir = lambda wildcards: join(outdir, "tmp_" + random_generator(size=6))
    resources: 
        time = 1,
        mem = 8
    shell: """
        mkdir -p {params.tmpdir}
        cd {params.tmpdir}

        # ensure the resulting contigs actaully have lines
        if [[ $(wc -l <{input.file1}) -ge 1 ]] &&  [[ $(wc -l <{input.file2}) -ge 1 ]]; then
            nucmer --prefix={wildcards.name} {input.file1} {input.file2}
            show-coords -rcl {wildcards.name}.delta > {wildcards.name}.coords

            # show-aligns {wildcards.name}.delta pre post > {wildcards.name}.aligns
            
            # ensure we have data in the alignments before moving along with the pipeline
            if [[ $(wc -l <{wildcards.name}.delta) -ge 3 ]]; then
                mummerplot -f -l -png --large -p {wildcards.name}_dotplot {wildcards.name}.delta
                mummerplot -f -l -png --large -S -p {wildcards.name}_dotplot_snps {wildcards.name}.delta
                # snp detection
                show-snps -Clr {wildcards.name}.delta > {wildcards.name}.snps
                # reports
                dnadiff -d {wildcards.name}.delta -p {wildcards.name}
            else
                # make some fakes
                touch {wildcards.name}_dotplot.png
                touch {wildcards.name}_dotplot_snps.png
                echo -e "AvgIdentity 0\nTotalLength 0\n" > {wildcards.name}.report
            fi
        else
            # make some fakes
            touch {wildcards.name}_dotplot.png
            touch {wildcards.name}_dotplot_snps.png
            echo -e "AvgIdentity 0\nTotalLength 0\n" > {wildcards.name}.report
        fi

        # cleanup
        mkdir -p {params.outdir_reports}
        mkdir -p {params.outdir_dotplots}
        mkdir -p {params.outdir_dotplots_snps}
        mv {wildcards.name}.report {params.outdir_reports}
        mv {wildcards.name}_dotplot.png {params.outdir_dotplots}
        mv {wildcards.name}_dotplot_snps.png {params.outdir_dotplots_snps}
        rm -f {wildcards.name}* 
        cd {params.outdir_original}
        rm -r {params.tmpdir}
    """

rule collect_reports_contigs:
    input:
        expand(join(outdir, "02_nucmer_pairwise/contig_reports/{name}.report"), name=pairwise_combination_names),
    output:
        join(outdir, "02_nucmer_pairwise/contig_report_stats.txt")
    params:
        report_dir = join(outdir, '02_nucmer_pairwise/contig_reports')

    shell: """
        if [ -f {output} ]; then rm {output}; fi 
        for f in {params.report_dir}/*.report; do
            idty=$(grep "AvgIdentity" "$f" | head -n 1 | awk '$1=$1' | cut -f2 -d " ")
            ttl=$(grep "TotalLength" "$f" | head -n 1 | awk '$1=$1' | cut -f2 -d " ")
            comp_name1=$(basename $f)
            comp_name=$(echo $comp_name1 | sed "s/.report//g")
            printf '%s\t%s\t%s\n' "$comp_name" "$idty" "$ttl" >> {output}
        done
    """

rule nucmer_heatmaps_contigs:
    input:
        join(outdir, "02_nucmer_pairwise/contig_report_stats.txt")
    output:
        join(outdir, "02_nucmer_pairwise/contig_plots/nucmer_identity_heatmap_complete.pdf")
    params:
        outdir = join(outdir, "02_nucmer_pairwise/contig_plots"),
        sample_names = sample_names,
        min_identity_plot = 90
    script:
        "scripts/heatmaps_nucmer.R"
