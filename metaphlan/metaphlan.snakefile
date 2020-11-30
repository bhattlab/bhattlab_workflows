# Snakefile for metaphaln3 classification
# works with paired end read files
from os.path import join
import sys
import snakemake
import time
# output base directory
outdir = config['outdir']
db = config['db']
localrules: 

# function to get the sample reads from the tsv
def get_sample_reads(sample_file):
    sample_reads = {}
    paired_end = ''
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            sample = s[0]
            # paired end specified
            if (len(s)==3):
                reads = [s[1],s[2]]
                if paired_end != '' and not paired_end:
                    sys.exit('All samples must be paired or single ended.')
                paired_end = True
            # single end specified
            elif len(s)==2:
                reads=s[1]
                if paired_end != '' and paired_end:
                    sys.exit('All samples must be paired or single ended.')
                paired_end = False
            if sample in sample_reads:
                raise ValueError("Non-unique sample encountered!")
            sample_reads[sample] = reads
    return (sample_reads, paired_end)


# read in sample info and reads from the sample_file
sample_reads, paired_end = get_sample_reads(config['sample_file'])
sample_names = list(sample_reads.keys())
# print(sample_names)

rule all:
    input:
        expand(join(outdir, "results/{samp}.txt"), samp=sample_names),
        join(outdir, "merged_abundance_table_species.txt"),


rule metaphlan3:
    input:
        r1 = lambda wildcards: sample_reads[wildcards.samp][0],
        r2 = lambda wildcards: sample_reads[wildcards.samp][1],
    output:
        join(outdir, "results", "{samp}.txt")
    params:
        bowtie2out = join(outdir, "bowtie2files", "{samp}.out"),
        db = db
    threads: 8
    resources: 
        mem = 64, 
        time = 24
    shell: """
        # need to remove bowtie2 outfile if it exists already
        rm -f {params.bowtie2out}

        metaphlan {input.r1},{input.r2} --nproc {threads} --input_type fastq \
            --bowtie2db {params.db} --bowtie2out {params.bowtie2out} -o {output}
    """

rule merge_tables:
    input:
        expand(join(outdir, "results/{samp}.txt"), samp=sample_names)
    output:
        t1 = join(outdir, "merged_abundance_table.txt"),
        t2 = join(outdir, "merged_abundance_table_species.txt"),
    params:
        results_dir = join(outdir, "results")

    shell: """
        merge_metaphlan_tables.py  {params.results_dir}/*.txt > {output.t1}
        grep -E "s__|clade" merged_abundance_table.txt | sed 's/^.*s__//g' | cut -f2 --complement | sed -e 's/clade_name/species/g' > {output.t1}
    """

