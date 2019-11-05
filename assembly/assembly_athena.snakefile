from os.path import join, abspath
import sys
import gzip

localrules: prepare_athena_config

# script for 10x athena assembly
# takes in the output of spades (or other assembly) 
# and preprocessed reads
# does the athena assembly on the results

# helper function for opening gzipped files
def _open(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    return open(filename, mode)

# specify project directories
PROJECT_DIR = config["output_directory"]
# convert PROJECT_DIR to absolute path
if PROJECT_DIR[0] == '~':
    PROJECT_DIR = expanduser(PROJECT_DIR)
PROJECT_DIR = abspath(PROJECT_DIR)

# get input samples from sample table
with open(config["sample_table"]) as inf:
    insamps = [i for i in inf.readlines() if i != '\n']
    sample_dict = {sample: files.split(",") for sample, files in [l.strip().split("\t") for l in insamps]}
# ensure no comment lines
sample_dict = {k:sample_dict[k] for k in sample_dict.keys() if k[0] != '#'}

# get list of samples
sample_list = list(sample_dict.keys())

# directory where each assembly resides 
# relies on details of each assember used 
if config['seed_assembler_used'] == 'spades':
    seed_contig_dict = {s: join(PROJECT_DIR, '02_assembly/02_metaspades/', s, 'contigs.fasta') for s in sample_list}
elif config['seed_assembler_used'] == 'megahit':
    seed_contig_dict = {s: join(PROJECT_DIR, '02_assembly/01_megahit/', s, s + 'contigs.fa') for s in sample_list}
else: 
    sys.exit('Wrong seed_assembler_used specified in config')

rule all:
    input:
        expand(join(PROJECT_DIR, "02_assembly/03_athena/{sample}/athena_config.json"), sample=sample_list),
        expand(join(PROJECT_DIR, "02_assembly/03_athena/{sample}/athena_asm.fa"), sample=sample_list),
        expand(join(PROJECT_DIR, "02_assembly/03_athena/{sample}/quast/report.tsv"), sample=sample_list),
        join(PROJECT_DIR, "02_assembly/03_athena/quast_report_merged.tsv")


################################################
##### ATHENA ASSEMBLY ##########################
################################################

rule interleave_reads:
    input:
        fwd = lambda wildcards: sample_dict[wildcards.sample][0],
        rev = lambda wildcards: sample_dict[wildcards.sample][1]
    output:
        unsorted = temp(join(PROJECT_DIR, "02_assembly/03_athena/{sample}/preprocessed_reads_interleaved.fq"))
    threads: 1 
    resources:
        time = lambda wildcards, attempt: attempt * 12
    run: 
        # original author: SĂŠbastien Boisvert
        # part of Ray distribution

        with _open(input.fwd) as f1:
            with _open(input.rev) as f2:
                with _open(output.unsorted, 'w') as f3:
                    while True:
                        line = f1.readline()
                        if line.strip() == "":
                            break
                        print(line.strip(), file=f3)
                        for i in range(3):
                            print(f1.readline().strip(), file=f3)
                        for i in range(4):
                            print(f2.readline().strip(), file=f3)

rule sort_interleaved:
    input:
        rules.interleave_reads.output
    output:
        sorted = join(PROJECT_DIR, "02_assembly/03_athena/{sample}/preprocessed_reads_interleaved_sorted.fq")
    threads: 1
    resources:
        time = lambda wildcards, attempt: attempt * 12
    shell:
        """
        cat {input} | paste - - - - -d "\r" | sort -k2,2 | tr "\r" "\n" > {output}
        """

rule align_reads_contigs:
    input:
        interleaved_reads = rules.sort_interleaved.output.sorted,
        contigs = lambda wildcards: seed_contig_dict[wildcards.sample]
    output:
        bam = join(PROJECT_DIR, "02_assembly/03_athena/{sample}/align_reads_seed_contigs.bam"),
        bam_bai = join(PROJECT_DIR, "02_assembly/03_athena/{sample}/align_reads_seed_contigs.bam.bai"),
    threads: 8
    resources:
        time = 24,
        mem = 64
    singularity: "shub://bhattlab/bhattlab_workflows:align"
    shell: """
        bwa index {input.contigs}
        bwa mem -t {threads} -C -p {input.contigs} {input.interleaved_reads} | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """

rule prepare_athena_config:
    input:
        interleaved_reads = rules.sort_interleaved.output.sorted,
        contigs = lambda wildcards: seed_contig_dict[wildcards.sample],
        bwa_aligned = rules.align_reads_contigs.output.bam,
        bwa_index = rules.align_reads_contigs.output.bam_bai
    params:
        athena_cores = config['athena_cores']
    output:
        json = join(PROJECT_DIR, "02_assembly/03_athena/{sample}/athena_config.json")
    run:
        options = {
        "ctgfasta_path": input.contigs,
        "reads_ctg_bam_path": input.bwa_aligned,
        "input_fqs": input.interleaved_reads,
        "cluster_settings":{
            "cluster_type": "multiprocessing",
            "processes": config['athena_cores'],
            }
        }
        # print(options)
        with open(output.json, 'w') as outfile:
            json.dump(options, outfile, indent=1)

rule athena_meta:
    input:
        rules.prepare_athena_config.output
    output:
        join(PROJECT_DIR, "02_assembly/03_athena/{sample}/athena_asm.fa")
    params:
        orig_out = join(PROJECT_DIR, "02_assembly/03_athena/{sample}/results/olc/athena.asm.fa")
    singularity: "shub://bsiranosian/bens_1337_workflows:athena"
    threads:  config['athena_cores']
    resources:
        time = 24,
        mem = 400
    shell: """
        athena-meta --config {input}
        cp {params.orig_out} {output}
        """

rule quast_athena:
    input:
        rules.athena_meta.output
    output:
        join(PROJECT_DIR, "02_assembly/03_athena/{sample}/quast/report.tsv")
    params:
        outdir = join(PROJECT_DIR, "02_assembly/03_athena/{sample}/quast/")
    threads: 1
    singularity: "shub://bsiranosian/bens_1337_workflows:athena"
    shell: """
        quast -o {params.outdir} {input} --threads {threads}
        """

rule combine_quast_reports_R:
    input:
        expand(join(PROJECT_DIR, "02_assembly/03_athena/{sample}/quast/report.tsv"), sample=sample_list),
    output:
        join(PROJECT_DIR, "02_assembly/03_athena/quast_report_merged.tsv")
    params:
        sample_names = sample_list,
        assembly_dir = join(PROJECT_DIR, "02_assembly/03_athena/")
    script: "scripts/combine_quast_reports.R"
