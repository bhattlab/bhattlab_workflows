from os.path import join
import sys

################################################################################
# specify project directories
PROJECT_DIR = config["output_directory"]
# get input samples from sample table
with open(config["sample_table"]) as inf:
    insamps = inf.readlines()
    sample_dict = {sample: files.split(",") for sample, files in [l.strip().split("\t") for l in insamps]}
# ensure no comment lines
sample_dict = {k:sample_dict[k] for k in sample_dict.keys() if k[0] != '#'}

# get list of samples
sample_list = list(sample_dict.keys())

# can do more than one assembler here
assemblers = config['assemblers'] 
# ensure at least one valid option
if not ('megahit' in assemblers or 'spades' in assemblers):
    sys.exit('Must have at least one valid assembler in config!')

# define output files depending on which assemblers
outfiles_megahit_assembly = expand(join(PROJECT_DIR, "02_assembly/01_megahit/{sample}/{sample}.contigs.fa"), sample=sample_list)
outfiles_megahit_quast = expand(join(PROJECT_DIR, "02_assembly/01_megahit/{sample}/quast/report.txt"), sample=sample_list)
outfiles_spades_assembly = expand(join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/contigs.fasta"), sample=sample_list)
outfiles_spades_quast = expand(join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/quast/report.txt"), sample=sample_list)

outfiles_all = []
if 'megahit' in assemblers:
    outfiles_all.append(outfiles_megahit_assembly)
    outfiles_all.append(outfiles_megahit_quast)
if 'spades' in assemblers:
    outfiles_all.append(outfiles_spades_assembly)
    outfiles_all.append(outfiles_spades_quast)

################################################################################

rule all:
    input:
        outfiles_all
    
rule megahit:
    input: lambda wildcards: sample_dict[wildcards.sample]
    output: join(PROJECT_DIR, "02_assembly/01_megahit/{sample}/{sample}.contigs.fa")
    resources:
        time = lambda wildcards, attempt: 24 * attempt,
        mem = lambda wildcards, attempt: 100 * attempt
    threads: 8
    params:
        outdir = join(PROJECT_DIR, "02_assembly/01_megahit/{sample}/")
    run:
        reads = sample_dict[wildcards.sample]
        if len(reads) == 3: #paired plus orphans
            cmd = "megahit -1 {reads[0]} -2 {reads[1]} -r {reads[2]}"
        elif len(reads) == 2: #paired
            cmd = "megahit -1 {reads[0]} -2 {reads[1]}"
        elif len(reads) == 1: #single-ended
            cmd = "megahit --12 {reads[0]}"
        cmd += " -o {params.outdir} -t {threads} --out-prefix {wildcards.sample}"
        shell("rm -r {params.outdir}")
        shell(cmd)


rule quast_megahit:
        input:
            rules.megahit.output
        output:
            join(PROJECT_DIR,"02_assembly/01_megahit/{sample}/quast/report.txt")
        params:
            outdir = join(PROJECT_DIR,"02_assembly/01_megahit/{sample}/quast")
        shell: """
            quast.py -o {params.outdir} {input} --fast
            """

rule spades:
    input: lambda wildcards: sample_dict[wildcards.sample]
    output: join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/contigs.fasta")
    resources:
        time = lambda wildcards, attempt: 24 * attempt,
        mem = lambda wildcards, attempt: 100 * attempt
    threads: 8
    params:
        outdir = join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/")
    run:
        reads = sample_dict[wildcards.sample]
        if len(reads) == 3: #paired plus orphans
            cmd = "spades.py --meta -1 {reads[0]} -2 {reads[1]} -s {reads[2]}"
        elif len(reads) == 2: #paired
            cmd = "spades.py --meta -1 {reads[0]} -2 {reads[1]}"
        elif len(reads) == 1: #single-ended
            cmd = "spades.py --meta -12 {reads[0]}"
        cmd += " -o {params.outdir} -m {resources.mem} -t {threads}"
        shell(cmd)

rule quast_spades:
    input:
        rules.spades.output
    output:
        join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/quast/report.txt")
    params:
        outdir = join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/quast")
    resources:
        mem=8,
        time=1
    shell: """
        quast.py -o {params.outdir} {input} --fast
        """
