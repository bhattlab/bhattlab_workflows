from os.path import join, abspath, expanduser
import sys

################################################################################
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

# can do more than one assembler here
assemblers = config['assemblers'] 
# ensure at least one valid option
if not ('megahit' in assemblers or 'spades' in assemblers):
    sys.exit('Must have at least one valid assembler in config!')

# define output files depending on which assemblers
outfiles_megahit_assembly = expand(join(PROJECT_DIR, "02_assembly/01_megahit/{sample}/{sample}.contigs.fa"), sample=sample_list)
outfiles_megahit_quast = expand(join(PROJECT_DIR, "02_assembly/01_megahit/{sample}/quast/report.tsv"), sample=sample_list)
outfiles_spades_assembly = expand(join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/contigs.fasta"), sample=sample_list)
outfiles_spades_quast = expand(join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/quast/report.tsv"), sample=sample_list)

outfiles_all = []
if 'megahit' in assemblers:
    outfiles_all.append(outfiles_megahit_assembly)
    outfiles_all.append(outfiles_megahit_quast)
    outfiles_all.append(join(PROJECT_DIR, "02_assembly/01_megahit/quast_report_merged.tsv"))
if 'spades' in assemblers:
    outfiles_all.append(outfiles_spades_assembly)
    outfiles_all.append(outfiles_spades_quast)
    outfiles_all.append(join(PROJECT_DIR, "02_assembly/02_metaspades/quast_report_merged.tsv"))

# helper functions to build command for assembly so we can use singularity
def get_spades_reads_command(reads):
    if len(reads) == 3: #paired plus orphans
        return("-1 {0} -2 {1} -s {2}".format(reads[0], reads[1], reads[2]))
    elif len(reads) == 2: #paired
        return("-1 {0} -2 {1}".format(reads[0], reads[1]))
    elif len(reads) == 1: #single-ended
        return("-12 {0}".format(reads[0]))

def get_megahit_reads_command(reads):
        if len(reads) == 3: #paired plus orphans
            cmd = "-1 {0} -2 {1} -r {2}".format(reads[0], reads[1], reads[2])
        elif len(reads) == 2: #paired
            cmd = "-1 {0} -2 {1}".format(reads[0], reads[1])
        elif len(reads) == 1: #single-ended
            cmd = "--12 {0}".format(reads[0])
        return(cmd)

################################################################################

rule all:
    input:
        outfiles_all

################################################################################
######## MEGAHIT ASSEMBLY ######################################################
################################################################################
rule megahit:
    input: lambda wildcards: sample_dict[wildcards.sample]
    output: join(PROJECT_DIR, "02_assembly/01_megahit/{sample}/{sample}.contigs.fa")
    resources:
        time = lambda wildcards, attempt: 24 * attempt,
        mem = lambda wildcards, attempt: 100 * attempt
    threads: 8
    params:
        outdir = join(PROJECT_DIR, "02_assembly/01_megahit/{sample}/"),
        reads_command = lambda wildcards: get_megahit_reads_command(sample_dict[wildcards.sample])

    singularity: "docker://quay.io/biocontainers/megahit:1.2.9--h8b12597_0"
    shell: """
        rm -rf {params.outdir}
        megahit {params.reads_command} -o {params.outdir} -t {threads} --out-prefix {wildcards.sample}
    """

rule quast_megahit:
    input:
        rules.megahit.output
    output:
        join(PROJECT_DIR,"02_assembly/01_megahit/{sample}/quast/report.tsv")
    singularity: "docker://quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0"
    params:
        outdir = join(PROJECT_DIR,"02_assembly/01_megahit/{sample}/quast")
    shell: """
        quast.py -o {params.outdir} {input} --fast
        """

rule combine_megahit_quast_reports_R:
    input:
        expand(join(PROJECT_DIR, "02_assembly/01_megahit/{sample}/quast/report.tsv"), sample=sample_list),
    output:
        join(PROJECT_DIR, "02_assembly/01_megahit/quast_report_merged.tsv")
    params:
        sample_names = sample_list,
        assembly_dir = join(PROJECT_DIR, "02_assembly/01_megahit/")
    script: "scripts/combine_quast_reports.R"

################################################################################
######## SPADES ASSEMBLY #######################################################
################################################################################
rule spades:
    input: lambda wildcards: sample_dict[wildcards.sample]
    output: join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/contigs.fasta")
    resources:
        time = lambda wildcards, attempt: 24 * attempt,
        mem = lambda wildcards, attempt: 100 * attempt
    threads: 8
    singularity: "docker://quay.io/biocontainers/spades:3.14.0--h2d02072_0"
    params:
        outdir = join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/"),
        reads_command = lambda wildcards: get_spades_reads_command(sample_dict[wildcards.sample])
    shell: """
        spades.py --meta {params.reads_command} -o {params.outdir} -m {resources.mem} -t {threads}
    """
# cmd += " -o {params.outdir} -m {resources.mem} -t {threads} --only-assembler"
# add in --only-assembler if spades gets stuck on read error correction

rule quast_spades:
    input:
        rules.spades.output
    output:
        join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/quast/report.tsv")
    singularity: "docker://quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0"
    params:
        outdir = join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/quast")
    resources:
        mem=8,
        time=1
    shell: """
        quast.py -o {params.outdir} {input} --fast
        """

rule combine_spades_quast_reports_R:
    input:
        expand(join(PROJECT_DIR, "02_assembly/02_metaspades/{sample}/quast/report.tsv"), sample=sample_list),
    output:
        join(PROJECT_DIR, "02_assembly/02_metaspades/quast_report_merged.tsv")
    params:
        sample_names = sample_list,
        assembly_dir = join(PROJECT_DIR, "02_assembly/02_metaspades/")
    script: "scripts/combine_quast_reports.R"

