from os.path import join, abspath, expanduser
import sys

################################################################################
# specify project directories
PROJECT_DIR = config["output_directory"]

GCP = False
if GCP:
    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider(keep_local=True)
    # GS = GSRemoteProvider()
    GS_PREFIX = "gbsc-gcp-lab-bhatt_user-bsiranos"
else:
    # convert PROJECT_DIR to absolute path
    if PROJECT_DIR[0] == '~':
        PROJECT_DIR = expanduser(PROJECT_DIR)
    PROJECT_DIR = abspath(PROJECT_DIR)
    GS_PREFIX=''

def get_sample_group_reads(sample_file):
    sample_reads = {}
    group_reads1 = {}
    group_reads2 = {}
    group_reads_orp = {}
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            sample = s[0]
            if sample in sample_reads:
                print(sample)
                raise ValueError("Non-unique sample encountered!")
            group = s[1]
            reads_split = s[2].split(',')
            # get read pairs and singles from read specification
            if (len(reads_split) == 3):
                reads1 = reads_split[0]
                reads2 = reads_split[1]
                reads_orp = reads_split[2]
                group_reads1.setdefault(group, []).append(reads1)
                group_reads2.setdefault(group, []).append(reads2)
                group_reads_orp.setdefault(group, []).append(reads_orp)
            elif (len(reads_split) == 2):
                reads1 = reads_split[0]
                reads2 = reads_split[1]
                group_reads1.setdefault(group, []).append(reads1)
                group_reads2.setdefault(group, []).append(reads2)
                group_reads_orp.setdefault(group, [])
            elif (len(reads_split) == 1):
                reads_orp = reads_split[0]
                group_reads1.setdefault(group, [])
                group_reads2.setdefault(group, [])
                group_reads_orp.setdefault(group, []).append(reads_orp)
            
            # add these to dictionaries
            sample_reads[sample] = reads_split

    return (sample_reads, group_reads1, group_reads2, group_reads_orp)

# get sample and group info from the table
sample_reads, group_reads1, group_reads2, group_reads_orp = get_sample_group_reads(config['sample_table'])
sample_list = list(sample_reads.keys())
group_list = list(group_reads1.keys())

# print(sample_reads)
# print('     --------------------------------------       ')
# print(group_reads1)
# print('     --------------------------------------       ')
# print(group_reads2)
# print('     --------------------------------------       ')
# print(group_reads_orp)
# print('     --------------------------------------       ')
print(group_list)
print('Total groups: ' + str(len(group_list)))

# can do more than one assembler here
assemblers = config['assemblers'] 
# haven't implemented spades coassembly yet
if ('spades' in assemblers):
    sys.exit('spades coassembly not implemented yet')
# ensure at least one valid option
if not ('megahit' in assemblers or 'spades' in assemblers):
    sys.exit('Must have at least one valid assembler in config!')

# define output files depending on which assemblers
outfiles_megahit_assembly = expand(join(PROJECT_DIR, "01_megahit/{group}/{group}.contigs.fa"), group=group_list)
outfiles_megahit_quast = expand(join(PROJECT_DIR, "01_megahit/{group}/quast/report.tsv"), group=group_list)
outfiles_spades_assembly = expand(join(PROJECT_DIR, "02_metaspades/{group}/contigs.fasta"), group=group_list)
outfiles_spades_quast = expand(join(PROJECT_DIR, "02_metaspades/{group}/quast/report.tsv"), group=group_list)

outfiles_all = []
if 'megahit' in assemblers:
    outfiles_all.append(outfiles_megahit_assembly)
    outfiles_all.append(outfiles_megahit_quast)
    outfiles_all.append(join(GS_PREFIX, PROJECT_DIR, "01_megahit/quast_report_merged.tsv"))
if 'spades' in assemblers:
    outfiles_all.append(outfiles_spades_assembly)
    outfiles_all.append(outfiles_spades_quast)
    outfiles_all.append(join(GS_PREFIX, PROJECT_DIR, "02_metaspades/quast_report_merged.tsv"))

# helper functions to build command for assembly so we can use singularity
def get_spades_reads_command(reads):
    if len(reads) == 3: #paired plus orphans
        return("-1 {0} -2 {1} -s {2}".format(reads[0], reads[1], reads[2]))
    elif len(reads) == 2: #paired
        return("-1 {0} -2 {1}".format(reads[0], reads[1]))
    elif len(reads) == 1: #single-ended
        return("-12 {0}".format(reads[0]))

def get_megahit_reads_command(group):
    # get paired end reads
    reads1 = [join(GS_PREFIX, a) for a in group_reads1[group]]
    reads2 = [join(GS_PREFIX, a) for a in group_reads2[group]]
    reads_orp = [join(GS_PREFIX, a) for a in group_reads_orp[group]]
    if (len(reads1) >0 and len(reads2) >0 and len(reads_orp) >0):
        cmd = "-1 {0} -2 {1} -r {2}".format(','.join(reads1), ','.join(reads2), ','.join(reads_orp))
    elif (len(reads1) >0 and len(reads2) >0):
        cmd = "-1 {0} -2 {1}".format(','.join(reads1), ','.join(reads2))
    elif len(reads_orp) > 0: #single-ended
        cmd = "--12 {0}".format(','.join(reads_orp))
    else: 
        sys.exit('Bad read specification for group: ' + group)
    return(cmd)

################################################################################

rule all:
    input:
        outfiles_all

################################################################################
######## MEGAHIT ASSEMBLY ######################################################
################################################################################
rule megahit:   
    input:
        r1 = lambda wildcards: group_reads1[wildcards.group],
        r2 = lambda wildcards: group_reads2[wildcards.group],
        r_orp = lambda wildcards: group_reads_orp[wildcards.group],
    output: join(PROJECT_DIR, "01_megahit/{group}/{group}.contigs.fa")
    resources:
        time = lambda wildcards, attempt: 24 * attempt,
        mem = lambda wildcards, attempt: 100 * attempt, 
        mem_mb = lambda wildcards, attempt: 100000 * attempt
    threads: 16
    params:
        outdir = join(GS_PREFIX, PROJECT_DIR, "01_megahit/{group}/"),
        reads_command = lambda wildcards: get_megahit_reads_command(wildcards.group),
        mem_max = "100e9"
    singularity: "docker://quay.io/biocontainers/megahit:1.2.9--h8b12597_0"
    shell: """
        rm -rf {params.outdir}
        megahit {params.reads_command} -o {params.outdir} -t {threads} --out-prefix {wildcards.group} -m {params.mem_max}
    """

rule quast_megahit:
    input:
        rules.megahit.output
    output:
        join(PROJECT_DIR,"01_megahit/{group}/quast/report.tsv")
    singularity: "docker://quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0"
    params:
        outdir = join(GS_PREFIX, PROJECT_DIR,"01_megahit/{group}/quast")
    shell: """
        quast.py -o {params.outdir} {input} --fast
        """

rule combine_megahit_quast_reports_R:
    input:
        expand(join(PROJECT_DIR, "01_megahit/{group}/quast/report.tsv"), group=group_list),
    output:
        join(PROJECT_DIR, "01_megahit/quast_report_merged.tsv")
    params:
        sample_names = group_list,
        assembly_dir = join(GS_PREFIX, PROJECT_DIR, "01_megahit/")
    script: "scripts/combine_quast_reports.R"

################################################################################
######## SPADES ASSEMBLY #######################################################
################################################################################
rule spades:
    input: lambda wildcards: sample_dict[wildcards.sample]
    output: join(PROJECT_DIR, "02_metaspades/{group}/contigs.fasta")
    resources:
        time = lambda wildcards, attempt: 24 * attempt,
        mem = lambda wildcards, attempt: 100 * attempt
    threads: 8
    singularity: "docker://quay.io/biocontainers/spades:3.15.3--h95f258a_1"
    params:
        outdir = join(PROJECT_DIR, "02_metaspades/{group}/"),
        reads_command = lambda wildcards: get_spades_reads_command(sample_dict[wildcards.sample])
    shell: """
        spades.py --meta {params.reads_command} -o {params.outdir} -m {resources.mem} -t {threads} --phred-offset 33 --only-assembler
    """
# cmd += " -o {params.outdir} -m {resources.mem} -t {threads} --only-assembler"
# add in --only-assembler if spades gets stuck on read error correction

rule quast_spades:
    input:
        rules.spades.output
    output:
        join(PROJECT_DIR, "02_metaspades/{group}/quast/report.tsv")
    singularity: "docker://quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0"
    params:
        outdir = join(PROJECT_DIR, "02_metaspades/{group}/quast"),
        min_contig = 500
    resources:
        mem=8,
        time=1
    threads: 1
    shell: """
        rm -r {params.outdir}
        mkdir {params.outdir}
        quast.py -o {params.outdir} {input} --fast --min-contig {params.min_contig} -t {threads}
        """

rule combine_spades_quast_reports_R:
    input:
        expand(join(PROJECT_DIR, "02_metaspades/{group}/quast/report.tsv"), group=group_list),
    output:
        join(PROJECT_DIR, "02_metaspades/quast_report_merged.tsv")
    params:
        sample_names = group_list,
        assembly_dir = join(PROJECT_DIR, "02_metaspades/")
    script: "scripts/combine_quast_reports.R"

