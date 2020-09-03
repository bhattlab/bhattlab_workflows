# snakefile for running RGI on a set of contigs 
from os.path import join, abspath, expanduser

################################################################################
# specify project directories
PROJECT_DIR = config["outdir_base"]
# convert PROJECT_DIR to absolute path
if PROJECT_DIR[0] == '~':
    PROJECT_DIR = expanduser(PROJECT_DIR)
PROJECT_DIR = abspath(PROJECT_DIR)

# get input samples from sample table
with open(config["sample_table"]) as inf:
    insamps = [i for i in inf.readlines() if i != '\n']
    sample_dict = {sample: [read1, read2] for sample, read1, read2 in [l.strip().split("\t") for l in insamps]}
# ensure no comment lines
sample_dict = {k:sample_dict[k] for k in sample_dict.keys() if k[0] != '#'}
# get list of samples
SAMPLES = list(sample_dict.keys())
# set the aligner to use
aligner=config["aligner"]
if (aligner not in ["bwa", "bowtie2"]):
    quit('aligner must be one of bwa or bowtie2')

print('##################################################################')
print(' SAMPLE LIST ')
print(SAMPLES)
print('##################################################################')


rule all:
    input:
        expand(join(PROJECT_DIR, "{sample}", "{sample}.allele_mapping_data.txt"), sample = SAMPLES)

rule rgi_read:
    input:
        fwd = lambda wildcards: sample_dict[wildcards.sample][0],
        rev = lambda wildcards: sample_dict[wildcards.sample][1],
    output: 
        join(PROJECT_DIR, "{sample}", "{sample}.allele_mapping_data.txt")
    threads: 8
    resources:
        mem=16,
        time = lambda wildcards, attempt: 2 * attempt
    benchmark: join(PROJECT_DIR, "{sample}", "{sample}_time.txt")
    params: 
        outdir = join(PROJECT_DIR, "{sample}"),
        out_base = join(PROJECT_DIR, "{sample}", "{sample}"),
        aligner = "bwa"
    shell: """
        echo "Start RGI Read: {wildcards.sample}"
        rgi bwt --read_one {input.fwd} --read_two {input.rev} --aligner {params.aligner} --output_file {params.out_base} \
        --threads {threads} --include_wildcard --clean
        echo "Completed RGI"
    """

# commands to load database from Bens folder
# db_folder="/oak/stanford/scg/lab_asbhatt/bsiranos/databases/CARD"
# db_version="3.1.0"
# rgi load --card_json "$db_folder"/"$db_version"/card.json 
# rgi load -i "$db_folder"/"$db_version"/card.json --card_annotation "$db_folder"/"$db_version"/card_database_v"$db_version".fasta 
# rgi load --wildcard_annotation "$db_folder"/"$db_version"/wildcard/wildcard_database_v"$db_version".fasta --wildcard_index "$db_folder"/"$db_version"/wildcard/index-for-model-sequences.txt --card_annotation "$db_folder"/"$db_version"/card_database_v"$db_version".fasta 
# rgi load --kmer_database "$db_folder"/"$db_version"/wildcard/61_kmer_db.json --amr_kmers "$db_folder"/"$db_version"/wildcard/all_amr_61mers.txt --kmer_size 61 --debug > kmer_load.61.log 2>&1
# rgi database -v