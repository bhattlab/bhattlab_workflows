# script to download a dataset from SRA, using the
# parallel-fastq-dump package
import sys
from os.path import join

if config['srr']=='' or config['srr'][0:3] != 'SRR':
    sys.exit('Must specify SRR id on the commad line like: --config srr=SRR000000')
srrid = config['srr']

# unsure of what the output will be in pairs, etc
# so make a completed flag
rule all:
    input:
        join(config['srr'], 'completed.txt')

rule dump:
    output:
        join(config['srr'], 'completed.txt')
    resources:
        time = 6,
        mem = 64
    threads: 4
    singularity: "shub://bsiranosian/bens_1337_workflows:fastq-dump"
    shell: """
    rm -f {output}
    parallel-fastq-dump --threads {threads} --outdir {config[srr]} --sra-id {config[srr]} \
    --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip --readids
    touch {output}
    """

