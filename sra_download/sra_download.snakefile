# script to download a dataset from SRA, using the
# parallel-fastq-dump package
import sys
from os.path import join

if config['srr']=='' or config['srr'][0:3] not in ['SRR', 'ERR']:
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
    threads: 1
    singularity: "shub://bsiranosian/bens_1337_workflows:fastq-dump"
    shell: """
    rm -f {output}
    rm -rf {config[srr]}
    parallel-fastq-dump --threads {threads} --outdir {config[srr]} --sra-id {config[srr]} \
    --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip
    # prefetch {config[srr]}
    # fastq-dump --outdir {config[srr]} --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip {config[srr]} 
    touch {output}
    """

# removed readids command form this