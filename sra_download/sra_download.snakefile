# script to download a dataset from SRA, using the
# parallel-fastq-dump package
import sys
from os.path import join

def get_srr(sample_list):
    srr_list = []
    with open(sample_list) as sf:
        for l in sf.readlines():
            s = l.strip().split()
            if s == []:
                continue
            if s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            if len(s) > 1:
                sys.exit('Incorrect format')
            if s[0]=='' or s[0][0:3] not in ['SRR', 'ERR']:
                sys.exit('SRA ids must begin with SRR or ERR')
            sample = s[0]
            srr_list.append(sample)
    return srr_list

# ensure srr_list is specified on the command line
if 'srr_list' not in config:
    sys.exit('Specify list of IDs to download in a newline-delimited file, and run snakemake with the \
        argument: --config srr_list=my_list.txt')

list_file = config['srr_list']
srr_list = get_srr(list_file)

print('srr_list: ' + str(srr_list))

# unsure of what the output will be in pairs, etc
# so make a completed flag once the download is successful
rule all:
    input:
        expand(join("{foo}", 'completed.txt'), foo=srr_list)

rule dump:
    output:
        join("{sample}", 'completed.txt')
    resources:
        time = 6,
        mem = 64
#    singularity: "docker://quay.io/biocontainers/sra-tools:2.10.7--pl526haddd2b5_1"
    threads: 4
    log: 
        join("{sample}", 'log.txt')
    shell: """
    # remove previous directory
    rm -rf {wildcards.sample}
    mkdir -p {wildcards.sample}

    prefetch {wildcards.sample} -O {wildcards.sample} -f yes
    fasterq-dump {wildcards.sample} -O {wildcards.sample} --threads {threads} > {log} 2>&1

    # ensure we have a valid output by checking the logfile
    if $(grep -q "reads written" {log}); then 
        touch {output}
    else
        echo "This is a dumb error in fasterq-dump, trying again......"
        exit 1
    fi
    """
