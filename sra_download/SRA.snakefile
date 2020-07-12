# script to download a dataset from SRA, using the
# parallel-fastq-dump package
import sys
from os.path import join

def get_srr(sample_list):
    srr_list = []
    with open(sample_list) as sf:
        for l in sf.readlines():
            s = l.strip().split()
    	    if s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
       	    if len(s) > 1:
                sys.exit('Incorrect format')
            if s[0]=='' or s[0][0:3] not in ['SRR', 'ERR']:
                sys.exit('Must specify SRR/ERR id on the command line like: --config srr=SRR000000')
            sample = s[0]
       	    srr_list.append(sample)
    print(srr_list)
    return srr_list

list2 = config['srr_list']
srr_list = get_srr(list2)

print(expand("{foo}", foo=srr_list))

# unsure of what the output will be in pairs, etc
# so make a completed flag
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
#    params:
#        folder = "fastq_files/{sample}"  -O {params.folder} -t tmp/scratch
    shell: """
    fasterq-dump {wildcards.sample}
    touch {output}
    """
