# snakefile for running RGI on a set of contigs 
from os.path import join, abspath, expanduser

def get_sample_assemblies(sample_file):
    sample_assemblies = {}
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            if len(s) < 2:
                sys.exit('Badly formatted sample_file')
            sample = s[0]
            assembly = s[1]
            if sample in sample_assemblies:
                print(sample)
                raise ValueError("Non-unique sample encountered!")
            sample_assemblies[sample] = assembly

    return sample_assemblies

# Read in sample and outdir from config file
sample_file = config['sample_file']
outdir = config['outdir_base']
sample_assemblies = get_sample_assemblies(sample_file)
sample_list = list(sample_assemblies.keys())
# convert outdir to absolute path
if outdir[0] == '~':
    outdir = expanduser(outdir)
outdir = abspath(outdir)

print('##################################################################')
print(' SAMPLE LIST ')
print(sample_list)
print('##################################################################')


rule all:
    input: expand(join(outdir, '{sample}.rgi.txt'), sample=sample_list)

rule rgi:
    input:
        asm = lambda wildcards: sample_assemblies[wildcards.sample]
    output:
        join(outdir, '{sample}.rgi.txt')
    params:
        outstring = join(outdir, '{sample}.rgi')
    threads: 8
    resources:
        mem = 16, 
        time = lambda wildcards, attempt: 2 * attempt 
    shell: """
        echo "Starting RGI: {wildcards.sample}" 
        rgi main --input_sequence {input} --output_file {params.outstring} --input_type contig --clean -n {threads}
        echo "Completed RGI" 
    """
