from os.path import join, abspath, expanduser

# get samples and assemblies 
def get_sample_assemblies(sample_file):
    sample_assemblies = {}
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            if len(s) != 2:
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
# sample_reads, sample_assemblies = get_sample_assemblies_reads(sample_file)
# sample_list = list(sample_reads.keys())
sample_assemblies = get_sample_assemblies(sample_file)
sample_list = list(sample_assemblies.keys())

print('##################################################################')
print(' SAMPLE LIST ')
print(sample_list)
print('##################################################################')

# convert outdir to absolute path
if outdir[0] == '~':
    outdir = expanduser(outdir)
outdir = abspath(outdir)

# settings of database and file directories for VIBRANT
dbdir="/labs/asbhatt/bsiranos/databases/vibrant-1.2.0/databases"
fdir="/labs/asbhatt/bsiranos/databases/vibrant-1.2.0/files"

rule all:
    input:
        expand(join(outdir, "{sample}/VIBRANT_{sample}/VIBRANT_results_{sample}/VIBRANT_summary_results_{sample}.tsv"), sample=sample_list)


rule run_vibrant:
    input:
        lambda wildcards: sample_assemblies[wildcards.sample]
    output:
        join(outdir, "{sample}/VIBRANT_{sample}/VIBRANT_results_{sample}/VIBRANT_summary_results_{sample}.tsv")
    params:
        dbdir=dbdir, 
        fdir=fdir, 
        outdir = join(outdir, "{sample}")
    threads: 8
    resources: 
        mem=64,
        time=lambda wildcards, attempt: attempt * 4
    shell: """
        VIBRANT_run.py -i {input} -t {threads} -folder {params.outdir} -d {params.dbdir} -m {params.fdir}
    """
