# inStrain profile pipeline, for comparing against a single
# reference genome, rather than a complex set of geneomes 
# from dRep. Won't do any of the complex cluster finding, etc
# can just use the cluster to fasta mapping to define the name
# for the genomes
from os.path import join, abspath, expanduser, exists, basename
from os import makedirs
from pathlib import Path
import sys
localrules: copy_reference, faidx, idxstats, aggregate_idxstats, faidx, bwa_index_contigs

# get mapping of clusters to fasta files 
# that define the genome
def get_cluster_fasta(cluster_file):
    cluster_dict = {}
    with open(cluster_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'cluster' or s[0].startswith('#'):
                continue
            cluster = s[0]
            fasta = s[1]
            if cluster in cluster_dict:
                print(cluster)
                raise ValueError("Non-unique cluster encountered!")
            cluster_dict[cluster] = fasta
    return cluster_dict

# get mapping from sample to read files
def get_sample_reads(sample_file):
    sample_reads = {}
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            if len(s) != 2:
                sys.exit('Badly formatted sample_file')
            sample = s[0]
            reads_split = s[1].split(',')
            if sample in sample_reads:
                print(sample)
                raise ValueError("Non-unique sample encountered!")
            # get read pairs and singles from read specification
            if (len(reads_split) == 3) or (len(reads_split) == 2):
                sample_reads[sample] = reads_split[0:2]
            elif len(reads_split)==1:
                sample_reads[sample] = reads_split[0]
                # sys.exit('must be paired end reads')
    return sample_reads


###################################################################
### CONFIGURATION #################################################
###################################################################
sample_file = config['sample_reads']
cluster_file = config['cluster_file']
limit_clusters = config['limit_clusters']
outdir = config['outdir']
barcodes = config['barcodes']

# read filters
instrain_min_reads = config['instrain_min_reads']
instrain_max_reads = config['instrain_max_reads']
###################################################################

# get the mapping from the cluter to a single 
# fasta file which will be used for mapping
cluster_to_fasta = get_cluster_fasta(cluster_file)
cluster_list = list(cluster_to_fasta.keys())
ref_fasta = list(cluster_to_fasta.values())[0]
ref_cluster = cluster_list[0]
if len(cluster_list) >1:
    sys.exit("Can only use one reference genome in this pipeline") 
# get sample reads and list
sample_reads = get_sample_reads(sample_file)
sample_list = list(sample_reads.keys())

print()
print("#########################################")
print("# cluster list: " + str(cluster_list))
print("# instrain_max_reads: " + str(instrain_max_reads))
print("# instrain_min_reads: " + str(instrain_min_reads))
print("#########################################")
print()
print(' ONLY CONDUCTING INSTRAIN **PROFILE** STEPS IN THIS PIELINE ')
print()

rule all:
    input:
        expand(join(outdir, "map_reads_ref_sorted/{sample}.bam"), sample=sample_list),
        expand(join(outdir, "map_reads_ref_sorted/{sample}.idxstats.agg"), sample=sample_list),
        expand(join(outdir, "map_reads_ref_sorted/{sample}.bam.bai"),
            cluster=cluster_list, sample=sample_list),
        expand(join(outdir, "instrain_profile/{sample}/output/{sample}_scaffold_info.tsv"), 
            cluster=cluster_list, sample=sample_list),

# copy that fasta over to a central location
rule copy_reference:
    input: 
        ref_fasta
    output: 
        join(outdir, 'ref.fa')
    shell: """
        cp {input} {output}
    """

# make fai if it doesn't exist
rule faidx: 
    input:
      rules.copy_reference.output
    output:
        join(outdir, 'ref.fa.fai')
    shell: """
        samtools faidx {input}
    """

rule bwa_index_contigs:
    input:
        rules.copy_reference.output
    output:
        join(outdir, 'ref.fa.bwt')
    threads: 1
    resources:
        mem=32,
        time=48
    shell: """
        bwa index {input}
    """

rule map_reads_contigs:
    input:
        contigs = rules.copy_reference.output,
        index = rules.bwa_index_contigs.output, 
        reads=lambda wildcards: sample_reads[wildcards.sample],
    output: 
        bam=join(outdir, "map_reads_ref_sorted/{sample}.bam"),
        bai=join(outdir, "map_reads_ref_sorted/{sample}.bam.bai"),
    threads: 8 
    resources:
        mem=32,
        time=lambda wildcards, attempt: 24 * attempt
    params: 
        samtools_threads=7,
        barcode_map= '-C' if barcodes else ''
    shell: """
        bwa mem {params.barcode_map} -t {threads} {input.contigs} {input.reads} | \
        samtools view -@ {params.samtools_threads} -b -F 4 | \
        samtools sort -@ {params.samtools_threads} > {output.bam}
        samtools index -@ {params.samtools_threads} {output.bam}
    """

rule idxstats:
    input:
        bam=join(outdir, "map_reads_ref_sorted/{sample}.bam"),
    output:
        idxstats=join(outdir, "map_reads_ref_sorted/{sample}.idxstats")
    shell: """
        samtools idxstats {input} > {output}
    """

rule aggregate_idxstats:
    input:
        idxstats=join(outdir, "map_reads_ref_sorted/{sample}.idxstats")
    output:
        idxstats=join(outdir, "map_reads_ref_sorted/{sample}.idxstats.agg")
    params:
        rl=int(config['read_length'])
    shell: """
        rl={params.rl}
        echo -e "bin\tlength\tmapped\tunmapped\tcoverage" > {output}
        # this is happening on a single genome so should just add 
        # all the contigs together 
        # But its hard to know how the genome and contigs are differentiated
        # just use last underscore and take everything before it

        cat {input} | grep -v "\*" | sed "s/_[^_\t]*\t/\t\t/g" | cut -f 1,3,4,5 | \
            awk 'BEGIN {{FS=OFS="\t"}}  {{ b[$1]; for(i=2;i<=NF;i++)a[$1,i]+=$i }} END {{for( i in b) {{printf("%s",i);for(j=2;j<=NF;j++) {{printf("%s%s",OFS,a[i,j])}} print ""}}}}' | \
            awk -v rl=$rl 'BEGIN {{FS=OFS="\t"}} {{print $1,$2,$3,$4,$3/$2*rl}}' | \
            awk 'BEGIN {{FS=OFS="\t"}} {{$4=sprintf("%.5f",$4)}}7' >> {output}


    """

# get number of reads mapping to this cluster, and subsample if necessary
rule subsample_bam:
    input:
        bam = join(outdir, "map_reads_ref_sorted/{sample}.bam")
    output: 
        readcounts = join(outdir, "readcounts/{sample}.txt")
    params:
        max_reads=instrain_max_reads,
    threads: 4
    shell: """
        if [ -s {input} ]; then
            samtools flagstat {input} | grep "read1" | cut -f 1 -d " " > {output.readcounts}
        else
            echo "0" > {output.readcounts}
        fi

        rc=$(cat {output.readcounts})
        echo "RC: " $rc
        # if greater than 2 million reads, subsample
        if [ $rc -gt {params.max_reads} ]; then
            frac=$(echo "scale=4; {params.max_reads} / $rc" | bc)
            echo $frac
            mv {input} {input}.bak
            samtools view -@ 3 -s $frac -b {input}.bak > {input}
            samtools index -@ 3 {input}
            rm {input}.bak
        fi
    """


# profile for all samples, for each cluster
rule instrain_profile:
    input:
        bam = join(outdir, "map_reads_ref_sorted/{sample}.bam"),
        bai = join(outdir, "map_reads_ref_sorted/{sample}.bam.bai"),
        readcounts = join(outdir, "readcounts/{sample}.txt"),
        ref = rules.copy_reference.output
    output:
        scaffold = join(outdir, "instrain_profile/{sample}/output/{sample}_scaffold_info.tsv"),
        mapping = join(outdir, "instrain_profile/{sample}/output/{sample}_mapping_info.tsv")
    threads: 4
    resources: 
        time = lambda wildcards, attempt: 4 ** attempt,
        mem=64
    params:
        outdir = join(outdir, "instrain_profile/{sample}"),
        logfile = join(outdir, "instrain_profile/{sample}/log/log.log"),
        drep_fasta_folder = join(outdir, "dereplicated_genomes"),
        min_reads = instrain_min_reads
    shell: """
        # skip if readcounts are less than the min reads
        rc="$(cat {input.readcounts})"
        if [ $rc -gt {params.min_reads} ]; then
            inStrain profile {input.bam} {input.ref} -o {params.outdir} -p {threads} || true
            
            # catch error with zero read pairs remaining, and then fake output if thats the case
            if grep -q "no read pairs remain" {params.logfile}; then
                echo "NO FILTERED READS REMAIN, FAKING OUTPUT"
                touch {output}
            else
                inStrain genome_wide -i {params.outdir} -p {threads}
            fi
        else
            echo "# Not more than {params.min_reads} in the input bamfile, sample skipped" > {output.scaffold}
            echo "# Not more than {params.min_reads} in the input bamfile, sample skipped" > {output.mapping}
            fi
    """
