# Snakefile for running instrain on many samples
# mapping to a single set of genomes
# takes the input of dRep and goes from there
# works by mapping all reads from all samples against
# the combined contig set, extracting bam files for each cluster, 
# and running instrain from that set of bams. 
from os.path import join, abspath, expanduser, exists, basename
from os import makedirs
from pathlib import Path
localrules: make_cluster_bed, index_cluster_bam, idxstats, aggregate_idxstats, drep_fai, bwa_index_contigs, decide_clusters, genome_coverage_count, profile_completed_aggregate

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
drep_folder = config['drep_folder']
drep_contigs = join(drep_folder, 'all_dereplicated_genomes.fa')
drep_fasta_dir = join(drep_folder, "dereplicated_genomes")
# read filters
instrain_min_reads = config['instrain_min_reads']
instrain_max_reads = config['instrain_max_reads']
###################################################################

# if no clusters file specified, will use the determination method from 
# this pipeline which looks at coverage
if cluster_file is not None:
    cluster_to_fasta = get_cluster_fasta(cluster_file)
    cluster_list = list(cluster_to_fasta.keys())
    # limit clusters if desired
    if limit_clusters is not None:
        limit_list = [int(a) for a in str(limit_clusters).split(",")]
        cluster_list = [cluster_list[a] for a in limit_list]
else: 
    cluster_list =[]

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
        expand(join(outdir, "map_reads_drep_sorted/{sample}.bam"), sample=sample_list),
        expand(join(outdir, "map_reads_drep_sorted/{sample}.idxstats.agg"), sample=sample_list),
        join(outdir, "genome_coverage_count.txt"),
        join(outdir, "top_clusters.txt"),
        expand(join(drep_folder, 'cluster_bedfiles/{cluster}.bed'), cluster=cluster_list),
        expand(join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/{sample}.bam.bai"),
            cluster=cluster_list, sample=sample_list),
        expand(join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_scaffold_info.tsv"), 
            cluster=cluster_list, sample=sample_list),
        join(outdir, 'completed.txt')
        # expand(join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_mapping_info.tsv"), 
        #     cluster=cluster_list, sample=sample_list),
        # expand(join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/figures/{cluster}_inStrainCompare_dendrograms.pdf"), cluster=cluster_list),
        # expand(join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/heatmaps/popANI/popANI_heatmap_unfiltered_complete.pdf"), cluster=cluster_list),
        # expand(join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/figures/{cluster}_inStrainCompare_dendrograms.pdf"), cluster=cluster_list),
        # expand(join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/heatmaps/popANI/popANI_heatmap_unfiltered_complete.pdf"), cluster=cluster_list),

# make fai if it doesn't exist
rule drep_fai: 
    input:
        drep_contigs=drep_contigs
    output:
        fai=drep_contigs + '.fai'
    shell: """
        samtools faidx {input}
    """

rule bwa_index_contigs:
    input:
        drep_contigs
    output:
        drep_contigs + '.bwt'
    threads: 1
    resources:
        mem=32,
        time=48
    shell: """
        bwa index {input}
    """

rule map_reads_contigs:
    input:
        drep_contigs = drep_contigs,
        index = rules.bwa_index_contigs.output, 
        reads=lambda wildcards: sample_reads[wildcards.sample],
    output: 
        bam=join(outdir, "map_reads_drep_sorted/{sample}.bam"),
        bai=join(outdir, "map_reads_drep_sorted/{sample}.bam.bai"),
    threads: 8 
    resources:
        mem=32,
        time=lambda wildcards, attempt: 24 * attempt
    params: 
        samtools_threads=7,
        barcode_map= '-C' if barcodes else ''
    shell: """
        bwa mem {params.barcode_map} -t {threads} {input.drep_contigs} {input.reads} | \
        samtools view -@ {params.samtools_threads} -b -F 4 | \
        samtools sort -@ {params.samtools_threads} > {output.bam}
        samtools index -@ {params.samtools_threads} {output.bam}
    """

rule idxstats:
    input:
        bam=join(outdir, "map_reads_drep_sorted/{sample}.bam"),
    output:
        idxstats=join(outdir, "map_reads_drep_sorted/{sample}.idxstats")
    shell: """
        samtools idxstats {input} > {output}
    """

rule aggregate_idxstats:
    input:
        idxstats=join(outdir, "map_reads_drep_sorted/{sample}.idxstats")
    output:
        idxstats=join(outdir, "map_reads_drep_sorted/{sample}.idxstats.agg")
    params:
        rl=int(config['read_length'])
    shell: """
        rl={params.rl}
        echo -e "bin\tlength\tmapped\tunmapped\tcoverage" > {output}
        # need something to check if were using the output of the drep pipeline
        # where contigs of the same genome are indicated with __ designation
        # or a single reference genome, where they will probbaly not have that spec
        if grep -q "__" {input}; then 
            cat {input} | grep -v "\*" | sed "s/.fa__/.fa\t/g" | sed "s/.fasta__/.fasta\t/g" | cut -f 1,3,4,5 | \
                awk 'BEGIN {{FS=OFS="\t"}}  {{ b[$1]; for(i=2;i<=NF;i++)a[$1,i]+=$i }} END {{for( i in b) {{printf("%s",i);for(j=2;j<=NF;j++) {{printf("%s%s",OFS,a[i,j])}} print ""}}}}' | \
                awk -v rl=$rl 'BEGIN {{FS=OFS="\t"}} {{print $1,$2,$3,$4,$3/$2*rl}}' | \
                awk 'BEGIN {{FS=OFS="\t"}} {{$4=sprintf("%.5f",$4)}}7' >> {output}
        else 
            cat {input} | grep -v "\*" | \
                awk 'BEGIN {{FS=OFS="\t"}}  {{ b[$1]; for(i=2;i<=NF;i++)a[$1,i]+=$i }} END {{for( i in b) {{printf("%s",i);for(j=2;j<=NF;j++) {{printf("%s%s",OFS,a[i,j])}} print ""}}}}' | \
                awk -v rl=$rl 'BEGIN {{FS=OFS="\t"}} {{print $1,$2,$3,$4,$3/$2*rl}}' | \
                awk 'BEGIN {{FS=OFS="\t"}} {{$4=sprintf("%.5f",$4)}}7' >> {output}
        fi
    """

# genomes that have greater than 1x coverage in greater than two samples
rule genome_coverage_count: 
    input:
        expand(join(outdir, "map_reads_drep_sorted/{sample}.idxstats.agg"), sample=sample_list),
    output: 
        coverage = join(outdir, "genome_coverage_count.txt"),
        clusters = join(outdir, "top_clusters.txt")
    params:
        idxstats_string = join(outdir, "map_reads_drep_sorted/*.idxstats.agg"),
        min_coverage = 1,
        min_samples = 2,
        drep_folder = join(drep_folder),
        drep_genome_folder = join(drep_folder, "drep_actual/dereplicated_genomes"),
        Wdb = join(drep_folder, "drep_actual/data_tables/Wdb.csv")
    shell: """
        echo -e "genome\tsamples_with_coverage" > {output.coverage}
        tail -n +2 -q {params.idxstats_string} | \
            awk 'BEGIN {{FS=OFS="\t"}} {{if ($5 > {params.min_coverage}) print}}' | cut -f 1 | sort | uniq -c |\
            sort -rn | tr -s " " | awk 'BEGIN {{OFS="\t"}} {{print $2,$1}}' >> {output.coverage}

        # get the drep clusters that correspond to these genomes
        rm -f {output.clusters}
        echo -e "cluster\tgenome_fasta" > {output.clusters}
        while read line; do 
            ctg=$(echo $line | cut -f1 -d " ")
            n=$(echo $line | cut -f2 -d " ")
            if [ $ctg != "genome" ]; then
                if [ $n -ge {params.min_samples} ]; then
                    clus=$(grep -E  "$ctg.fa|$ctg," {params.Wdb} | cut -f2 -d ",")
                    echo -e "$clus\t{params.drep_genome_folder}/$ctg" >> {output.clusters}
                fi
            fi
        done < {output.coverage}
    """

checkpoint decide_clusters:
    input:
        rules.genome_coverage_count.output[1]
    output: 
        cluster_dir = directory(join(outdir, 'drep_alignment_comparison/cluster_files'))
    run:
        cluster_to_fasta = get_cluster_fasta(input[0])
        cluster_list = list(cluster_to_fasta.keys())
        makedirs(output.cluster_dir)
        for c in cluster_list:
            Path(join(output.cluster_dir, c + '.txt')).touch()

# snakemake function to get the clusters and re-evaluate the DAG
def aggregate_clusters(wildcards):
    if (len(cluster_list)) > 0:
        print('DEBUG: returning cluster_list directly')
        # print(cluster_list)
        return(cluster_list)

    checkpoint_output = checkpoints.decide_clusters.get(**wildcards).output[0]
    # print(checkpoint_output)
    # to_return = expand(join(outdir, 'drep_alignment_comparison/cluster_bams/{cluster}'),
    #    cluster=glob_wildcards(join(checkpoint_output, "{cluster}.txt")).cluster)
    clusters = glob_wildcards(join(checkpoint_output, "{cluster}.txt")).cluster
    # clusters = clusters[0:25]
    # print(clusters)
    return clusters

# make a bedfile for each cluster
rule make_cluster_bed:
    input:
        fai = rules.drep_fai.output,
        cluster_file = rules.genome_coverage_count.output[1]
    output: 
        bedfile = join(drep_folder, 'cluster_bedfiles/{cluster}.bed')
    params: 
        bed_folder = join(drep_folder, 'cluster_bedfiles')
    shell: """
        f1=$(grep -P "^{wildcards.cluster}\t" {input.cluster_file} | cut -f2)
        f=$(basename $f1 | sed 's/.fasta//g' | sed 's/.fna//g'| sed 's/.fa//g')
        grep "$f" {input.fai} | awk 'BEGIN {{OFS=FS="\t"}}  {{ print $1,1,$2 }}' > {output.bedfile}
    """

# sneakily use more threads for this because CPU usage was much below 100%
rule extract_cluster_bam: 
    input: 
        bam = join(outdir, "map_reads_drep_sorted/{sample}.bam"),
        bedfile = join(drep_folder, 'cluster_bedfiles/{cluster}.bed'),
        idxstats_agg = join(outdir, "map_reads_drep_sorted/{sample}.idxstats.agg"),
        cluster_file = rules.genome_coverage_count.output[1]
    output:
        bam = join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/{sample}.bam"),
    threads: 1
    resources: 
        time = lambda wildcards, attempt: attempt + 6
    params:
        samtools_threads=3,
        min_reads = instrain_min_reads
    shell: """
        f1=$(grep -P "^{wildcards.cluster}\t" {input.cluster_file} | cut -f2)
        f=$(basename $f1 | sed 's/.fasta//g' | sed 's/.fna//g'| sed 's/.fa//g')
        echo $f
        rc=$(grep $f {input.idxstats_agg} | cut -f 3)
        echo "Readcounts in idxstats: {wildcards.sample} {wildcards.cluster} $rc"
        if [ $rc -gt {params.min_reads} ]; then
            samtools view -@ {params.samtools_threads} {input.bam} -L {input.bedfile} -f 3 -b > {output.bam}
        else
            # fake the output
            touch {output}
        fi 
        # rc=$(samtools flagstat {input} | grep "read1" | cut -f 1 -d " ")
        # if [ $rc -eq 0 ]; then
        #     samtools view {input.bam} -f 3 | head -n 1 | samtools view -b > {output.bam}
        # fi
    """

# get number of reads mapping to this cluster, and subsample if necessary
rule subsample_bam:
    input:
        bam = join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/{sample}.bam")
    output: 
        readcounts = join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/readcounts/{sample}.txt")
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
            rm {input}.bak
        fi
    """

rule index_cluster_bam: 
    input: 
        bam = join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/{sample}.bam"),
        readcounts = join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/readcounts/{sample}.txt")
    output:
        bai = join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/{sample}.bam.bai"),
    threads: 1
    params:
        samtools_threads=3
    shell: """
        # skip if readcounts are zero
        rc="$(cat {input.readcounts})"
        if [ $rc -gt 0 ]; then
            samtools index -@ {params.samtools_threads} {input.bam}
        else
            touch {output}
        fi
    """



# profile for all samples, for each cluster
rule instrain_profile:
    input:
        bam = join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/{sample}.bam"),
        bai = join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/{sample}.bam.bai"),
        readcounts = join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/readcounts/{sample}.txt"),
        cluster_file = rules.genome_coverage_count.output[1]
    output:
        scaffold = join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_scaffold_info.tsv"),
        mapping = join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_mapping_info.tsv")
    threads: 4
    resources: 
        time = lambda wildcards, attempt: 4 ** attempt,
        mem=64
    params:
        outdir = join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}"),
        logfile = join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/log/log.log"),
        drep_fasta_folder = join(drep_folder, "dereplicated_genomes"),
        min_reads = instrain_min_reads
    shell: """
        f=$(grep -P "^{wildcards.cluster}\t" {input.cluster_file} | cut -f2)
        # skip if readcounts are less than the min reads
        rc="$(cat {input.readcounts})"
        if [ $rc -gt {params.min_reads} ]; then
            inStrain profile {input.bam} $f -o {params.outdir} -p {threads} || true
            
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


rule profile_completed_aggregate:
    input: 
        lambda wildcards: expand(join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_scaffold_info.tsv"),
            sample=sample_list, cluster=aggregate_clusters(wildcards))
    output: 
        join(outdir, 'completed.txt')
    shell: """
        touch {output}
    """
