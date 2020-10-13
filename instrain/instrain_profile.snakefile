# Snakefile for running instrain on many samples
# mapping to a single set of genomes
localrules: make_cluster_bed, index_cluster_bam, idxstats, drep_fai, bwa_index_contigs

from os.path import join, abspath, expanduser, exists, basename

# starts with...
 # reads 
 # drep contig set

# Does ...
 # maps reads against drep contig set
 # extract mapping read for each cluster set

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
sample_file=config['sample_reads']
cluster_file=config['cluster_file']
limit_clusters=config['limit_clusters']
outdir=config['outdir']
barcodes=config['barcodes']
drep_folder=config['drep_folder']
drep_contigs=config['drep_contigs']
drep_fasta_dir = join(drep_folder, "dereplicated_genomes")
###################################################################

cluster_to_fasta = get_cluster_fasta(cluster_file)
cluster_list = list(cluster_to_fasta.keys())
# limit clusters if desired
if limit_clusters is not None:
    limit_list = [int(a) for a in str(limit_clusters).split(",")]
    cluster_list = [cluster_list[a] for a in limit_list]
sample_reads = get_sample_reads(sample_file)
sample_list = list(sample_reads.keys())


# max reads that we'll subsample down to if there's more than this many
instrain_max_reads = 2000000
# minimum number of filtered reads for a sample to be included in the "filtered"
#   instrain compare version
# ALSO minimum number of reads in the bamfile for a sample to be profiled
instrain_filtered_reads = 20000


print()
print("#########################################")
print("# cluster list: " + str(cluster_list))
print("# instrain_max_reads: " + str(instrain_max_reads))
print("# instrain_filtered_reads: " + str(instrain_filtered_reads))
print("#########################################")
print()
print(' ONLY CONDUCTING INSTRAIN **PROFILE** STEPS IN THIS PIELINE ')
print()

rule all:
    input:
        expand(join(outdir, "map_reads_drep_sorted/{sample}.bam"), sample=sample_list),
        # expand(join(outdir, "map_reads_drep_sorted/{sample}.idxstats"), sample=sample_list),
        expand(join(drep_folder, 'cluster_bedfiles/{cluster}.bed'), cluster=cluster_list),
        expand(join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/{sample}.bam.bai",
            cluster=cluster_list, sample=sample_list),
        expand(join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_scaffold_info.tsv"), 
            cluster=cluster_list, sample=sample_list),
        # expand(join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_mapping_info.tsv"), 
        #     cluster=cluster_list, sample=sample_list),
        # expand(join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/figures/{cluster}_inStrainCompare_dendrograms.pdf"), cluster=cluster_list),
        # expand(join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/heatmaps/popANI/popANI_heatmap_unfiltered_complete.pdf"), cluster=cluster_list),
        # expand(join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/figures/{cluster}_inStrainCompare_dendrograms.pdf"), cluster=cluster_list),
        # expand(join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/heatmaps/popANI/popANI_heatmap_unfiltered_complete.pdf"), cluster=cluster_list),

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
        rl=150
    shell: """
        rl={params.rl}
        echo -e "bin\tlength\tmapped\tunmapped\tcoverage" > {output}
        cat {input} | grep -v "\*" | sed "s/.fa__/.fa\t/g" | cut -f 1,3,4,5 | \
        awk 'BEGIN {{FS=OFS="\t"}}  {{ b[$1]; for(i=2;i<=NF;i++)a[$1,i]+=$i }} END {{for( i in b) {{printf("%s",i);for(j=2;j<=NF;j++) {{printf("%s%s",OFS,a[i,j])}} print ""}}}}' | \
        awk -v rl=$rl 'BEGIN {{FS=OFS="\t"}} {{print $1,$2,$3,$4,$3/$2*rl}}' | \
        awk 'BEGIN {{FS=OFS="\t"}} {{$4=sprintf("%.5f",$4)}}7' >> {output}
    """

# making the clsuter file
# to just run outside of snakemake, its easier
'''
cd instrain/map_reads_drep_sorted/
ls *.bam | sed "s/.bam//g" | xargs -P 16 -I {} sh -c "echo {}; samtools idxstats {}.bam > {}.idxstats"
for i in $(ls *.bam | sed "s/.bam//g"); do
    echo $i
    samtools idxstats "$i".bam > "$i".idxstats
done
# awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$3/$2}' | awk 'BEGIN {FS=OFS="\t"} {$4=sprintf("%.5f",$4)}7' > "$i".idxstats

# needs to be aggregated across the whole genome/bin, if the bin is split 
# read length for coverage calculation
rl=150
# for i in *.idxstats; do echo $i; 

i=map_reads_drep/idxstats/bas_bb_spri_10x.idxstats
echo -e "bin\tlength\tmapped\tunmapped\tcoverage" > test.agg
cat "$i" | grep -v "\*" | sed "s/.fa__/.fa\t/g" | cut -f 1,3,4,5 | \
awk 'BEGIN {FS=OFS="\t"}  { b[$1]; for(i=2;i<=NF;i++)a[$1,i]+=$i } END {for( i in b){printf("%s",i);for(j=2;j<=NF;j++){printf("%s%s",OFS,a[i,j])} print ""}}' | \
awk -v rl=$rl 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$3/$2*rl}' | \
awk 'BEGIN {FS=OFS="\t"} {$4=sprintf("%.5f",$4)}7' >> test.agg


# get clusters with X samples with over 1x coverage
cat *.idxstats.agg | awk 'BEGIN {FS=OFS="\t"} {if ($4 > 1) print}' | cut -f 1 | sort | uniq -c | sort -rn | tr -s " " | awk 'BEGIN {OFS="\t"} {print $2,$1}' > ../../contig_coverage_count.txt
cd ../../

# get drep clusters these belong to
drep_folder=drep
rm -f contig_clusters.txt
while read line; do ctg=$(echo $line | cut -f1 -d " "); grep -E  "$ctg.fa|$ctg," "$drep_folder"/data_tables/Wdb.csv | cut -f2 -d "," >> contig_clusters.txt; done  < contig_coverage_count.txt

# put it together in a top_clusters file
genome_folder="$drep_folder/dereplicated_genomes"
paste contig_clusters.txt contig_coverage_count.txt | awk -v f=$genome_folder 'BEGIN {FS=OFS="\t"} {print $1,f"/"$2}' > top_clusters_all.txt
head -n 50 top_clusters_all.txt > top_clusters_50.txt
head -n 10 top_clusters_all.txt  > top_clusters.txt

Could potentially aggregate acoss memebers of the same primary cluster 

# getting the genome information from each of these clusters 
echo "cluster" > tmp1
cat contig_clusters.txt >> tmp1
# from the binning file and contig_coverage_count.txt
head -n 1 binning_das_tool/binning_table_all_simple.tsv > tmp2
while read line; do s=$(echo $line| cut -f 1 -d" " |sed "s/__/\t/g" | sed "s/.fa//g"); grep -w "$s" binning_das_tool/binning_table_all_simple.tsv >> tmp2; done < contig_coverage_count.txt 
paste tmp1 tmp2 > cluster_matching_bins.txt
rm tmp1 tmp2

# for the classification    


'''

# make fai if it doesn't exist
rule drep_fai: 
    input:
        drep_contigs=drep_contigs
    output:
        fai=drep_contigs + '.fai'
    shell: """
        samtools faidx {input}
    """

# make a bedfile for each cluster
# need a mapping from cluster to fasta
rule make_cluster_bed:
    input:
        fai = rules.drep_fai.output
    output: 
        bedfile = join(drep_folder, 'cluster_bedfiles/{cluster}.bed')
    params: 
        fasta_search = lambda wildcards: cluster_to_fasta[wildcards.cluster],
        bed_folder = join(drep_folder, 'cluster_bedfiles')
    shell: """
        f=$(basename {params.fasta_search} | sed 's/.fa//g'| sed 's/.fna//g')
        grep "$f" {input.fai} | awk 'BEGIN {{OFS=FS="\t"}}  {{ print $1,1,$2 }}' > {output.bedfile}
    """

# sneakily use more threads for this because CPU usage was much below 100%
rule extract_cluster_bam: 
    input: 
        bam = join(outdir, "map_reads_drep_sorted/{sample}.bam"),
        bedfile = join(drep_folder, 'cluster_bedfiles/{cluster}.bed'),
        idxstats_agg = join(outdir, "map_reads_drep_sorted/{sample}.idxstats.agg")
    output:
        bam = join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/{sample}.bam"),
    threads: 1
    params:
        samtools_threads=3,
        fasta_search = lambda wildcards: cluster_to_fasta[wildcards.cluster],
        min_reads = instrain_filtered_reads
    shell: """
        f=$(basename {params.fasta_search} | sed 's/.fa//g'| sed 's/.fna//g')
        # f="NC_017491.1"
        rc=$(grep $f {input.idxstats_agg} | cut -f 3)
        # echo $f
        echo $rc
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
        samtools flagstat {input} | grep "read1" | cut -f 1 -d " " > {output.readcounts}
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
        readcounts = join(outdir, "drep_alignment_comparison/cluster_bams/{cluster}/readcounts/{sample}.txt")
    output:
        scaffold = join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_scaffold_info.tsv"),
        # mapping = join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_mapping_info.tsv")
    threads: 4
    resources: 
        time = lambda wildcards, attempt: 4 ** attempt,
        mem=64
    params:
        fasta = lambda wildcards: cluster_to_fasta[wildcards.cluster],
        outdir = join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}"),
        logfile = join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/log/log.log"),
        drep_fasta_folder = join(drep_folder, "dereplicated_genomes"),
        min_reads = instrain_filtered_reads
    shell: """
        # skip if readcounts are less than the min reads
        rc="$(cat {input.readcounts})"
        if [ $rc -gt {params.min_reads} ]; then
            inStrain profile {input.bam} {params.fasta} -o {params.outdir} -p {threads} || true
            
            # catch error with zero read pairs remaining, and then fake output if thats the case
            if grep -q "no read pairs remain" {params.logfile}; then
                echo "NO FILTERED READS REMAIN, FAKING OUTPUT"
                touch {output}
            else
                inStrain genome_wide -i {params.outdir} -p {threads}
            fi
        else
            echo "# Not more than {params.min_reads} in the input bamfile, sample skipped" > {output.scaffold}
            fi
    """

            # echo "# Not more than {params.min_reads} in the input bamfile, sample skipped" > {output.mapping}