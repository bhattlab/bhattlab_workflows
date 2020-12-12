from os.path import join, abspath, expanduser, exists, basename
# running the instrain compare part of the pipeline, separate from the 
# profile part to speed things up
localrules: instrain_heatmaps_filtered, instrain_heatmaps_all, instrain_filter_reads, instrain_plot_filtered, instrain_genome_wide_filtered

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
drep_contigs = join(drep_folder, 'all_dereplicated_genomes.fa')
drep_fasta_dir = join(drep_folder, "dereplicated_genomes")
# sample groups added later, keep support for old configfiles
if sample_groups in config.keys():
    sample_groups = config['sample_groups']
else:
    sample_groups = None
###################################################################

# if no clusters file specified, will use the determination method from 
# this pipeline which looks at coverage
# print(cluster_file)
if cluster_file is None:
    cluster_file = join(outdir, "top_clusters.txt")

cluster_to_fasta = get_cluster_fasta(cluster_file)
cluster_list = list(cluster_to_fasta.keys())
sample_reads = get_sample_reads(sample_file)
sample_list = list(sample_reads.keys())

# limit clusters if desired
if limit_clusters is not None:
    limit_list = [int(a) for a in str(limit_clusters).split(",")]
    cluster_list = [cluster_list[a] for a in limit_list]

# minimum number of filtered reads for a sample to be included in the "filtered"
# instrain version
instrain_min_reads = config['instrain_min_reads']


print()
print("#########################################")
print("# cluster list: " + str(cluster_list))
print("# instrain_min_reads: " + str(instrain_min_reads))
print("#########################################")
print()
print(' ONLY CONDUCTING INSTRAIN **COMPARE** STEPS IN THIS PIELINE ')
print(' ON SAMPLES THAT HAVE AT LEAST ' + str(instrain_min_reads) + ' READS AFTER FILTERING')
print()
 
rule all:
    input:
        expand(join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/{cluster}_filtered_samples.txt"), cluster=cluster_list),
        expand(join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/figures/{cluster}_inStrainCompare_dendrograms.pdf"), cluster=cluster_list),
        expand(join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/heatmaps/popANI/popANI_heatmap_unfiltered_complete.pdf"), cluster=cluster_list),
        expand(join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/output/{cluster}_comparisonsTable.tsv", cluster=cluster_list)),
        join(outdir, "drep_alignment_comparison/instrain_compare_filtered/instrain_compare_compiled.tsv")
        # expand(join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/figures/{cluster}_inStrainCompare_dendrograms.pdf"), cluster=cluster_list),
        # expand(join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/heatmaps/popANI/popANI_heatmap_unfiltered_complete.pdf"), cluster=cluster_list),

rule instrain_compare_all:
    input:
        lambda wildcards: expand(join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_scaffold_info.tsv"),
            sample=sample_list, cluster=wildcards.cluster)
    output:
        join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/output/{cluster}_comparisonsTable.tsv")
    threads: 16
    params: 
        files = lambda wildcards: expand(join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}"),
            sample=sample_list, cluster=wildcards.cluster),
        outdir = join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}"),
        threads_actual = 16 
    shell: """
        inStrain compare -i {params.files} -o {params.outdir} -p {params.threads_actual}
    """

rule instrain_genome_wide_all:
    input:
        join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/output/{cluster}_comparisonsTable.tsv")
    output:
        join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/output/{cluster}_genomeWide_compare.tsv")
    params: 
        fasta_name = lambda wildcards: cluster_to_fasta[wildcards.cluster],
        fasta_dir = drep_fasta_dir,
        outdir = join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}")
    shell: """
        inStrain genome_wide -i {params.outdir} -s {params.fasta_name}
    """

rule instrain_plot_all:
    input:
        join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/output/{cluster}_genomeWide_compare.tsv")
    output:
        join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/figures/{cluster}_inStrainCompare_dendrograms.pdf")
    params: 
        outdir = join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}")
    shell: """
        inStrain plot -i {params.outdir}
    """

rule instrain_filter_reads:
    input:
        lambda wildcards: expand(join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_mapping_info.tsv"),
            sample=sample_list, cluster=wildcards.cluster)
    output:
        filtered_samples_file = join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/{cluster}_filtered_samples.txt"),
        all_samples_file = join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/{cluster}_all_samples.txt"),
    params: 
        tail_string = join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/*/output/*_genomeWide_read_report.tsv"),
        outdir = join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}"),
        min_reads = instrain_min_reads,
        filtered_reads_file = join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/{cluster}_filtered_reads.txt"),
    run:
        filtered_reads = {}
        for s in sample_list:
            f = join(outdir, "drep_alignment_comparison/cluster_instrain/", 
                wildcards.cluster, s, "output", s + "_mapping_info.tsv")
            # print(f)
            with open(f, 'r') as inf:
                lines = [l for l in inf.readlines() if l[0]!="#"]
                # print(lines)  
                if len(lines) < 2:
                    filtered_reads[s] = '0'
                else: 
                    line = lines[1].strip().split('\t')
                    filtered_reads[s] = line[2]
                    # print(filtered_reads[s])
        # write out the results if they pass the threshold
        outf = output[0]
        with open(outf, 'w') as outf:
            for k,v in filtered_reads.items():
                if (int(float(v)) >= instrain_min_reads):
                    outf.write('\t'.join([k,v]) + '\n')
        # write out all samples
        outf = output[1]
        with open(outf, 'w') as outf:
            for k,v in filtered_reads.items():
                outf.write('\t'.join([k,v]) + '\n')

                
rule instrain_compare_filtered:
    input:
        lambda wildcards: expand(join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/{sample}/output/{sample}_scaffold_info.tsv"),
            sample=sample_list, cluster=wildcards.cluster),
        filtered_samples_file = join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/{cluster}_filtered_samples.txt"),
    output:
        join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/output/{cluster}_comparisonsTable.tsv")
    threads: 8
    resources:
        time= 72,
        mem=256
    params: 
        outdir = join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}"),
        sed_string = join(outdir, "drep_alignment_comparison/cluster_instrain/{cluster}/"),
    shell: """
        if [[ $(wc -l <{input.filtered_samples_file}) -ge 2 ]]; then
            files="$(cut -f1 {input.filtered_samples_file} | sed "s#^#{params.sed_string}#g" | tr "\n" " ")"
            # echo $files
            inStrain compare -i $files -o {params.outdir} -p {threads} --store_mismatch_locations
        else
            # not enough samples to run this pipeline
            echo "NOT ENOUGH FILTERD SAMPLES PRESENT" > {output}
        fi
    """

rule instrain_genome_wide_filtered:
    input:
        join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/output/{cluster}_comparisonsTable.tsv")
    output:
        join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/output/{cluster}_genomeWide_compare.tsv")
    params: 
        fasta_name = lambda wildcards: cluster_to_fasta[wildcards.cluster],
        fasta_dir = drep_fasta_dir,
        outdir = join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}")
    shell: """
        if [[ $(wc -l <{input}) -ge 2 ]]; then
            inStrain genome_wide -i {params.outdir} -s {params.fasta_name}
        else
            # not enough samples to run this pipeline
            echo "NOT ENOUGH FILTERD SAMPLES PRESENT" > {output}
        fi

    """

rule instrain_plot_filtered:
    input:
        join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/output/{cluster}_genomeWide_compare.tsv")
    output:
        join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/figures/{cluster}_inStrainCompare_dendrograms.pdf")
    params: 
        outdir = join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}")
    shell: """
        if [[ $(wc -l <{input}) -ge 2 ]]; then
            inStrain plot -i {params.outdir}
        else
            # not enough samples to run this pipeline
            echo "NOT ENOUGH FILTERD SAMPLES PRESENT" > {output}
        fi
    """

# custom heatmaps
rule instrain_heatmaps_all:
    input:
        rules.instrain_genome_wide_all.output
    output: 
        join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/heatmaps/popANI/popANI_heatmap_unfiltered_complete.pdf")
    params: 
        outdir = join(outdir, "drep_alignment_comparison/instrain_compare_all/{cluster}/heatmaps/"),
        min_frac_compared = 0.2, 
        min_identity_plot = 0,
        cluster_name = lambda wildcards: wildcards.cluster
    script: "scripts/heatmaps_instrain.R"

rule instrain_heatmaps_filtered:
    input:
        rules.instrain_genome_wide_filtered.output
    output: 
        join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/heatmaps/popANI/popANI_heatmap_unfiltered_complete.pdf")
    params: 
        outdir = join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/heatmaps/"),
        min_frac_compared = 0.2, 
        min_identity_plot = 0,
        cluster_name = lambda wildcards: wildcards.cluster
    script: "scripts/heatmaps_instrain.R"

# compilation of tables from all clusters
rule instrain_compile_tables:
    input:
        expand(join(outdir, "drep_alignment_comparison/instrain_compare_filtered/{cluster}/output/{cluster}_comparisonsTable.tsv", cluster=cluster_list)
    output:
        join(outdir, "drep_alignment_comparison/instrain_compare_filtered/instrain_compare_compiled.tsv")
    params:
        outdir: join(outdir, "drep_alignment_comparison/instrain_compare_filtered/"),
        sample_groups: sample_groups,
    scripts: "scripts/compile_compare_tables.R"
