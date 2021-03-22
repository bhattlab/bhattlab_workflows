from os.path import join, abspath, expanduser, exists, basename
# isntrain compare for the version of the pipeline that uses 
# a single refrence genome
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
print(cluster_file)
limit_clusters=config['limit_clusters']
outdir=config['outdir']
barcodes=config['barcodes']
###################################################################

# if no clusters file specified, will use the determination method from 
# this pipeline which looks at coverage
if cluster_file is None:
    cluster_file = join(outdir, "top_clusters.txt")

cluster_to_fasta = get_cluster_fasta(cluster_file)
cluster_list = list(cluster_to_fasta.keys())
sample_reads = get_sample_reads(sample_file)
sample_list = list(sample_reads.keys())
ref_fasta = list(cluster_to_fasta.values())[0]
ref_cluster = cluster_list[0]

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
        expand(join(outdir, "instrain_compare_filtered/filtered_samples.txt")),
        expand(join(outdir, "instrain_compare_filtered/figures/inStrainCompare_dendrograms.pdf")),
        expand(join(outdir, "instrain_compare_filtered/heatmaps/popANI/popANI_heatmap_unfiltered_complete.pdf")),

rule instrain_filter_reads:
    input:
        lambda wildcards: expand(join(outdir, "instrain_profile/{sample}/output/{sample}_mapping_info.tsv"),
            sample=sample_list)
    output:
        filtered_samples_file = join(outdir, "instrain_compare_filtered/filtered_samples.txt"),
        all_samples_file = join(outdir, "instrain_compare_filtered/all_samples.txt"),
    params: 
        tail_string = join(outdir, "instrain_profile/*/output/*_genomeWide_read_report.tsv"),
        outdir = join(outdir, "instrain_compare_filtered"),
        min_reads = instrain_min_reads,
        filtered_reads_file = join(outdir, "instrain_compare_filtered/filtered_reads.txt"),
    run:
        filtered_reads = {}
        for s in sample_list:
            f = join(outdir, "instrain_profile", s, "output", s + "_mapping_info.tsv")
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
        lambda wildcards: expand(join(outdir, "instrain_profile/{sample}/output/{sample}_scaffold_info.tsv"),
            sample=sample_list),
        filtered_samples_file = join(outdir, "instrain_compare_filtered/filtered_samples.txt"),
    output:
        join(outdir, "instrain_compare_filtered/output/comparisonsTable.tsv")
    threads: 8
    resources:
        time= 48,
        mem=256
    params: 
        outdir = join(outdir, "instrain_compare_filtered"),
        sed_string = join(outdir, "instrain_profile/")
    singularity: "quay.io/biocontainers/instrain:1.5.2--py_0"
    shell: """
        if [[ $(wc -l <{input.filtered_samples_file}) -ge 2 ]]; then
            files="$(cut -f1 {input.filtered_samples_file} | sed "s#^#{params.sed_string}#g" | tr "\n" " ")"
            # echo $files
            inStrain compare -i $files -o {params.outdir} -p {threads} --store_mismatch_locations
            mv {params.outdir}/output/instrain_compare_filtered_comparisonsTable.tsv {params.outdir}/output/comparisonsTable.tsv
            mv {params.outdir}/output/instrain_compare_filtered_pairwise_SNP_locations.tsv {params.outdir}/output/pairwise_SNP_locations.tsv
        else
            # not enough samples to run this pipeline
            echo "NOT ENOUGH FILTERD SAMPLES PRESENT" > {output}
        fi
    """

rule instrain_genome_wide_filtered:
    input:
        join(outdir, "instrain_compare_filtered/output/comparisonsTable.tsv")
    output:
        join(outdir, "instrain_compare_filtered/output/genomeWide_compare.tsv")
    params: 
        fasta_name = join(outdir, "ref.fa"),
        outdir = join(outdir, "instrain_compare_filtered")
    singularity: "quay.io/biocontainers/instrain:1.5.2--py_0"
    shell: """
        if [[ $(wc -l <{input}) -ge 2 ]]; then
            inStrain genome_wide -i {params.outdir} -s {params.fasta_name}
            mv {params.outdir}/output/instrain_compare_filtered_genomeWide_compare.tsv {params.outdir}/output/genomeWide_compare.tsv
        else
            # not enough samples to run this pipeline
            echo "NOT ENOUGH FILTERD SAMPLES PRESENT" > {output}
        fi

    """

rule instrain_plot_filtered:
    input:
        join(outdir, "instrain_compare_filtered/output/genomeWide_compare.tsv")
    output:
        join(outdir, "instrain_compare_filtered/figures/inStrainCompare_dendrograms.pdf")
    params: 
        outdir = join(outdir, "instrain_compare_filtered")
    singularity: "quay.io/biocontainers/instrain:1.5.2--py_0"
    shell: """
        if [[ $(wc -l <{input}) -ge 2 ]]; then
            inStrain plot -i {params.outdir}
            mv {params.outdir}/figures/instrain_compare_filtered_inStrainCompare_dendrograms.pdf {output}
        else
            # not enough samples to run this pipeline
            echo "NOT ENOUGH FILTERD SAMPLES PRESENT" > {output}
        fi
    """

# custom heatmaps
rule instrain_heatmaps_filtered:
    input:
        rules.instrain_genome_wide_filtered.output
    output: 
        join(outdir, "instrain_compare_filtered/heatmaps/popANI/popANI_heatmap_unfiltered_complete.pdf")
    params: 
        outdir = join(outdir, "instrain_compare_filtered/heatmaps/"),
        min_frac_compared = 0.2, 
        min_identity_plot = 0,
        cluster_name = ref_cluster
    conda: "envs/r_processing.yaml"
    script: "scripts/heatmaps_instrain.R"

