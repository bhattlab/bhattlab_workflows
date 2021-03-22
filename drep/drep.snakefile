# Running dRep on the ouptut of DAS_Tool binning, automatically. 
# Compares all bins above a quality threshold together
# and outputs dereplicated genome set
from os.path import join, abspath, expanduser, exists, basename
localrules: rename_bins, filter_bin_quality, filter_bin_fasta, make_top_clusters

# Config file specifications 
outdir = config['outdir_base']
das_tool_folder = config['das_tool_folder']
min_completeness = config['min_completeness']
max_contamination = config['max_contamination']

# convert outdir to absolute path
if outdir[0] == '~':
    outdir = expanduser(outdir)
outdir = abspath(outdir)
# convert bin to absolute path
if das_tool_folder[0] == '~':
    das_tool_folder = expanduser(das_tool_folder)
das_tool_folder = abspath(das_tool_folder)


rule all:
    input:
        join(outdir, "all_dereplicated_genomes.fa"),
        join(outdir, 'top_clusters.txt')

# rename the bins as sample__bin.fa
# and create symlinks in a new folder
rule rename_bins:
    input:
        join(das_tool_folder, "binning_table_all_full.tsv")
    output:
        join(outdir, "symlinks_completed.txt")
    params:
        symlink_dir_unfiltered = join(outdir, "bins_das_tool_symlinks"),
        das_tool_folder = das_tool_folder
    shell: """
        # cleanup
        if [ -d {params.symlink_dir_unfiltered} ]; then
            rm -r {params.symlink_dir_unfiltered}
        fi
        mkdir {params.symlink_dir_unfiltered}

        cd  {das_tool_folder}
        for i in */; do 
            i=$(echo $i | sed 's#/##g')
            for j in $i/DAS_tool_bins/*.fa; do
                ln -s $(pwd -P)/"$j" {params.symlink_dir_unfiltered}/"$i"__$(basename $j)
            done
        done 
        touch {output}
    """

# filter for completeness and contamination 
# based on the standard binning output
# output quality metrics file
rule filter_bin_quality:
    input:
        join(das_tool_folder, "binning_table_all_full.tsv")
    output:
        qm = join(outdir, "quality_metrics.csv"),
        qmf = join(outdir, "quality_metrics_filtered.csv"),
    params:
        min_completeness=min_completeness, 
        max_contamination=max_contamination, 
    shell: """
        paste <(tail -n +2 {input} | cut -f 1,2 | sed "s/\t/__/g" | \
           sed "s/$/.fa/g") <(tail -n +2 {input} | \
           cut -f 5,6) | tr "\t" "," > {output.qm}
        awk -F"," '{{ if ($2 >={params.min_completeness} && $3 <={params.max_contamination}) print }}' {output.qm} > {output.qmf}
    """

# copy symlinks that meet quality metrics
rule filter_bin_fasta:
    input:
        quality = join(outdir, "quality_metrics_filtered.csv"),
        symlinks = join(outdir, "symlinks_completed.txt")
    output:
        join(outdir, "symlinks_filtered_completed.txt")
    params:
        symlink_dir_unfiltered = join(outdir, "bins_das_tool_symlinks"),
        symlink_dir_filtered = join(outdir, "bins_das_tool_symlinks_filtered")
    shell: """
        mkdir {params.symlink_dir_filtered}
        cut -f 1 -d"," {input.quality} | xargs -I{{}} cp -P {params.symlink_dir_unfiltered}/{{}} {params.symlink_dir_filtered}
        touch {output}
    """

#run drep, mostly default settings
rule drep:
    input:
        quality = join(outdir, "quality_metrics_filtered.csv"),
        symlinks = join(outdir, "symlinks_filtered_completed.txt")
    output:
        join(outdir, "all_dereplicated_genomes.fa")
    params:
        drep_outdir = join(outdir, "drep_actual"), 
        drep_genome_dir = join(outdir, "drep_actual", 'dereplicated_genomes'),
        symlink_dir_filtered = join(outdir, "bins_das_tool_symlinks_filtered"),
        min_completeness=min_completeness, 
        max_contamination=max_contamination, 
        secondary_clustering_algorithm = config['secondary_clustering_algorithm']
    threads: 16
    resources: 
        mem=128, 
        time=24
    singularity: "docker://quay.io/biocontainers/drep:2.6.2--py_0"
    shell: """
        if [ -d {params.drep_outdir}/data ]; then
            rm -rf {params.drep_outdir}/data
            rm -rf {params.drep_outdir}/data_tables
            rm -rf {params.drep_outdir}/dereplicated_genomes
            rm -rf {params.drep_outdir}/log
            rm -rf {params.drep_outdir}/figures
        fi

        dRep dereplicate {params.drep_outdir} \
            --ignoreGenomeQuality \
            -comp {params.min_completeness} \
            -con {params.max_contamination} \
            --genomeInfo {input.quality} \
            -g {params.symlink_dir_filtered}/*.fa \
            --S_algorithm {params.secondary_clustering_algorithm} \
            -p {threads}

        # add fasta headers to each from the genome
        cd {params.drep_genome_dir}
        for i in *.fa; do echo $i; sed -i "s/>/>"$i"__/g" $i; done
        cat *.fa > {output}
    """

# find the top clusters to use instrain on 
# that's the clusters with the most genomes in it
rule make_top_clusters:
    input:
        join(outdir, "all_dereplicated_genomes.fa")
    output:
        join(outdir, 'top_clusters.txt')
    params:
        outdir = outdir,
        keep_clusters = 20
    shell: """
        set +o pipefail
        cat {params.outdir}/drep_actual/data_tables/Cdb.csv | cut -f 2 -d "," | sort | uniq -c | sort -nr | head -n {params.keep_clusters} | tr -s " " | cut -f 3  -d " " > top_cluster_names.txt
        awk '{{ print ","$1"," }}' top_cluster_names.txt | grep {params.outdir}/data_tables/Wdb.csv -f - | cut -f 1 -d "," > top_fasta.txt
        paste top_cluster_names.txt  top_fasta.txt > {output}
    """
