#!/usr/bin/env python
from os.path import join, abspath, expanduser, exists
import sys
localrules: bwa_index_setup, postprocess, label_bins, extract_DAStool, concoct_extract_bins

# TODO
# bin_idxstats and bin_coverage need to be done with average of group
# is there a threshold in das tool I'm not aware of yet

def get_sample_group_reads(sample_file):
    sample_reads = {}
    group_reads1 = {}
    group_reads2 = {}
    group_to_samples = {}
    group_assemblies = {}
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            sample = s[0]
            if sample in sample_reads:
                print(sample)
                raise ValueError("Non-unique sample encountered!")
            group = s[1]
            assembly = s[2]
            reads_split = s[3].split(',')
            # get read pairs and singles from read specification
            if (len(reads_split) == 3) or (len(reads_split) == 2):
                group_reads1.setdefault(group, []).append(reads_split[0])
                group_reads2.setdefault(group, []).append(reads_split[1])
            else: 
                sys.exit('Bad read specification')
            # add these to dictionaries
            sample_reads[sample] = reads_split
            group_to_samples.setdefault(group, []).append(sample)
            if group in group_assemblies:
                if group_assemblies[group] != assembly:
                    sys.exit('Different assemblies specified for the same group')
            else: 
                group_assemblies[group] = assembly

    return (sample_reads, group_reads1, group_reads2, group_to_samples, group_assemblies)

# Read in sample and outdir from config file
outdir = config['outdir_base']
sample_file = config['sample_file']
sample_reads, group_reads1, group_reads2, group_to_samples, group_assemblies = get_sample_group_reads(sample_file)
group_list = list(group_to_samples.keys())
sample_list = list(sample_reads.keys())
# convert outdir to absolute path
if outdir[0] == '~':
    outdir = expanduser(outdir)
outdir = abspath(outdir)

# how are we going to do this: 
#   map all sample reads file to the same contig set
#   

# create set of output bams that we need to map
# each sample to the set of contigs
bam_align_files = []
for group in group_list:
    for sample in group_to_samples[group]:
        bam_align_files.append(join(outdir, group, sample + ".bam"))


# are we using a non-standard (non ncbi) taxonomy
if 'custom_taxonomy' in config:
    custom_taxonomy = config['custom_taxonomy']
else:
    custom_taxonomy = False

# Determine if long reads
if 'long_read' in config and config['long_read']:
    long_read = True
    sys.exit('Untested in this implementation')
else:
    long_read = False


# to speedup execution - remove samples that are done completely
skip_finished = True
if skip_finished:
    # sample_reads_new = [s for s in sample_reads if (!(exists(join(outdir, s, "final" + s + ".tsv")))) ]
    group_list_new = [s for s in group_list if not exists(join(outdir, s, "final", s + ".tsv")) ]
    group_list = group_list_new

print('##################################################################')
print(' GROUP LIST ')
print(group_list)
print('##################################################################')

def get_metabat_bins(wildcards):
    outputs = checkpoints.metabat.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(outputs, "{metabat_bin}.fa")).metabat_bin

def get_maxbin_bins(wildcards):
    outputs = checkpoints.maxbin.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(outputs, "{maxbin_bin}.fasta")).maxbin_bin

def get_mycc_bins(wildcards):
    outputs = checkpoints.mycc.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(outputs, "{mycc_bin}.fasta")).mycc_bin

def get_concoct_bins(wildcards):
    outputs = checkpoints.concoct_extract_bins.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(outputs, "{concoct_bin}.fasta")).concoct_bin

def get_DAStool_bins(wildcards):
    outputs = checkpoints.extract_DAStool.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(outputs, "{bin}.fa")).bin

rule all:
    input:
        config['kraken2db'],
        bam_align_files,
        expand(join(outdir, "{group}/{group}.fa.depth.txt"), group=group_list),
        # checkm - not necessary to run for all binners
        # expand(join(outdir, "{sample}/metabat/checkm/checkm.tsv"), sample = config['sample']),
        # expand(join(outdir, "{sample}/maxbin/checkm/checkm.tsv"), sample = config['sample']),
        # expand(join(outdir, "{sample}/mycc/checkm/checkm.tsv"), sample = config['sample']),
        # Post-processing
        expand(join(outdir, "{group}/classify/bin_species_calls.tsv"), group = group_list),
        expand(join(outdir, "{group}/final/{group}.tsv"), group = group_list),
        expand(join(outdir, "{group}/final/{group}_simple.tsv"), group = group_list),
        join(outdir, "binning_table_all_full.tsv")

##########################################
####### Prepping input for binning #######
##########################################

# copy assembly fasta file to working directory
rule bwa_index_setup:
    input:
        lambda wildcards: group_assemblies[wildcards.group]
    output:
        join(outdir, "{group}/idx/{group}.fa")
    resources:
        mem=1,
        time=1
    threads: 1
    shell: """
        cp {input} {output}
        """

# index assembly file
rule bwa_index:
    input:
        join(outdir, "{group}/idx/{group}.fa")
    output:
        join(outdir, "{group}/idx/{group}.fa.amb"),
        join(outdir, "{group}/idx/{group}.fa.ann"),
        join(outdir, "{group}/idx/{group}.fa.bwt"),
        join(outdir, "{group}/idx/{group}.fa.pac"),
        join(outdir, "{group}/idx/{group}.fa.sa")
    log:
        join(outdir, "{group}/logs/bwa_index.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem=8,
        time=2
    threads: 1
    shell: """
        bwa index {input}
        """

# Align reads to assembly
rule bwa_align:
    input:
        fwd = lambda wildcards: sample_reads[wildcards.sample][0],
        rev = lambda wildcards: sample_reads[wildcards.sample][1],
        asm = join(outdir, "{group}/idx/{group}.fa"),
        amb = join(outdir, "{group}/idx/{group}.fa.amb"),
        ann = join(outdir, "{group}/idx/{group}.fa.ann"),
        bwt = join(outdir, "{group}/idx/{group}.fa.bwt"),
        pac = join(outdir, "{group}/idx/{group}.fa.pac"),
        sa = join(outdir, "{group}/idx/{group}.fa.sa")
    output:
        join(outdir, "{group}/{sample}.bam")
    log:
        join(outdir, "{group}/logs/{sample}_bwa_mem.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem=16,
        time=12
    threads: 8
    shell: """
        bwa mem -t {threads} {input.asm} {input.fwd} {input.rev} | samtools sort --threads {threads} > {output}
        samtools index {output}
        """
'''
rule align_lr:
    input:
        join(outdir, "{sample}/idx/{sample}.fa"),
        reads
    log:
        join(outdir, "{sample}/logs/align_lr.log")
    output:
        join(outdir, "{sample}/{sample}_lr.bam")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem=48,
        time=6
    threads: 16
    shell: """
        minimap2 -t {threads} -ax map-ont {input} | samtools sort --threads {threads} > {output}
        """
'''

# Generate a depth file from BAM file for MetaBat input
rule metabat_pre:
    input:
        lambda wildcards: expand(join(outdir, "{group}/{sample}.bam"), 
            group=wildcards.group, sample=group_to_samples[wildcards.group])
    output:
        single = join(outdir, "{group}/{group}.fa.depth.txt"),
        paired = join(outdir, "{group}/{group}.fa.paired.txt"),
    singularity:
        "docker://quay.io/biocontainers/metabat2:2.15--h137b6e9_0"
    shell: """
        jgi_summarize_bam_contig_depths --outputDepth {output.single} \
        --pairedContigs {output.paired} --minContigLength 1000 \
        --minContigDepth 1 --percentIdentity 50 \
        {input} 
        """

#####################################################
################ Binning methods ####################
#####################################################

# Run MetaBat binner
checkpoint metabat:
    input:
        asm = join(outdir, "{group}/idx/{group}.fa"),
        depth = join(outdir, "{group}/{group}.fa.depth.txt"),
    output:
        directory(join(outdir, "{group}/metabat/bins/"))
    singularity:
        "docker://quay.io/biocontainers/metabat2:2.15--h137b6e9_0"
    resources:
        mem=64,
        time=24
    threads: 4
    params:
        outstring = join(outdir, "{group}/metabat/bins/bin")
    shell: """
        metabat2 --seed 1 -t {threads} --unbinned --inFile {input.asm} --outFile {params.outstring} --abdFile {input.depth}

        # if no bins produced, copy contigs to bin.unbinned
        if [ $(ls {output} | wc -l ) == "0" ]; then
            cp {input.asm} {output}/bin.unbinned.fa

        # check for bin.tooShort.fa thats empty
        if [ -f {output}/bin.tooShort.fa ]; then
            echo "Found bin.tooShort.fa"
            if [ $(cat {output}/bin.tooShort.fa | wc -l ) == "0" ]; then
                echo "Removing bin.tooShort.fa"
                rm {output}/bin.tooShort.fa
            fi
        fi
        """

# Run MaxBin2 binner
checkpoint maxbin:
    input:
        contigs = join(outdir, "{group}/idx/{group}.fa"),
        depth = join(outdir, "{group}/{group}.fa.depth.txt"),
    output:
        directory(join(outdir, "{group}/maxbin/bins"))
    singularity:
        "docker://quay.io/biocontainers/maxbin2:2.2.7--he1b5a44_1"
    params:
        outfolder = join(outdir, "{group}/maxbin/"),
        logfile = join(outdir, "{group}/maxbin/maxbin.log"),
        abundance_folder = join(outdir, "{group}/maxbin/depth_files"),
        abundance_list = join(outdir, "{group}/maxbin/abundance_list.txt"),
        num_samples = lambda wildcards: len(group_to_samples[wildcards.group]),
    resources:
        time=lambda wildcards, attempt: attempt * 2
    threads: 8
    shell: """
        if [ -d {params.outfolder} ]; then rm -r {params.outfolder}; fi
        mkdir -p {params.outfolder}
        cd {params.outfolder}
        # create abundance file, which we already have from the metabat rule
        # need to get out just the relevant columns, which is every
        # even column after 4
        mkdir {params.abundance_folder}
        for i in $(seq 1 {params.num_samples}); do
            col=$((2 + i * 2))
            cut -f 1,$col {input.depth} | tail -n +2 > {params.abundance_folder}/$i.txt
        done
        # make a list of abundance files
        ls {params.abundance_folder}/*.txt > {params.abundance_list}
        run_MaxBin.pl -contig {input.contigs} -out maxbin \
        -abund_list {params.abundance_list} -thread {threads} || true

        # look in log  if dataset cant be binned, just copy contigs over
        if $(grep -q "cannot be binned" {params.logfile}); then
            echo 'DATASET CANNOT BE BINNED'
            mkdir {params.outfolder}/bins/
            cp {input.contigs} {params.outfolder}/bins/maxbin.unbinned.fasta
        elif ls {params.outfolder}/*.fasta 1> /dev/null 2>&1; then
            echo "FOUND BINS"
            mkdir {params.outfolder}/bins/
            mv {params.outfolder}/*.fasta {params.outfolder}/bins/
        else 
            echo "PROGRAM MUST HAVE FAILED"
            exit 1
        fi

        """

# Run myCC binner
checkpoint mycc:
    input:
        contigs = join(outdir, "{group}/idx/{group}.fa"),
    output:
        directory(join(outdir, "{group}/mycc/bins"))
    singularity:
        "docker://990210oliver/mycc.docker:v1"
    params:
        outfolder = join(outdir, "{group}/mycc"),
        workdir = join(outdir, "{group}")
    resources:
        time=lambda wildcards, attempt: attempt * 12,
        mem=lambda wildcards, attempt: attempt * 128
    shell: """
        if [ -d {params.outfolder} ]; then rm -r {params.outfolder}; fi
        cd {params.workdir}
        MyCC.py {input.contigs} -meta
        mv 20*/ {params.outfolder}/  # change variable folder name to mycc
        mkdir {params.outfolder}/bins/
        mv {params.outfolder}/*.fasta {params.outfolder}/bins/
        """

# Run concoct binner in multiple steps
rule concoct_cut:
    input:
        join(outdir, "{group}/idx/{group}.fa"),
    output:
        cut_contigs = join(outdir, '{group}/concoct/{group}_contigs_10K.fa'),
        bedfile = join(outdir, '{group}/concoct/{group}_contigs_10K.bed')
    resources:
        mem = 64,
        time = 12
    singularity:
        'docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0'
    shell: """
        cut_up_fasta.py {input} -c 10000 -o 0 --merge_last -b {output.bedfile} > {output.cut_contigs}
    """

rule concoct_coverage:
    input:
        bedfile = rules.concoct_cut.output.bedfile,
        depth = join(outdir, "{group}/{group}.fa.depth.txt"),
    output:
        join(outdir, '{group}/concoct/coverage_table.tsv')
    resources:
        mem = 64,
        time = 12
    params:
        bam_search = join(outdir, "{group}/*.bam")
    singularity:
        'docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0'
    shell: """
        concoct_coverage_table.py {input.bedfile} {params.bam_search} > {output}
    """

rule concoct_run:
    input:
        cut_contigs = rules.concoct_cut.output.cut_contigs,
        coverage = rules.concoct_coverage.output
    output:
        join(outdir, '{group}/concoct/clustering_gt1000.csv')
    params:
        outdir = join(outdir, '{group}/concoct')
    threads: 4
    resources:
        mem = 64, 
        time = 24
    singularity:
        'docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0'
    shell: """
        concoct --composition_file {input.cut_contigs} \
        --coverage_file {input.coverage} -b {params.outdir} --threads {threads}
    """

rule concoct_merge:
    input:
        clustering = rules.concoct_run.output
    output:
        join(outdir, '{group}/concoct/clustering_merged.csv')
    resources:
        mem = 64,
        time = 12
    singularity:
        'docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0'
    shell: """
        merge_cutup_clustering.py {input} > {output}
    """

checkpoint concoct_extract_bins:
    input:
        original_contigs = join(outdir, "{group}/idx/{group}.fa"),
        clustering_merged = rules.concoct_merge.output
    output:
        directory(join(outdir, "{group}/concoct/bins/")) #the number of bins is unknown prior to execution
    params:
        outdir = join(outdir, '{group}/concoct/bins/')
    resources:
        mem = 64,
        time = 12
    singularity:
        'docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0'
    shell: """
    extract_fasta_bins.py {input.original_contigs} {input.clustering_merged} \
    --output_path {params.outdir}
    """

# Aggregate binning results using DAStool
rule DAStool:
    input:
        lambda wildcards: expand(join(outdir, "{group}/metabat/bins/{metabat_bin}.fa"), metabat_bin = get_metabat_bins(wildcards), group = wildcards.group),
        lambda wildcards: expand(join(outdir, "{group}/maxbin/bins/{maxbin_bin}.fasta"), maxbin_bin = get_maxbin_bins(wildcards), group = wildcards.group),
        # lambda wildcards: expand(join(outdir, "{group}/mycc/bins/{mycc_bin}.fasta"), mycc_bin = get_mycc_bins(wildcards), group = wildcards.group),
        lambda wildcards: expand(join(outdir, "{group}/concoct/bins/{concoct_bin}.fasta"), concoct_bin = get_concoct_bins(wildcards), group = wildcards.group),
        contigs = join(outdir, "{group}/idx/{group}.fa"),
    output: 
        join(outdir, "{group}/DAS_tool/completed.txt")
    # singularity:
        # "docker://quay.io/biocontainers/das_tool:1.1.1--py36r351_1"
    conda: "envs/das_tool.yaml"
    params:
        outfolder = join(outdir, "{group}/DAS_tool/"),
        metabat_dir = join(outdir, "{group}/metabat/bins/"),
        maxbin_dir = join(outdir, "{group}/maxbin/bins/"),
        mycc_dir = join(outdir, "{group}/mycc/bins/"),
        concoct_dir = join(outdir, "{group}/concoct/bins/"),
        metabat_tsv = join(outdir, "{group}/DAS_tool/metabat_scaffold2bin.tsv"),
        maxbin_tsv = join(outdir, "{group}/DAS_tool/maxbin_scaffold2bin.tsv"),
        mycc_tsv = join(outdir, "{group}/DAS_tool/mycc_scaffold2bin.tsv"),
        concoct_tsv = join(outdir, "{group}/DAS_tool/concoct_scaffold2bin.tsv"),
        logfile = join(outdir, "{group}/DAS_tool/_DASTool.log"),
    threads: 8
    resources:
        time=lambda wildcards, attempt: attempt * 6
    shell: """
        # Prepare scaffold2bin file for each set of bins
        Fasta_to_Scaffolds2Bin.sh -e fa -i {params.metabat_dir} > {params.metabat_tsv}
        Fasta_to_Scaffolds2Bin.sh -e fasta -i {params.maxbin_dir} > {params.maxbin_tsv}
        Fasta_to_Scaffolds2Bin.sh -e fa -i {params.concoct_dir} > {params.concoct_tsv}
        # Fasta_to_Scaffolds2Bin.sh -e fasta -i {params.mycc_dir} > {params.mycc_tsv}
        # Run DAS_Tool
        DAS_Tool -i {params.metabat_tsv},{params.maxbin_tsv},{params.concoct_tsv} \
        -l metabat,maxbin,concoct -c {input.contigs} -o {params.outfolder} \
        --search_engine diamond --threads {threads} --write_bins 1 --write_unbinned 1 || true

        if $(grep -q "No bins with bin-score >0.5 found" {params.logfile}) || $(grep -q "single copy gene prediction using diamond failed" {params.logfile}); then
            echo 'DATASET CANNOT BE BINNED'
            mkdir {params.outfolder}/_DASTool_bins
            cp {input.contigs} {params.outfolder}/_DASTool_bins/unbinned.fa
            touch {output}
        elif ls {params.outfolder}/_DASTool_bins/*.fa 1> /dev/null 2>&1; then
            echo "FOUND BINS"
            touch {output}
        else 
            echo "PROGRAM MUST HAVE FAILED IN SOME OTHER WAY"
            exit 1
        fi

        """

# extract the bins as a separate rule to speed up execution
checkpoint extract_DAStool:
    input: rules.DAStool.output
    output:
        directory(join(outdir, "{group}/DAS_tool_bins"))
    params:
        old_binfolder = join(outdir, "{group}/DAS_tool/_DASTool_bins"),
        new_binfolder = join(outdir, "{group}/DAS_tool_bins"),
    shell: """
        cp -r {params.old_binfolder} {params.new_binfolder}
    """


#####################################################
###################### CheckM #######################
#####################################################

# checkm for metabat output
rule checkm_metabat:
    input:
        lambda wildcards: expand(join(outdir, "{group}/metabat/bins/{metabat_bin}.fa"), metabat_bin = get_metabat_bins(wildcards), group = wildcards.group)
    output:
        join(outdir, "{group}/metabat/checkm/checkm.tsv")
    log:
        join(outdir, "{group}/metabat/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{group}/metabat/bins/"),
        checkmfolder = join(outdir, "{group}/metabat/checkm/"),
    shell: """
        rm -rf {group}/metabat/checkm/*
        checkm lineage_wf -t {threads} -x fa --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

# checkm for maxbin output
rule checkm_maxbin:
    input:
        lambda wildcards: expand(join(outdir, "{group}/maxbin/bins/{maxbin_bin}.fasta"), maxbin_bin = get_maxbin_bins(wildcards), group = wildcards.group)
    output:
        join(outdir, "{group}/maxbin/checkm/checkm.tsv")
    log:
        join(outdir, "{group}/maxbin/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{group}/maxbin/bins"),
        checkmfolder = join(outdir, "{group}/maxbin/checkm/")
    shell: """
        rm -rf {group}/maxbin/checkm/*
        checkm lineage_wf -t {threads} -x fasta --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

rule checkm_mycc:
    input:
        lambda wildcards: expand(join(outdir, "{group}/mycc/bins/{mycc_bin}.fasta"), mycc_bin = get_mycc_bins(wildcards), group = wildcards.group)
    output:
        join(outdir, "{group}/mycc/checkm/checkm.tsv")
    log:
        join(outdir, "{group}/mycc/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{group}/mycc/bins/"),
        checkmfolder = join(outdir, "{group}/mycc/checkm/"),
    shell: """
        rm -rf {group}/mycc/checkm/*
        checkm lineage_wf -t {threads} -x fasta --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

# checkm for DAStool output
rule checkm_DAStool:
    input:
        lambda wildcards: expand(join(outdir, "{group}/DAS_tool_bins/{bin}.fa"), bin = get_DAStool_bins(wildcards), group = wildcards.group)
    output:
        join(outdir, "{group}/DAS_tool/checkm/checkm.tsv")
    log:
        join(outdir, "{group}/DAS_tool/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{group}/DAS_tool_bins"),
        checkmfolder = join(outdir, "{group}/DAS_tool/checkm/"),
        bin_ex = ".fa"
    shell: """
        rm -rf {params.checkmfolder}
        checkm lineage_wf -t {threads} -x {params.bin_ex} --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

#####################################################
############ ANALYSIS OF DAS_TOOL BINS ##############
#####################################################

# Use aragorn to detect tRNA and tmRNA genes
rule aragorn:
    input:
        join(outdir, "{group}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{group}/rna/trna/{bin}.fa.txt")
    log:
        join(outdir, "{group}/logs/aragorn_{bin}.log")
    singularity:
        "docker://quay.io/biocontainers/prokka:1.14.5--pl526_0"
    resources:
        mem = 8,
        time = 1
    shell: """
        aragorn -t {input} -o {output}
        """

# Use barrnap to predict ribosomal RNA
rule barrnap:
    input:
        join(outdir, "{group}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{group}/rna/rrna/{bin}.fa.txt")
    log:
        join(outdir, "{group}/logs/barrnap_{bin}.log")
    singularity:
        "docker://quay.io/biocontainers/prokka:1.14.5--pl526_0"
    resources:
        mem = 8,
        time = 1
    shell: """
        barrnap {input} > {output}
        """

# Use quast to evaluate genome assemblies
rule quast:
    input:
        join(outdir, "{group}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{group}/quast/{bin}.fa/report.tsv")
    log:
        join(outdir, "{group}/logs/quast_{bin}.log")
    singularity:
        "docker://quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0"
    resources:
        mem=8,
        time=1
    params:
        quastfolder = join(outdir, "{group}/quast/{bin}.fa/"),
        thresholds = "0,10000,50000,100000,250000,500000,1000000,2000000,3000000"
    shell: """
        quast.py -o {params.quastfolder} {input} --contig-thresholds {params.thresholds} --fast
        """

rule pull_prokka:
    input: join(outdir, "{group}/idx/{group}.fa"),
    output: join(outdir, "{group}/prokka/pulled.txt") 
    singularity:
        "docker://quay.io/biocontainers/prokka:1.14.5--pl526_0"
    shell: """
        touch {output}
    """

# Use prokka to annotate metagenomes
rule prokka:
    input:
        join(outdir, "{group}/DAS_tool_bins/{bin}.fa"),
        rules.pull_prokka.output
    output:
        join(outdir, "{group}/prokka/{bin}.fa/{group}_{bin}.fa.gff")
    log:
        join(outdir, "{group}/logs/prokka_{bin}.log")
    singularity:
        "docker://quay.io/biocontainers/prokka:1.14.5--pl526_0"
    resources:
        mem = 48,
        time = lambda wildcards, attempt: 4 * attempt,
    threads: 8
    params:
        prokkafolder = join(outdir, "{group}/prokka/{bin}.fa"),
        prefix = "{group}_{bin}.fa"
    shell: """
        # don't run this on unbinned contigs, takes forever
        if [[ {wildcards.bin} =~ "unbinned" ]]; then
            touch {output}
            touch {params.prokkafolder}/prokka_skipped.out
        else
            prokka {input} --outdir {params.prokkafolder} --prefix {params.prefix} \
            --centre X --compliant --force --cpus {threads} --noanno
        fi
        """

# Index bam file to prepare for bamidx
rule bam_idx:
    input:
        join(outdir, "{group}/{group}_lr.bam") if long_read else join(outdir, "{group}/{group}.bam") #choose a long read alignment or short read alignment
    output:
        join(outdir, "{group}/{group}_lr.bam.bai") if long_read else join(outdir, "{group}/{group}.bam.bai") #choose a long read alignment or short read alignment
    log:
        join(outdir, "{group}/logs/bamidx.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem = 2,
        time = 2
    shell:
        "samtools index {input}"

# Retrieve stats on mapping from whole metagenome group
rule bam_idxstats:
    input:
        join(outdir, "{group}/{group}.bam"),
        join(outdir, "{group}/{group}.bam.bai"),
    output:
        join(outdir, "{group}/{group}.bam.bai.tsv"),
    log:
        join(outdir, "{group}/logs/bamidxstats.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem = 2,
        time = 2
    shell:
        "samtools idxstats {input[0]} > {output}"

# do this with the depth file from metabat now
rule bin_idxstats:
    input:
        bin_fasta = join(outdir, "{group}/DAS_tool_bins/{bin}.fa"),
        depth = rules.metabat_pre.output.single
    output:
        join(outdir, "{group}/coverage/raw/{bin}.tsv")
    log:
        join(outdir, "{group}/logs/coverage_idxstats_{bin}.log")
    resources:
        mem = 2,
        time = 6
    shell: """
        # some dangerous execution here because ths was failing on unbinned...
        set +e; set +o pipefail; grep '>' {input.bin_fasta} | tr -d '>' | cut -d " " -f 1 | xargs -I foo -n 1 grep -P 'foo\t' {input.depth} > {output} || true
        """
# Determine coverage for each bin
rule bin_coverage:
    input:
        rules.bin_idxstats.output
    output:
        join(outdir, "{group}/coverage/{bin}.txt")
    log:
        join(outdir, "{group}/logs/coverage_{bin}.log")
    resources:
        mem = 2,
        time = 1
    params:
        read_length = config['read_length'],
        sample = lambda wildcards: wildcards.group
    run:
        length = 0
        coverage = 0
        with open(input[0], 'r') as inf:
            for l in inf.readlines():
                ls = l.strip().split('\t')
                length += float(ls[1])
                coverage += float(ls[1]) * float(ls[2])
        final_coverage = coverage /length
        with open(output[0], 'w') as outf:
            outf.write('\t'.join([params['sample'],
                wildcards['bin'], str(final_coverage)]) + '\n')
        # do a weighted average of the contigs 
        # "scripts/bin_coverage.py"

# index DAStool bins
rule fasta_index:
    input:
        join(outdir, "{group}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{group}/DAS_tool_bins/{bin}.fa.fai")
    log:
        join(outdir, "{group}/logs/faidx_{bin}.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem = 8,
        time = 1
    threads: 1
    shell:
        "samtools faidx {input}"

# Use kraken2 to assign taxonomic labels
rule kraken2:
    input:
        join(outdir, "{group}/idx/{group}.fa")
    output:
        krak = join(outdir, "{group}/classify/{group}.krak"),
        krak_report = join(outdir, "{group}/classify/{group}.krak.report")
    log:
        join(outdir, "{group}/logs/kraken_class.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    params:
        db = config['kraken2db']
    resources:
        mem = 256,
        time = 6
    threads: 16
    shell: """
        kraken2 --db {params.db} --db {params.db} --threads {threads} \
        --output {output.krak} --report {output.krak_report} {input}
        """

rule label_bins:
    input:
        krak = join(outdir, "{group}/classify/{group}.krak"),
        bins = lambda wildcards: expand(join(outdir, "{group}/DAS_tool_bins/{bin}.fa.fai"),
                                        bin = get_DAStool_bins(wildcards), group = wildcards.group)
    output:
        join(outdir, "{group}/classify/bin_species_calls.tsv")
    log:
        join(outdir, "{group}/logs/assign_species.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    params:
        binfolder = join(outdir, "{group}/DAS_tool_bins/"),
        custom_taxonomy = custom_taxonomy
    script:
        "scripts/assign_species.py"

rule postprocess:
    input:
        prokka = lambda wildcards: expand(join(outdir, "{group}/prokka/{bin}.fa/{group}_{bin}.fa.gff"), bin = get_DAStool_bins(wildcards), group = wildcards.group),
        quast = lambda wildcards: expand(join(outdir, "{group}/quast/{bin}.fa/report.tsv"), bin = get_DAStool_bins(wildcards), group = wildcards.group),
        checkm = join(outdir, "{group}/DAS_tool/checkm/checkm.tsv"),
        trna = lambda wildcards: expand(join(outdir, "{group}/rna/trna/{bin}.fa.txt"), bin = get_DAStool_bins(wildcards), group = wildcards.group),
        rrna = lambda wildcards: expand(join(outdir, "{group}/rna/rrna/{bin}.fa.txt"), bin = get_DAStool_bins(wildcards), group = wildcards.group),
        classify = rules.label_bins.output,
        coverage = lambda wildcards: expand(join(outdir, "{group}/coverage/{bin}.txt"), bin = get_DAStool_bins(wildcards), group = wildcards.group),
    output:
        full = join(outdir, "{group}/final/{group}.tsv"),
        simple = join(outdir, "{group}/final/{group}_simple.tsv")
    params:
        bins = lambda wildcards: get_DAStool_bins(wildcards),
        sample = lambda wildcards: wildcards.group
    script: "scripts/postprocess.R"

rule combine_final_reports:
    input:
        all_full = expand(join(outdir, "{group}/final/{group}.tsv"), group=group_list),
        single_full = expand(join(outdir, "{group}/final/{group}.tsv"), group=group_list[0]),
        all_simple = expand(join(outdir, "{group}/final/{group}_simple.tsv"), group=group_list),
        single_simple = expand(join(outdir, "{group}/final/{group}.tsv"), group=group_list[0]),
    output:
        full = join(outdir, "binning_table_all_full.tsv"),
        simple = join(outdir, "binning_table_all_simple.tsv"),
    shell: """
        head -n 1 {input.single_full} > {output.full}
        tail -n +2 -q {input.all_full} >> {output.full}
        head -n 1 {input.single_simple} > {output.simple}
        tail -n +2 -q {input.all_simple} >> {output.simple}
    """

