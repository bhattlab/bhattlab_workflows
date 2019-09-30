#!/usr/bin/env python
# metagenomic binning with DASTool
# Written by Alyssa Benjamin (https://github.com/ambenj) June 2019

from os.path import join, abspath, expanduser

localrules: bwa_index_setup, postprocess, label_bins, extract_DAStool, concoct_extract_bins

# Read in sample and outdir from config file
samp = config['sample']
outdir = config['outdir_base']

# are we using a non-standard (non ncbi) taxonomy
if 'custom_taxonomy' in config:
    custom_taxonomy = config['custom_taxonomy']
else:
    custom_taxonomy = False

# convert outdir to absolute path
if outdir[0] == '~':
    outdir = expanduser(outdir)
outdir = abspath(outdir)

# Read in fastq files
if 'reads2' in config and not config['reads2'] == '':
    reads = [config['reads1'], config['reads2']]
else:
    reads = [config['reads1']]

# Determine if long reads
if 'long_read' in config and config['long_read']:
    long_read = True
else:
    long_read = False


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
    # print(glob_wildcards(os.path.join(outputs, "{bin}.fa")).bin)
    return glob_wildcards(os.path.join(outputs, "{bin}.fa")).bin

rule all:
    input:
        reads,
        config['assembly'],
        config['kraken2db'],
        # checkm - not necessary to run for all binners
        # expand(join(outdir, "{samp}/metabat/checkm/checkm.tsv"), samp = config['sample']),
        # expand(join(outdir, "{samp}/maxbin/checkm/checkm.tsv"), samp = config['sample']),
        # expand(join(outdir, "{samp}/mycc/checkm/checkm.tsv"), samp = config['sample']),
        # Post-processing
        expand(join(outdir, "{samp}/classify/bin_species_calls.tsv"), samp = config['sample']),
        expand(join(outdir, "{samp}/final/{samp}.tsv"), samp = config['sample']),
        expand(join(outdir, "{samp}/final/{samp}_simple.tsv"), samp = config['sample']),

##########################################
####### Prepping input for binning #######
##########################################

# copy assembly fasta file to working directory
rule bwa_index_setup:
    input:
        config['assembly']
    output:
        join(outdir, "{samp}/idx/{samp}.fa")
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
        join(outdir, "{samp}/idx/{samp}.fa")
    output:
        join(outdir, "{samp}/idx/{samp}.fa.amb"),
        join(outdir, "{samp}/idx/{samp}.fa.ann"),
        join(outdir, "{samp}/idx/{samp}.fa.bwt"),
        join(outdir, "{samp}/idx/{samp}.fa.pac"),
        join(outdir, "{samp}/idx/{samp}.fa.sa")
    log:
        join(outdir, "{samp}/logs/bwa_index.log")
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
        asm = join(outdir, "{samp}/idx/{samp}.fa"),
        reads = reads,
        amb = join(outdir, "{samp}/idx/{samp}.fa.amb"),
        ann = join(outdir, "{samp}/idx/{samp}.fa.ann"),
        bwt = join(outdir, "{samp}/idx/{samp}.fa.bwt"),
        pac = join(outdir, "{samp}/idx/{samp}.fa.pac"),
        sa = join(outdir, "{samp}/idx/{samp}.fa.sa")
    output:
        join(outdir, "{samp}/{samp}.bam")
    log:
        join(outdir, "{samp}/logs/bwa_mem.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem=16,
        time=12
    threads: 8
    shell: """
        bwa mem -t {threads} {input.asm}  {input.reads} |samtools sort --threads {threads} > {output}
        """

# Align long reads (not sure if I should keep this)
rule align_lr:
    input:
        join(outdir, "{samp}/idx/{samp}.fa"),
        reads
    log:
        join(outdir, "{samp}/logs/align_lr.log")
    output:
        join(outdir, "{samp}/{samp}_lr.bam")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem=48,
        time=6
    threads: 16
    shell: """
        minimap2 -t {threads} -ax map-ont {input} | samtools sort --threads {threads} > {output}
        """

# Generate a depth file from BAM file for MetaBat input
rule metabat_pre:
    input:
        join(outdir, "{samp}/{samp}_lr.bam") if long_read else join(outdir, "{samp}/{samp}.bam") #choose a long read alignment or short read alignment
    output:
        single = join(outdir, "{samp}/{samp}.fa.depth.txt"),
        paired = join(outdir, "{samp}/{samp}.fa.paired.txt"),
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    shell: """
        jgi_summarize_bam_contig_depths --outputDepth {output.single} --pairedContigs {output.paired} --minContigLength 1000 --minContigDepth 1  {input} --percentIdentity 50
        """

#####################################################
################ Binning methods ####################
#####################################################

# Run MetaBat binner
checkpoint metabat:
    input:
        asm = join(outdir, "{samp}/idx/{samp}.fa"),
        depth = join(outdir, "{samp}/{samp}.fa.depth.txt"),
    output:
        directory(join(outdir, "{samp}/metabat/bins/")) #the number of bins is unknown prior to execution
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem=64,
        time=24
    threads: 4
    params:
        outstring = join(outdir, "{samp}/metabat/bins/bin")
    shell: """
        metabat2 --seed 1 -t {threads} --unbinned --inFile {input.asm} --outFile {params.outstring} --abdFile {input.depth}
        """

# Run MaxBin2 binner
checkpoint maxbin:
    input:
        contigs = config['assembly'],
        reads1 = reads[0],
        reads2 = reads[1]
    output:
        directory(join(outdir, "{samp}/maxbin/bins"))
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    params:
        outfolder="{samp}/maxbin/"
    resources:
        time=lambda wildcards, attempt: attempt * 2
    threads: 8
    shell: """
        if [ -d {params.outfolder} ]; then rm -r {params.outfolder}; fi
        mkdir -p {params.outfolder}
        cd {params.outfolder}
        run_MaxBin.pl -contig {input.contigs} -out maxbin \
        -reads {input.reads1} -reads2 {input.reads2} -thread {threads}
        mkdir bins/
        mv *.fasta bins/
        """

# Run myCC binner
checkpoint mycc:
    input:
        contigs = join(outdir, "{samp}/idx/{samp}.fa"),
    output:
        directory(join(outdir, "{samp}/mycc/bins"))
    singularity:
        "docker://990210oliver/mycc.docker:v1"
    params:
        outfolder = join(outdir, "{samp}/mycc"),
        workdir = join(outdir, "{samp}")
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
        config['assembly']
    output:
        cut_contigs = join(outdir, '{samp}/{samp}_contigs_10K.fa'),
        bedfile = join(outdir, '{samp}/{samp}_contigs_10K.bed')
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
        bam = rules.bwa_align.output
    output:
        join(outdir, '{samp}/{samp}_coverage_table.tsv')
    resources:
        mem = 64,
        time = 12
    singularity:
        'docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0'
    shell: """
        concoct_coverage_table.py {input.bedfile} {input.bam} > {output}
    """

rule concoct_run:
    input:
        cut_contigs = rules.concoct_cut.output.cut_contigs,
        coverage = rules.concoct_coverage.output
    output:
        join(outdir, '{samp}/concoct/clustering_gt1000.csv')
    params:
        outdir = join(outdir, '{samp}/concoct')
    threads: 4
    resources:
        mem = 64, 
        time = 24
    singularity:
        'docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0'
    shell: """
        concoct --composition_file {input.cut_contigs} \
        --coverage_file {input.coverage} -b {params.outdir}
    """

rule concoct_merge:
    input:
        clustering = rules.concoct_run.output
    output:
        join(outdir, '{samp}/concoct/clustering_merged.csv')
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
        original_contigs = config['assembly'],
        clustering_merged = rules.concoct_merge.output
    output:
        directory(join(outdir, "{samp}/concoct/bins/")) #the number of bins is unknown prior to execution
    params:
        outdir = join(outdir, '{samp}/concoct/bins/')
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
        lambda wildcards: expand(join(outdir, "{samp}/metabat/bins/{metabat_bin}.fa"), metabat_bin = get_metabat_bins(wildcards), samp = wildcards.samp),
        lambda wildcards: expand(join(outdir, "{samp}/maxbin/bins/{maxbin_bin}.fasta"), maxbin_bin = get_maxbin_bins(wildcards), samp = wildcards.samp),
        lambda wildcards: expand(join(outdir, "{samp}/mycc/bins/{mycc_bin}.fasta"), mycc_bin = get_mycc_bins(wildcards), samp = wildcards.samp),
        lambda wildcards: expand(join(outdir, "{samp}/concoct/bins/{mycc_bin}.fasta"), mycc_bin = get_concoct_bins(wildcards), samp = wildcards.samp),
        contigs = join(outdir, "{samp}/idx/{samp}.fa"),
    output: 
        join(outdir, "{samp}/DAS_tool/completed.txt")
    singularity:
        "shub://ambenj/bin_das_tool:dastool"
    params:
        outfolder = join(outdir, "{samp}/DAS_tool"),
        outfolder_fourmethods = join(outdir, "{samp}/DAS_tool/fourmethods"),
        metabat_dir = join(outdir, "{samp}/metabat/bins/"),
        maxbin_dir = join(outdir, "{samp}/maxbin/bins/"),
        mycc_dir = join(outdir, "{samp}/mycc/bins/"),
        concoct_dir = join(outdir, "{samp}/concoct/bins/"),
        metabat_tsv = join(outdir, "{samp}/DAS_tool/metabat_scaffold2bin.tsv"),
        maxbin_tsv = join(outdir, "{samp}/DAS_tool/maxbin_scaffold2bin.tsv"),
        mycc_tsv = join(outdir, "{samp}/DAS_tool/mycc_scaffold2bin.tsv"),
        concoct_tsv = join(outdir, "{samp}/DAS_tool/concoct_scaffold2bin.tsv"),
    threads: 8
    resources:
        time=lambda wildcards, attempt: attempt * 6
    shell: """
        # Prepare scaffold2bin file for each set of bins
        Fasta_to_Scaffolds2Bin.sh -e fa -i {params.metabat_dir} > {params.metabat_tsv}
        Fasta_to_Scaffolds2Bin.sh -e fasta -i {params.maxbin_dir} > {params.maxbin_tsv}
        Fasta_to_Scaffolds2Bin.sh -e fasta -i {params.mycc_dir} > {params.mycc_tsv}
        Fasta_to_Scaffolds2Bin.sh -e fa -i {params.concoct_dir} > {params.concoct_tsv}
        # Run DAS_Tool
        DAS_Tool -i {params.metabat_tsv},{params.maxbin_tsv},{params.mycc_tsv},{params.concoct_tsv} \
        -l metabat,maxbin,mycc,concoct -c {input.contigs} -o {params.outfolder_fourmethods} \
        --search_engine diamond --threads {threads} --write_bins 1 --write_unbinned 1
        touch {output}
        """

# extract the bins as a separate rule to speed up execution
checkpoint extract_DAStool:
    input: rules.DAStool.output
    output:
        directory(join(outdir, "{samp}/DAS_tool_bins"))
    params:
        old_binfolder = join(outdir, "{samp}/DAS_tool/fourmethods_DASTool_bins"),
        new_binfolder = join(outdir, "{samp}/DAS_tool_bins"),
    shell: """
        mkdir -p {params.new_binfolder}
        cp {params.old_binfolder}/*.fa {params.new_binfolder}
    """


#####################################################
###################### CheckM #######################
#####################################################

# checkm for metabat output
rule checkm_metabat:
    input:
        lambda wildcards: expand(join(outdir, "{samp}/metabat/bins/{metabat_bin}.fa"), metabat_bin = get_metabat_bins(wildcards), samp = wildcards.samp)
    output:
        join(outdir, "{samp}/metabat/checkm/checkm.tsv")
    log:
        join(outdir, "{samp}/metabat/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{samp}/metabat/bins/"),
        checkmfolder = join(outdir, "{samp}/metabat/checkm/"),
    shell: """
        rm -rf {samp}/metabat/checkm/*
        checkm lineage_wf -t {threads} -x fa --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

# checkm for maxbin output
rule checkm_maxbin:
    input:
        lambda wildcards: expand(join(outdir, "{samp}/maxbin/bins/{maxbin_bin}.fasta"), maxbin_bin = get_maxbin_bins(wildcards), samp = wildcards.samp)
    output:
        join(outdir, "{samp}/maxbin/checkm/checkm.tsv")
    log:
        join(outdir, "{samp}/maxbin/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{samp}/maxbin/bins"),
        checkmfolder = join(outdir, "{samp}/maxbin/checkm/")
    shell: """
        rm -rf {samp}/maxbin/checkm/*
        checkm lineage_wf -t {threads} -x fasta --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

rule checkm_mycc:
    input:
        lambda wildcards: expand(join(outdir, "{samp}/mycc/bins/{mycc_bin}.fasta"), mycc_bin = get_mycc_bins(wildcards), samp = wildcards.samp)
    output:
        join(outdir, "{samp}/mycc/checkm/checkm.tsv")
    log:
        join(outdir, "{samp}/mycc/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{samp}/mycc/bins/"),
        checkmfolder = join(outdir, "{samp}/mycc/checkm/"),
    shell: """
        rm -rf {samp}/mycc/checkm/*
        checkm lineage_wf -t {threads} -x fasta --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

# checkm for DAStool output
rule checkm_DAStool:
    input:
        lambda wildcards: expand(join(outdir, "{samp}/DAS_tool_bins/{bin}.fa"), bin = get_DAStool_bins(wildcards), samp = wildcards.samp)
    output:
        join(outdir, "{samp}/DAS_tool/checkm/checkm.tsv")
    log:
        join(outdir, "{samp}/DAS_tool/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{samp}/DAS_tool_bins"),
        checkmfolder = join(outdir, "{samp}/DAS_tool/checkm/"),
        bin_ex = ".fa"
    shell: """
        rm -rf {samp}/DAS_tool/checkm/*
        checkm lineage_wf -t {threads} -x {params.bin_ex} --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

#####################################################
############ ANALYSIS OF DAS_TOOL BINS ##############
#####################################################

# Use aragorn to detect tRNA and tmRNA genes
rule aragorn:
    input:
        join(outdir, "{samp}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{samp}/rna/trna/{bin}.fa.txt")
    log:
        join(outdir, "{samp}/logs/aragorn_{bin}.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem = 8,
        time = 1
    shell: """
        aragorn -t {input} -o {output}
        """

# Use barrnap to predict ribosomal RNA
rule barrnap:
    input:
        join(outdir, "{samp}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{samp}/rna/rrna/{bin}.fa.txt")
    log:
        join(outdir, "{samp}/logs/barrnap_{bin}.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem = 8,
        time = 1
    shell: """
        barrnap {input} > {output}
        """

# Use quast to evaluate genome assemblies
rule quast:
    input:
        join(outdir, "{samp}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{samp}/quast/{bin}.fa/report.tsv")
    log:
        join(outdir, "{samp}/logs/quast_{bin}.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem=8,
        time=1
    params:
        quastfolder = join(outdir, "{samp}/quast/{bin}.fa/"),
        thresholds = "0,10000,50000,100000,250000,500000,1000000,2000000,3000000"
    shell: """
        quast.py -o {params.quastfolder} {input} --contig-thresholds {params.thresholds} --fast
        """

rule pull_prokka:
    input: config['assembly'],
    output: join(outdir, "{samp}/prokka/pulled.txt") 
    singularity:
        "shub://bsiranosian/bens_1337_workflows:prokka"
    shell: """
        touch {output}
    """

# Use prokka to annotate metagenomes
rule prokka:
    input:
        join(outdir, "{samp}/DAS_tool_bins/{bin}.fa"),
        rules.pull_prokka.output
    output:
        join(outdir, "{samp}/prokka/{bin}.fa/{samp}_{bin}.fa.gff")
    log:
        join(outdir, "{samp}/logs/prokka_{bin}.log")
    singularity:
        "shub://bsiranosian/bens_1337_workflows:prokka"
    resources:
        mem = 48,
        time = lambda wildcards, attempt: 4 * attempt,
    threads: 8
    params:
        prokkafolder = join(outdir, "{samp}/prokka/{bin}.fa"),
        prefix = "{samp}_{bin}.fa"
    shell: """
        # don't run this on unbinned contigs, takes forever
        if [ {wildcards.bin} = "bin.unbinned.contigs" ]; then
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
        join(outdir, "{samp}/{samp}_lr.bam") if long_read else join(outdir, "{samp}/{samp}.bam") #choose a long read alignment or short read alignment
    output:
        join(outdir, "{samp}/{samp}_lr.bam.bai") if long_read else join(outdir, "{samp}/{samp}.bam.bai") #choose a long read alignment or short read alignment
    log:
        join(outdir, "{samp}/logs/bamidx.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem = 2,
        time = 2
    shell:
        "samtools index {input}"

# Retrieve stats on mapping from whole metagenome sample
rule bam_idxstats:
    input:
        join(outdir, "{samp}/{samp}_lr.bam") if long_read else join(outdir, "{samp}/{samp}.bam"), #choose a long read alignment or short read alignment,
        join(outdir, "{samp}/{samp}_lr.bam.bai") if long_read else join(outdir, "{samp}/{samp}.bam.bai"), #choose a long read alignment or short read alignment,
    output:
        join(outdir, "{samp}/{samp}_lr.bam.bai.tsv") if long_read else join(outdir, "{samp}/{samp}.bam.bai.tsv"), #choose a long read alignment or short read alignment,
    log:
        join(outdir, "{samp}/logs/bamidxstats.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem = 2,
        time = 2
    shell:
        "samtools idxstats {input[0]} > {output}"

# Get mapping stats for each bin
rule bin_idxstats:
    input:
        join(outdir, "{samp}/DAS_tool_bins/{bin}.fa"),
        join(outdir, "{samp}/{samp}_lr.bam.bai.tsv") if long_read else join(outdir, "{samp}/{samp}.bam.bai.tsv"), #choose a long read alignment or short read alignment,
    output:
        join(outdir, "{samp}/coverage/raw/{bin}.tsv")
    log:
        join(outdir, "{samp}/logs/coverage_idxstats_{bin}.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem = 2,
        time = 6
    shell:
        "grep '>' {input[0]} | tr -d '>' | xargs -I foo -n 1 grep -P 'foo\t' {input[1]} > {output}"

# Determine coverage for each bin
rule bin_coverage:
    input:
        rules.bin_idxstats.output
    output:
        join(outdir, "{samp}/coverage/{bin}.txt")
    log:
        join(outdir, "{samp}/logs/coverage_{bin}.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem = 2,
        time = 1
    params:
        read_length = config['read_length']
    script:
        "scripts/bin_coverage.py"

# index DAStool bins
rule fasta_index:
    input:
        join(outdir, "{samp}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{samp}/DAS_tool_bins/{bin}.fa.fai")
    log:
        join(outdir, "{samp}/logs/faidx_{bin}.log")
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
        join(outdir, "{samp}/idx/{samp}.fa")
    output:
        krak = join(outdir, "{samp}/classify/{samp}.krak"),
        krak_report = join(outdir, "{samp}/classify/{samp}.krak.report")
    log:
        join(outdir, "{samp}/logs/kraken_class.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    params:
        db = config['kraken2db']
    resources:
        mem = 256,
        time = 6
    threads: 4
    shell: """
        kraken2 --db {params.db} --db {params.db} --threads {threads} \
        --output {output.krak} --report {output.krak_report} {input}
        """

rule label_bins:
    input:
        krak = join(outdir, "{samp}/classify/{samp}.krak"),
        bins = lambda wildcards: expand(join(outdir, "{samp}/DAS_tool_bins/{bin}.fa.fai"),
                                        bin = get_DAStool_bins(wildcards), samp = wildcards.samp)
    output:
        join(outdir, "{samp}/classify/bin_species_calls.tsv")
    log:
        join(outdir, "{samp}/logs/assign_species.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    params:
        binfolder = join(outdir, "{samp}/DAS_tool_bins/"),
        custom_taxonomy = custom_taxonomy
    script:
        "scripts/assign_species.py"

rule postprocess:
    input:
        prokka = lambda wildcards: expand(join(outdir, "{samp}/prokka/{bin}.fa/{samp}_{bin}.fa.gff"), bin = get_bins(wildcards), samp = wildcards.samp),
        quast = lambda wildcards: expand(join(outdir, "{samp}/quast/{bin}.fa/report.tsv"), bin = get_bins(wildcards), samp = wildcards.samp),
        checkm = join(outdir, "{samp}/checkm/checkm.tsv"),
        trna = lambda wildcards: expand(join(outdir, "{samp}/rna/trna/{bin}.fa.txt"), bin = get_bins(wildcards), samp = wildcards.samp),
        rrna = lambda wildcards: expand(join(outdir, "{samp}/rna/rrna/{bin}.fa.txt"), bin = get_bins(wildcards), samp = wildcards.samp),
        classify = rules.label_bins.output,
        coverage = lambda wildcards: expand(join(outdir, "{samp}/coverage/{bin}.txt"), bin = get_bins(wildcards), samp = wildcards.samp),
    output:
        full = join(outdir, "{samp}/final/{samp}.tsv"),
        simple = join(outdir, "{samp}/final/{samp}_simple.tsv")
    params:
        bins = lambda wildcards: get_DAStool_bins(wildcards),
        sample = samp
    script: "scripts/postprocess.R"
