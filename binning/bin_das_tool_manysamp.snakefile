from os.path import join, abspath, expanduser
localrules: bwa_index_setup, postprocess, label_bins, extract_DAStool, concoct_extract_bins

def get_sample_assemblies_reads(sample_file):
    sample_reads = {}
    sample_assemblies = {}
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            if len(s) < 3:
                sys.exit('Badly formatted sample_file')
            sample = s[0]
            assembly = s[1]
            reads_split = s[2].split(',')
            if sample in sample_reads:
                print(sample)
                raise ValueError("Non-unique sample encountered!")
            # get read pairs and singles from read specification
            if (len(reads_split) == 3) or (len(reads_split) == 2):
                sample_reads[sample] = reads_split[0:2]
            else:
                sys.exit('must be paired end reads')
            sample_assemblies[sample] = assembly

    return sample_reads, sample_assemblies

# Read in sample and outdir from config file
sample_file = config['sample_file']
outdir = config['outdir_base']
sample_reads, sample_assemblies = get_sample_assemblies_reads(sample_file)
sample_list = list(sample_reads.keys())

print(sample_list)
print(sample_reads)
print(sample_assemblies)

# are we using a non-standard (non ncbi) taxonomy
if 'custom_taxonomy' in config:
    custom_taxonomy = config['custom_taxonomy']
else:
    custom_taxonomy = False

# convert outdir to absolute path
if outdir[0] == '~':
    outdir = expanduser(outdir)
outdir = abspath(outdir)


# Determine if long reads
if 'long_read' in config and config['long_read']:
    long_read = True
    sys.exit('Untested in this implementation')
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
        config['kraken2db'],
        # checkm - not necessary to run for all binners
        # expand(join(outdir, "{sample}/metabat/checkm/checkm.tsv"), sample = config['sample']),
        # expand(join(outdir, "{sample}/maxbin/checkm/checkm.tsv"), sample = config['sample']),
        # expand(join(outdir, "{sample}/mycc/checkm/checkm.tsv"), sample = config['sample']),
        # Post-processing
        expand(join(outdir, "{sample}/classify/bin_species_calls.tsv"), sample = sample_list),
        expand(join(outdir, "{sample}/final/{sample}.tsv"), sample = sample_list),
        expand(join(outdir, "{sample}/final/{sample}_simple.tsv"), sample = sample_list),

##########################################
####### Prepping input for binning #######
##########################################

# copy assembly fasta file to working directory
rule bwa_index_setup:
    input:
        lambda wildcards: sample_assemblies[wildcards.sample]
    output:
        join(outdir, "{sample}/idx/{sample}.fa")
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
        join(outdir, "{sample}/idx/{sample}.fa")
    output:
        join(outdir, "{sample}/idx/{sample}.fa.amb"),
        join(outdir, "{sample}/idx/{sample}.fa.ann"),
        join(outdir, "{sample}/idx/{sample}.fa.bwt"),
        join(outdir, "{sample}/idx/{sample}.fa.pac"),
        join(outdir, "{sample}/idx/{sample}.fa.sa")
    log:
        join(outdir, "{sample}/logs/bwa_index.log")
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
        asm = join(outdir, "{sample}/idx/{sample}.fa"),
        amb = join(outdir, "{sample}/idx/{sample}.fa.amb"),
        ann = join(outdir, "{sample}/idx/{sample}.fa.ann"),
        bwt = join(outdir, "{sample}/idx/{sample}.fa.bwt"),
        pac = join(outdir, "{sample}/idx/{sample}.fa.pac"),
        sa = join(outdir, "{sample}/idx/{sample}.fa.sa")
    output:
        join(outdir, "{sample}/{sample}.bam")
    log:
        join(outdir, "{sample}/logs/bwa_mem.log")
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
        join(outdir, "{sample}/{sample}_lr.bam") if long_read else join(outdir, "{sample}/{sample}.bam") #choose a long read alignment or short read alignment
    output:
        single = join(outdir, "{sample}/{sample}.fa.depth.txt"),
        paired = join(outdir, "{sample}/{sample}.fa.paired.txt"),
    singularity:
        "docker://quay.io/biocontainers/metabat2:2.15--h137b6e9_0"
    shell: """
        jgi_summarize_bam_contig_depths --outputDepth {output.single} --pairedContigs {output.paired} --minContigLength 1000 --minContigDepth 1  {input} --percentIdentity 50
        """

#####################################################
################ Binning methods ####################
#####################################################

# Run MetaBat binner
checkpoint metabat:
    input:
        asm = join(outdir, "{sample}/idx/{sample}.fa"),
        depth = join(outdir, "{sample}/{sample}.fa.depth.txt"),
    output:
        directory(join(outdir, "{sample}/metabat/bins/")) #the number of bins is unknown prior to execution
    singularity:
        "docker://quay.io/biocontainers/metabat2:2.15--h137b6e9_0"
    resources:
        mem=64,
        time=24
    threads: 4
    params:
        outstring = join(outdir, "{sample}/metabat/bins/bin")
    shell: """
        metabat2 --seed 1 -t {threads} --unbinned --inFile {input.asm} --outFile {params.outstring} --abdFile {input.depth}
        """

# Run MaxBin2 binner
checkpoint maxbin:
    input:
        contigs = join(outdir, "{sample}/idx/{sample}.fa"),
        fwd = lambda wildcards: sample_reads[wildcards.sample][0],
        rev = lambda wildcards: sample_reads[wildcards.sample][1],
        depth = rules.metabat_pre.output.single
    output:
        directory(join(outdir, "{sample}/maxbin/bins"))
    singularity:
        "docker://quay.io/biocontainers/maxbin2:2.2.7--he1b5a44_1"
    params:
        outfolder = join(outdir, "{sample}/maxbin/"),
        abunance_file = join(outdir, "{sample}/maxbin/{sample}_maxbin_depth.txt")
    resources:
        time=lambda wildcards, attempt: attempt * 2
    threads: 8
    shell: """
        if [ -d {params.outfolder} ]; then rm -r {params.outfolder}; fi
        mkdir -p {params.outfolder}
        cd {params.outfolder}

        # create abundance file, which we already have from the metabat rule
        cut -f 1,4 {input.depth} | tail -n +2 > {params.abundance_file}
        run_MaxBin.pl -contig {input.contigs} -out maxbin \
        -abund {params.abundance_file} -thread {threads}
        mkdir bins/
        mv *.fasta bins/
        """

# Run myCC binner
checkpoint mycc:
    input:
        contigs = join(outdir, "{sample}/idx/{sample}.fa"),
    output:
        directory(join(outdir, "{sample}/mycc/bins"))
    singularity:
        "docker://990210oliver/mycc.docker:v1"
    params:
        outfolder = join(outdir, "{sample}/mycc"),
        outfolder_bins = join(outdir, "{sample}/mycc/bins"),
        workdir = join(outdir, "{sample}")
    resources:
        time=lambda wildcards, attempt: attempt * 12,
        mem=lambda wildcards, attempt: attempt * 128
    shell: """
        if [ -d {params.outfolder} ]; then rm -r {params.outfolder}; fi
        cd {params.workdir}
        MyCC.py {input.contigs} -meta
        mv 20*/ {params.outfolder}/  # change variable folder name to mycc
        mkdir {params.outfolder_bins}
        mv {params.outfolder}/*.fasta {params.outfolder_bins}
        """

# Run concoct binner in multiple steps
rule concoct_cut:
    input:
        join(outdir, "{sample}/idx/{sample}.fa"),
    output:
        cut_contigs = join(outdir, '{sample}/{sample}_contigs_10K.fa'),
        bedfile = join(outdir, '{sample}/{sample}_contigs_10K.bed')
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
        join(outdir, '{sample}/{sample}_coverage_table.tsv')
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
        join(outdir, '{sample}/concoct/clustering_gt1000.csv')
    params:
        outdir = join(outdir, '{sample}/concoct')
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
        join(outdir, '{sample}/concoct/clustering_merged.csv')
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
        original_contigs = join(outdir, "{sample}/idx/{sample}.fa"),
        clustering_merged = rules.concoct_merge.output
    output:
        directory(join(outdir, "{sample}/concoct/bins/")) #the number of bins is unknown prior to execution
    params:
        outdir = join(outdir, '{sample}/concoct/bins/')
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
        lambda wildcards: expand(join(outdir, "{sample}/metabat/bins/{metabat_bin}.fa"), metabat_bin = get_metabat_bins(wildcards), sample = wildcards.sample),
        lambda wildcards: expand(join(outdir, "{sample}/maxbin/bins/{maxbin_bin}.fasta"), maxbin_bin = get_maxbin_bins(wildcards), sample = wildcards.sample),
        lambda wildcards: expand(join(outdir, "{sample}/mycc/bins/{mycc_bin}.fasta"), mycc_bin = get_mycc_bins(wildcards), sample = wildcards.sample),
        lambda wildcards: expand(join(outdir, "{sample}/concoct/bins/{mycc_bin}.fasta"), mycc_bin = get_concoct_bins(wildcards), sample = wildcards.sample),
        contigs = join(outdir, "{sample}/idx/{sample}.fa"),
    output: 
        join(outdir, "{sample}/DAS_tool/completed.txt")
    # singularity:
        # "docker://quay.io/biocontainers/das_tool:1.1.1--py36r351_1"
    conda: "envs/das_tool.yaml"
    params:
        outfolder = join(outdir, "{sample}/DAS_tool"),
        outfolder_fourmethods = join(outdir, "{sample}/DAS_tool/fourmethods"),
        metabat_dir = join(outdir, "{sample}/metabat/bins/"),
        maxbin_dir = join(outdir, "{sample}/maxbin/bins/"),
        mycc_dir = join(outdir, "{sample}/mycc/bins/"),
        concoct_dir = join(outdir, "{sample}/concoct/bins/"),
        metabat_tsv = join(outdir, "{sample}/DAS_tool/metabat_scaffold2bin.tsv"),
        maxbin_tsv = join(outdir, "{sample}/DAS_tool/maxbin_scaffold2bin.tsv"),
        mycc_tsv = join(outdir, "{sample}/DAS_tool/mycc_scaffold2bin.tsv"),
        concoct_tsv = join(outdir, "{sample}/DAS_tool/concoct_scaffold2bin.tsv"),
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
        directory(join(outdir, "{sample}/DAS_tool_bins"))
    params:
        old_binfolder = join(outdir, "{sample}/DAS_tool/fourmethods_DASTool_bins"),
        new_binfolder = join(outdir, "{sample}/DAS_tool_bins"),
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
        lambda wildcards: expand(join(outdir, "{sample}/metabat/bins/{metabat_bin}.fa"), metabat_bin = get_metabat_bins(wildcards), sample = wildcards.sample)
    output:
        join(outdir, "{sample}/metabat/checkm/checkm.tsv")
    log:
        join(outdir, "{sample}/metabat/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{sample}/metabat/bins/"),
        checkmfolder = join(outdir, "{sample}/metabat/checkm/"),
    shell: """
        rm -rf {sample}/metabat/checkm/*
        checkm lineage_wf -t {threads} -x fa --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

# checkm for maxbin output
rule checkm_maxbin:
    input:
        lambda wildcards: expand(join(outdir, "{sample}/maxbin/bins/{maxbin_bin}.fasta"), maxbin_bin = get_maxbin_bins(wildcards), sample = wildcards.sample)
    output:
        join(outdir, "{sample}/maxbin/checkm/checkm.tsv")
    log:
        join(outdir, "{sample}/maxbin/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{sample}/maxbin/bins"),
        checkmfolder = join(outdir, "{sample}/maxbin/checkm/")
    shell: """
        rm -rf {sample}/maxbin/checkm/*
        checkm lineage_wf -t {threads} -x fasta --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

rule checkm_mycc:
    input:
        lambda wildcards: expand(join(outdir, "{sample}/mycc/bins/{mycc_bin}.fasta"), mycc_bin = get_mycc_bins(wildcards), sample = wildcards.sample)
    output:
        join(outdir, "{sample}/mycc/checkm/checkm.tsv")
    log:
        join(outdir, "{sample}/mycc/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{sample}/mycc/bins/"),
        checkmfolder = join(outdir, "{sample}/mycc/checkm/"),
    shell: """
        rm -rf {sample}/mycc/checkm/*
        checkm lineage_wf -t {threads} -x fasta --tab_table -f {output} {params.binfolder} {params.checkmfolder}
        """

# checkm for DAStool output
rule checkm_DAStool:
    input:
        lambda wildcards: expand(join(outdir, "{sample}/DAS_tool_bins/{bin}.fa"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample)
    output:
        join(outdir, "{sample}/DAS_tool/checkm/checkm.tsv")
    log:
        join(outdir, "{sample}/DAS_tool/logs/checkm.log")
    singularity:
        "shub://bsiranosian/bin_genomes:checkm"
    resources:
        mem = 128,
        time = 12
    threads: 4
    params:
        binfolder = join(outdir, "{sample}/DAS_tool_bins"),
        checkmfolder = join(outdir, "{sample}/DAS_tool/checkm/"),
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
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{sample}/rna/trna/{bin}.fa.txt")
    log:
        join(outdir, "{sample}/logs/aragorn_{bin}.log")
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
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{sample}/rna/rrna/{bin}.fa.txt")
    log:
        join(outdir, "{sample}/logs/barrnap_{bin}.log")
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
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{sample}/quast/{bin}.fa/report.tsv")
    log:
        join(outdir, "{sample}/logs/quast_{bin}.log")
    singularity:
        "docker://quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0"
    resources:
        mem=8,
        time=1
    params:
        quastfolder = join(outdir, "{sample}/quast/{bin}.fa/"),
        thresholds = "0,10000,50000,100000,250000,500000,1000000,2000000,3000000"
    shell: """
        quast.py -o {params.quastfolder} {input} --contig-thresholds {params.thresholds} --fast
        """

rule pull_prokka:
    input: join(outdir, "{sample}/idx/{sample}.fa"),
    output: join(outdir, "{sample}/prokka/pulled.txt") 
    singularity:
        "docker://quay.io/biocontainers/prokka:1.14.5--pl526_0"
    shell: """
        touch {output}
    """

# Use prokka to annotate metagenomes
rule prokka:
    input:
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa"),
        rules.pull_prokka.output
    output:
        join(outdir, "{sample}/prokka/{bin}.fa/{sample}_{bin}.fa.gff")
    log:
        join(outdir, "{sample}/logs/prokka_{bin}.log")
    singularity:
        "docker://quay.io/biocontainers/prokka:1.14.5--pl526_0"
    resources:
        mem = 48,
        time = lambda wildcards, attempt: 4 * attempt,
    threads: 8
    params:
        prokkafolder = join(outdir, "{sample}/prokka/{bin}.fa"),
        prefix = "{sample}_{bin}.fa"
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
        join(outdir, "{sample}/{sample}_lr.bam") if long_read else join(outdir, "{sample}/{sample}.bam") #choose a long read alignment or short read alignment
    output:
        join(outdir, "{sample}/{sample}_lr.bam.bai") if long_read else join(outdir, "{sample}/{sample}.bam.bai") #choose a long read alignment or short read alignment
    log:
        join(outdir, "{sample}/logs/bamidx.log")
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
        join(outdir, "{sample}/{sample}_lr.bam") if long_read else join(outdir, "{sample}/{sample}.bam"), #choose a long read alignment or short read alignment,
        join(outdir, "{sample}/{sample}_lr.bam.bai") if long_read else join(outdir, "{sample}/{sample}.bam.bai"), #choose a long read alignment or short read alignment,
    output:
        join(outdir, "{sample}/{sample}_lr.bam.bai.tsv") if long_read else join(outdir, "{sample}/{sample}.bam.bai.tsv"), #choose a long read alignment or short read alignment,
    log:
        join(outdir, "{sample}/logs/bamidxstats.log")
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
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa"),
        join(outdir, "{sample}/{sample}_lr.bam.bai.tsv") if long_read else join(outdir, "{sample}/{sample}.bam.bai.tsv"), #choose a long read alignment or short read alignment,
    output:
        join(outdir, "{sample}/coverage/raw/{bin}.tsv")
    log:
        join(outdir, "{sample}/logs/coverage_idxstats_{bin}.log")
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
        join(outdir, "{sample}/coverage/{bin}.txt")
    log:
        join(outdir, "{sample}/logs/coverage_{bin}.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    resources:
        mem = 2,
        time = 1
    params:
        read_length = config['read_length'],
        sample = lambda wildcards: wildcards.sample
    script:
        "scripts/bin_coverage.py"

# index DAStool bins
rule fasta_index:
    input:
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa")
    output:
        join(outdir, "{sample}/DAS_tool_bins/{bin}.fa.fai")
    log:
        join(outdir, "{sample}/logs/faidx_{bin}.log")
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
        join(outdir, "{sample}/idx/{sample}.fa")
    output:
        krak = join(outdir, "{sample}/classify/{sample}.krak"),
        krak_report = join(outdir, "{sample}/classify/{sample}.krak.report")
    log:
        join(outdir, "{sample}/logs/kraken_class.log")
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
        krak = join(outdir, "{sample}/classify/{sample}.krak"),
        bins = lambda wildcards: expand(join(outdir, "{sample}/DAS_tool_bins/{bin}.fa.fai"),
                                        bin = get_DAStool_bins(wildcards), sample = wildcards.sample)
    output:
        join(outdir, "{sample}/classify/bin_species_calls.tsv")
    log:
        join(outdir, "{sample}/logs/assign_species.log")
    singularity:
        "shub://bsiranosian/bin_genomes:binning"
    params:
        binfolder = join(outdir, "{sample}/DAS_tool_bins/"),
        custom_taxonomy = custom_taxonomy
    script:
        "scripts/assign_species.py"

rule postprocess:
    input:
        prokka = lambda wildcards: expand(join(outdir, "{sample}/prokka/{bin}.fa/{sample}_{bin}.fa.gff"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample),
        quast = lambda wildcards: expand(join(outdir, "{sample}/quast/{bin}.fa/report.tsv"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample),
        checkm = join(outdir, "{sample}/DAS_tool/checkm/checkm.tsv"),
        trna = lambda wildcards: expand(join(outdir, "{sample}/rna/trna/{bin}.fa.txt"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample),
        rrna = lambda wildcards: expand(join(outdir, "{sample}/rna/rrna/{bin}.fa.txt"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample),
        classify = rules.label_bins.output,
        coverage = lambda wildcards: expand(join(outdir, "{sample}/coverage/{bin}.txt"), bin = get_DAStool_bins(wildcards), sample = wildcards.sample),
    output:
        full = join(outdir, "{sample}/final/{sample}.tsv"),
        simple = join(outdir, "{sample}/final/{sample}_simple.tsv")
    params:
        bins = lambda wildcards: get_DAStool_bins(wildcards),
        sample = lambda wildcards: wildcards.sample
    script: "scripts/postprocess.R"
