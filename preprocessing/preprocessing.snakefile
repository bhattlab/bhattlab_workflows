import re,os,subprocess
from os.path import join, expanduser, abspath

################################################################################
# specify project directories
DATA_DIR    = config["raw_reads_directory"]
PROJECT_DIR = config["output_directory"]
READ_SUFFIX = config["read_specification"]
EXTENSION   = config["extension"]
# if gzipped, set this. otherwise not
gz_ext = '.gz' if EXTENSION.endswith('.gz') else ''
# print(gz_ext)

# convert PROJECT_DIR and DATA_DIR to absolute path
if PROJECT_DIR[0] == '~':
    PROJECT_DIR = expanduser(PROJECT_DIR)
PROJECT_DIR = abspath(PROJECT_DIR)
if DATA_DIR[0] == '~':
    DATA_DIR = expanduser(DATA_DIR)
DATA_DIR = abspath(DATA_DIR)

# get file names
FILES = [f for f in os.listdir(DATA_DIR) if f.endswith(EXTENSION)]
SAMPLE_PREFIX = list(set([re.split('|'.join(['_' + a + '\.' for a in READ_SUFFIX]), i)[0] for i in FILES]))

# config for trim_galore: only enable start_trim if >0
start_trim = config['trim_galore']['start_trim']
end_trim = config['trim_galore']['end_trim']
if (start_trim) >0: 
    start_trim_string = '--clip_R1 {a} --clip_R2 {a}'.format(a=str(start_trim))
else: 
    start_trim_string = ''
if (end_trim) >0: 
    end_trim_string = '--three_prime_clip_R1 {a} --three_prime_clip_R2 {a}'.format(a=str(end_trim))
else: 
    end_trim_string = ''

################################################################################
localrules: assembly_meta_file, pre_multiqc, post_multiqc, cleanup
rule all:
    input:
        expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX),
        expand(join(PROJECT_DIR, "01_processing/00_qc_reports/pre_fastqc/{sample}_{read}_fastqc.html"), sample=SAMPLE_PREFIX, read=READ_SUFFIX),
        expand(join(PROJECT_DIR, "01_processing/00_qc_reports/post_fastqc/{sample}_{read}_fastqc.html"), sample=SAMPLE_PREFIX, read=['1', '2', 'orphans']),
        join(PROJECT_DIR, "01_processing/00_qc_reports/pre_multiqc/multiqc_report.html"),
        join(PROJECT_DIR, "01_processing/00_qc_reports/post_multiqc/multiqc_report.html"),
        join(PROJECT_DIR, "01_processing/assembly_input.txt"),
        join(PROJECT_DIR, "01_processing/classification_input.txt"),
        join(PROJECT_DIR, "01_processing/readcounts.tsv"),
        join(PROJECT_DIR, "01_processing/readcounts.pdf")

################################################################################
rule pre_fastqc:
    input:  
        fwd = join(DATA_DIR, "{sample}_" + READ_SUFFIX[0] + EXTENSION),
        rev = join(DATA_DIR, "{sample}_" + READ_SUFFIX[1] + EXTENSION)
    output:
        fwd = join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{sample}_" + READ_SUFFIX[0] + "_fastqc.html"),
        rev = join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{sample}_" + READ_SUFFIX[1] + "_fastqc.html")
    params:
        outdir = join(PROJECT_DIR, "01_processing/00_qc_reports/pre_fastqc/")
    threads: min(2, len(READ_SUFFIX))
    resources:
            time = 6,
            mem = 32
    singularity: "docker://quay.io/biocontainers/fastqc:0.11.8--1"
    benchmark: join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{sample}_time.txt")
    shell: """
        mkdir -p {params.outdir}
        fastqc {input} --outdir {params.outdir} --threads {threads}
    """

rule pre_multiqc:
    input:
        expand(join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{sample}_{read}_fastqc.html"), sample=SAMPLE_PREFIX, read=READ_SUFFIX)
    output: join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_multiqc/multiqc_report.html")
    params:
        indir = join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc"),
        outdir = join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_multiqc/")
    singularity: "docker://quay.io/biocontainers/multiqc:1.7--py_2"
    shell: """
        multiqc --force {params.indir} -o {params.outdir}
    """

###############################################################################
rule deduplicate:
    input:
        fwd = join(DATA_DIR, "{sample}_" + READ_SUFFIX[0] + EXTENSION),
        rev = join(DATA_DIR, "{sample}_" + READ_SUFFIX[1] + EXTENSION)
    output:
        fwd = join(PROJECT_DIR, "01_processing/01_dedup/{sample}_1.fq.gz"),
        rev = join(PROJECT_DIR, "01_processing/01_dedup/{sample}_2.fq.gz"),
    params:
        tmp_fwd = '{sample}_R1.fastq.gz',
        tmp_rev = '{sample}_R2.fastq.gz',
        outdir = join(PROJECT_DIR, "01_processing/01_dedup/")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        mem = lambda wildcards, attempt: attempt * 16, 
        time = 24
    singularity: "docker://dzs74/htstream"
    benchmark: join(PROJECT_DIR, "01_processing/01_dedup/{sample}_time.txt")
    shell: """
        mkdir -p {params.outdir} && cd {params.outdir}
        hts_SuperDeduper -1 {input.fwd} -2 {input.rev} -f {wildcards.sample} -F
        mv {params.tmp_fwd} {output.fwd}
        mv {params.tmp_rev} {output.rev}
    """
###############################################################################

################################################################################
rule trim_galore:
    input:
        fwd = rules.deduplicate.output.fwd,
        rev = rules.deduplicate.output.rev,
    output:
        fwd = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_1_val_1.fq" + gz_ext),
        rev = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_2_val_2.fq" + gz_ext),
        orp = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_unpaired.fq" + gz_ext)
    threads: 2
    resources:
        mem=32,
        time=lambda wildcards, attempt: attempt * 24
    params:
        orp_fwd = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_1_unpaired_1.fq" + gz_ext),
        orp_rev = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_2_unpaired_2.fq" + gz_ext),
        # output_fwd_pre_gz = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_") + READ_SUFFIX[0] + "_val_1.fq",
        # output_rev_pre_gz = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_") + READ_SUFFIX[1] + "_val_2.fq",
        q_min   = config['trim_galore']['quality'],
        left    = config['trim_galore']['start_trim'],
        min_len = config['trim_galore']['min_read_length'],
        outdir  = join(PROJECT_DIR, "01_processing/02_trimmed/"),
        gz_output = str(gz_ext == '.gz').lower()
    singularity: "docker://quay.io/biocontainers/trim-galore:0.6.5--0"
    benchmark: join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_time.txt")
    shell: """
        mkdir -p {params.outdir}
        trim_galore --quality {params.q_min} \
            --length {params.min_len} \
            --output_dir {params.outdir} \
            --paired {input.fwd} {input.rev} \
            --retain_unpaired \
            --cores {threads} \
            {start_trim_string} \
            {end_trim_string}

        # if output is gz, merge unpaired and gzip
        if {params.gz_output}; then
            zcat -f {params.orp_fwd} {params.orp_rev} | pigz -b 32 -p {threads} > {output.orp}
        else
            zcat -f {params.orp_fwd} {params.orp_rev} > {output.orp}
        fi
        # delete intermediate files
        rm {params.orp_fwd} {params.orp_rev}
    """

################################################################################
rule rm_host_reads:
    input:
        index_amb = config['rm_host_reads']['host_genome'] + '.amb',
        index_ann = config['rm_host_reads']['host_genome'] + '.ann',
        index_bwt = config['rm_host_reads']['host_genome'] + '.bwt',
        index_pac = config['rm_host_reads']['host_genome'] + '.pac',
        index_sa = config['rm_host_reads']['host_genome'] + '.sa',
        fwd       = rules.trim_galore.output.fwd,
        rev       = rules.trim_galore.output.rev,
        orp       = rules.trim_galore.output.orp
    output:
        unmapped_1 = join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq.gz"),
        unmapped_2 = join(PROJECT_DIR, "01_processing/05_sync/{sample}_2.fq.gz"),
        unmapped_singletons = join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"),
    params:
        bwa_index_base = join(config['rm_host_reads']['host_genome']),
        singelton_temp_1 = join(PROJECT_DIR, "01_processing/04_host_align/{sample}_rmHost_singletons1.fq.gz"),
        singelton_temp_2 = join(PROJECT_DIR, "01_processing/04_host_align/{sample}_rmHost_singletons2.fq.gz"),
    threads: 4
    resources:
        mem_mb=32000,
        mem=32,
        time=24
    singularity: "shub://bsiranosian/bens_1337_workflows:align"
    # conda: "envs/align.yaml"
    benchmark: join(PROJECT_DIR, "01_processing/04_host_align/{sample}_time.txt")
    shell: """
        mkdir -p {PROJECT_DIR}/01_processing/04_host_align/
        # if an index needs to be built, use bwa index ref.fa
        # run on paired reads
        bwa mem -t {threads} {params.bwa_index_base} {input.fwd} {input.rev} | \
            samtools fastq -@ {threads} -t -T BX -f 4 -1 {output.unmapped_1} -2 {output.unmapped_2} -s {params.singelton_temp_1} - 
        # run on unpaired reads
        bwa mem -t {threads} {params.bwa_index_base} {input.orp} | \
            samtools fastq -@ {threads} -t -T BX -f 4 - > {params.singelton_temp_2}
        # combine singletons 
        zcat -f {params.singelton_temp_1} {params.singelton_temp_2} | pigz -p {threads} > {output.unmapped_singletons}
        rm {params.singelton_temp_1} {params.singelton_temp_2}
    """

################################################################################
rule post_fastqc:
    input:  join(PROJECT_DIR, "01_processing/05_sync/{sample}_{read}.fq.gz")
    output: join(PROJECT_DIR,  "01_processing/00_qc_reports/post_fastqc/{sample}_{read}_fastqc.html"),
    params:
        outdir = join(PROJECT_DIR, "01_processing/00_qc_reports/post_fastqc/")
    threads: 1
    resources:
            time = 6,
            mem = 32
    singularity: "docker://quay.io/biocontainers/fastqc:0.11.8--1"
    benchmark: join(PROJECT_DIR, "01_processing/00_qc_reports/post_fastqc/{sample}_{read}_time.txt")
    shell: """
        mkdir -p {params.outdir}
        if [ -z $(gzip -cd {input} | head -c1) ]; then
            echo EMPTY!
            touch {output}
        else
            fastqc {input} -f fastq --outdir {params.outdir}
        fi
    """

rule post_multiqc:
    input: expand(join(PROJECT_DIR,  "01_processing/00_qc_reports/post_fastqc/{sample}_{read}_fastqc.html"), sample=SAMPLE_PREFIX, read=['1', '2', 'orphans']),
    output: join(PROJECT_DIR,  "01_processing/00_qc_reports/post_multiqc/multiqc_report.html")
    params:
        indir = join(PROJECT_DIR,  "01_processing/00_qc_reports/post_fastqc"),
        outdir = join(PROJECT_DIR,  "01_processing/00_qc_reports/post_multiqc/")
    singularity: "docker://quay.io/biocontainers/multiqc:1.7--py_2"
    shell: """
        multiqc --force {params.indir} -o {params.outdir}
    """

################################################################################
rule assembly_meta_file:
    input: expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX)
    output: join(PROJECT_DIR, "01_processing/assembly_input.txt")
    run:
        outfile = str(output)
        if (os.path.exists(outfile)):
            os.remove(outfile)
        with open(outfile, 'w') as outf:
            outf.writelines(['# Sample\tReads1.fq[.gz][,Reads2.fq[.gz][,orphans.fq[.gz]]]\n'])
            for sample in SAMPLE_PREFIX:
                outline = [sample, ','.join([
                join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq.gz"),
                join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_2.fq.gz"),
                join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_orphans.fq.gz")])]
                outf.writelines('\t'.join(outline) + '\n')

################################################################################
rule classification_meta_file:
    input: expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX)
    output: join(PROJECT_DIR, "01_processing/classification_input.txt")
    run:
        outfile = str(output)
        if (os.path.exists(outfile)):
            os.remove(outfile)
        with open(outfile, 'w') as outf:
            outf.writelines(['# Sample\tr1\tr2\n'])
            for sample in SAMPLE_PREFIX:
                outline = [sample, '\t'.join([
                join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq.gz"),
                join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_2.fq.gz")])]
                outf.writelines('\t'.join(outline) + '\n')

################################################################################
def file_len(fname):
    p = subprocess.Popen('zcat -f ' + fname + ' | wc -l', stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE, shell=True)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

rule readcounts:
    input:
        raw = expand(join(DATA_DIR, "{sample}_" + READ_SUFFIX[0] + EXTENSION), sample=SAMPLE_PREFIX),
        dedup = expand(join(PROJECT_DIR, "01_processing/01_dedup/{sample}_1.fq.gz"), sample=SAMPLE_PREFIX),
        trimmed = expand(join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_1_val_1.fq" + gz_ext), sample=SAMPLE_PREFIX),
        rmhost = expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq.gz"), sample=SAMPLE_PREFIX),
        orphans = expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX)
    output:
        join(PROJECT_DIR, "01_processing/readcounts.tsv")
    resources:
        time = 24
    run:
        outfile = str(output)
        if (os.path.exists(outfile)):
            os.remove(outfile)
        with open(outfile, 'w') as outf:
            outf.writelines('Sample\traw_reads\tdedup_reads\tdedup_frac\ttrimmed_reads\ttrimmed_frac\thost_removed_reads\thost_removed_frac\torphan_reads\torphan_frac\n')
            for sample in SAMPLE_PREFIX:
                raw_file = join(DATA_DIR, sample + "_" + READ_SUFFIX[0] + EXTENSION)
                dedup_file = join(PROJECT_DIR, "01_processing/01_dedup/" + sample + "_1.fq.gz")
                trimmed_file = join(PROJECT_DIR, "01_processing/02_trimmed/" + sample + "_1_val_1.fq" + gz_ext)
                rmhost_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq")
                orphans_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_orphans.fq")

                raw_reads = int(file_len(raw_file) / 4)
                dedup_reads = int(file_len(dedup_file) / 4)
                trimmed_reads = int(file_len(trimmed_file) / 4)
                rmhost_reads = int(file_len(rmhost_file) / 4)
                orphans_reads = int(file_len(orphans_file) / 4)

                dedup_frac = round(dedup_reads / float(raw_reads), 3)
                trimmed_frac = round(trimmed_reads / float(raw_reads), 3)
                rmhost_frac = round(rmhost_reads / float(raw_reads), 3)
                orphans_frac = round(orphans_reads / float(raw_reads), 3)

                line = '\t'.join([sample, str(raw_reads),
                    str(dedup_reads), str(dedup_frac),
                    str(trimmed_reads), str(trimmed_frac),
                    str(rmhost_reads), str(rmhost_frac),
                    str(orphans_reads), str(orphans_frac)])
                outf.writelines(line+ '\n')

rule readcounts_graph:
    input:
        rules.readcounts.output
    output:
        join(PROJECT_DIR, "01_processing/readcounts.pdf")
    singularity: "shub://bhattlab/bhattlab_workflows:plotting"
    script:
        "scripts/plot_readcounts.R"

################################################################################
rule cleanup:
    input: expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX)
    output: join(PROJECT_DIR, "cleaned")
    params:
        rmdir_1 = join(PROJECT_DIR, '01_processing/01_dedup'),
        rmdir_2 = join(PROJECT_DIR, '01_processing/02_trimmed'),
    shell: """
        rm -f {params.rmdir_1}/*.fq.gz
        rm -f {params.rmdir_2}/*.fq.gz
        touch {output}
    """