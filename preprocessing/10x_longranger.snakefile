import re,os,subprocess
from os.path import join, expanduser, abspath
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()


GS_PREFIX = "gbsc-gcp-lab-bhatt_user-bsiranos"
# GENOME_FA =  GS.remote(f"{GS_PREFIX}/griffithlab_brain_vs_uhr/GRCh38_Ens87_chr22_ERCC/chr22_ERCC92.fa")
# GENOME_GTF = GS.remote(f"{GS_PREFIX}/griffithlab_brain_vs_uhr/GRCh38_Ens87_chr22_ERCC/genes_chr22_ERCC92.gtf")
# HISAT2_INDEX_PREFIX = "hisat2_index/chr22_ERCC92"
print('grobbing')
samples, *_ = GS.glob_wildcards(GS_PREFIX + '/raw_data_renamed/{sample}_1.fq.gz')
print(samples)

rule all:
    input:
        expand('barcoded_fastq/{sample}/barcoded_1.fq.gz', sample=samples)

rule longranger:
    input: 
        r1 = 'raw_data_renamed/{sample}_1.fq.gz',
        r2 = 'raw_data_renamed/{sample}_2.fq.gz'
    output: 'barcoded_fastq/{sample}/outs/barcoded.fastq.gz'
    singularity: "docker://biocontainers/longranger:v2.2.2_cv2"
    threads: 8
    resources:
        mem=30,
        time=12
    params:
        fq_dir = 'raw_data_renamed'
    shell: """
        longranger basic --fastqs {params.fq_dir}/{wildcards.sample} --id {wildcards.sample} --sample {wildcards.sample} --disable-ui --localcores={threads}
    """

rule deinterleave:
    input:
        rules.longranger.output
    output:
        r1 = 'barcoded_fastq/{sample}/barcoded_1.fq.gz',
        r2 = 'barcoded_fastq/{sample}/barcoded_2.fq.gz'
    params:
        sample_dir = 'barcoded_fastq/{sample}'
    threads: 4
    resources: 
        mem=8,
        time=12
    shell: """
        zcat {input} | paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes {threads} > {output.r1}) | cut -f 5-8 | tr "\t" "\n" | pigz --best --processes {threads} > {output.r2}
    """
'''
longranger basic  --fastqs /labs/asbhatt/bsiranos/transmission_10x/raw_data_renamed/bas_bb_spri --id bas_bb_spri --sample bas_bb_spri --disable-ui --localcores=16

cd /labs/asbhatt/bsiranos/transmission_10x/barcoded_fastq/bas_bb_spri/outs
; zcat barcoded.fastq.gz | /labs/asbhatt/bsiranos/transmission_10x/deinterleave_fastq.sh barcoded_1.fq.gz barcoded_2.fq.gz compress"

ln -s /labs/asbhatt/bsiranos/transmission_10x/barcoded_fastq/bas_bb_spri/outs/barcoded_1.fq.gz 
/labs/asbhatt/bsiranos/transmission_10x/raw_links/bas_bb_spri_1.fq.gz
ln -s /labs/asbhatt/bsiranos/transmission_10x/barcoded_fastq/bas_bb_spri/outs/barcoded_2.fq.gz /labs/asbhatt/bsiranos/transmission_10x/raw_links/bas_bb_spri_2.fq.gz


################################################################################
# specify project directories
DATA_DIR    = config["raw_reads_directory"]
PROJECT_DIR = config["output_directory"]
READ_SUFFIX = config["read_specification"]
EXTENSION   = config["extension"]
# if gzipped, set this. otherwise not
gz_ext = '.gz' if EXTENSION.endswith('.gz') else ''
# print(gz_ext)

# convert PROJECT_DIR to absolute path
if PROJECT_DIR[0] == '~':
    PROJECT_DIR = expanduser(PROJECT_DIR)
PROJECT_DIR = abspath(PROJECT_DIR)

# get file names
FILES = [f for f in os.listdir(DATA_DIR) if f.endswith(EXTENSION)]
SAMPLE_PREFIX = list(set([re.split('|'.join(['_' + a + '\.' for a in READ_SUFFIX]), i)[0] for i in FILES]))

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
    input:  join(DATA_DIR, "{sample}_{read}" + EXTENSION)
    output: join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{sample}_{read}_fastqc.html")
    params:
        outdir = join(PROJECT_DIR, "01_processing/00_qc_reports/pre_fastqc/")
    threads: 1
    resources:
            time = 6,
            mem = 32
    singularity: "docker://quay.io/biocontainers/fastqc:0.11.8--1"
    benchmark: join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{sample}_{read}_time.txt")
    shell: """
        mkdir -p {params.outdir}
        fastqc {input} --outdir {params.outdir}
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

################################################################################
rule trim_galore:
    input:
        fwd = join(DATA_DIR, "{sample}_" + READ_SUFFIX[0] + EXTENSION),
        rev = join(DATA_DIR, "{sample}_" + READ_SUFFIX[1] + EXTENSION)
    output:
        fwd = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_" + READ_SUFFIX[0] + "_val_1.fq" + gz_ext),
        rev = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_" + READ_SUFFIX[1] + "_val_2.fq" + gz_ext),
        orp = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_unpaired.fq" + gz_ext)
    threads: 8
    resources:
        mem=32,
        time=lambda wildcards, attempt: attempt * 24
    params:
        orp_fwd = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_" + READ_SUFFIX[0] + "_unpaired_1.fq" + gz_ext),
        orp_rev = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_" + READ_SUFFIX[1] + "_unpaired_2.fq" + gz_ext),
        # output_fwd_pre_gz = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_") + READ_SUFFIX[0] + "_val_1.fq",
        # output_rev_pre_gz = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_") + READ_SUFFIX[1] + "_val_2.fq",
        q_min   = config['trim_galore']['quality'],
        left    = config['trim_galore']['start_trim'],
        min_len = config['trim_galore']['min_read_length'],
        outdir  = join(PROJECT_DIR, "01_processing/01_trimmed/"),
        gz_output = str(gz_ext == '.gz').lower()
    singularity: "docker://quay.io/biocontainers/trim-galore:0.5.0--0"
    benchmark: join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_time.txt")
    shell: """
        mkdir -p {params.outdir}
        trim_galore --quality {params.q_min} \
            --clip_R1 {params.left} --clip_R2 {params.left} \
            --length {params.min_len} \
            --output_dir {PROJECT_DIR}/01_processing/01_trimmed/ \
            --paired {input.fwd} {input.rev} \
            --retain_unpaired

        # if output is gz, merge unpaired and gzip
        if {params.gz_output}; then
            zcat -f {params.orp_fwd} {params.orp_rev} | pigz -p {threads} > {output.orp}
        else
            zcat -f {params.orp_fwd} {params.orp_rev} > {output.orp}
        fi
        # delete intermediate files
        rm {params.orp_fwd} {params.orp_rev}
    """

################################################################################
rule dereplicate_fwd:
    input:
        fwd = rules.trim_galore.output.fwd,
    output:
        fwd = join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_nodup_PE1.fastq"),
    params:
        outdir = join(PROJECT_DIR, "01_processing/02_dereplicate/"),
        details_fwd = join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_1_details.txt"),
    threads: 2
    resources:
        mem=32,
        time=24
    singularity: "docker://quay.io/biocontainers/seqkit:0.10.1--1"
    benchmark: join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_PE1_time.txt")
    shell: """
        mkdir -p {params.outdir}
        seqkit rmdup --by-seq {input.fwd} -D {params.details_fwd} > {output.fwd}
    """

rule dereplicate_rev:
    input:
        rev = rules.trim_galore.output.rev,
    output:
        rev = join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_nodup_PE2.fastq"),
    params:
        outdir = join(PROJECT_DIR, "01_processing/02_dereplicate/"),
        details_rev = join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_2_details.txt"),
    threads: 2
    resources:
        mem=32,
        time=24
    singularity: "docker://quay.io/biocontainers/seqkit:0.10.1--1"
    benchmark: join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_PE2_time.txt")
    shell: """
        mkdir -p {params.outdir}
        seqkit rmdup --by-seq {input.rev} -D {params.details_rev} > {output.rev}
    """

rule dereplicate_orp:
    input:
        orp = rules.trim_galore.output.orp
    output:
        orp = join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_nodup_unpaired.fastq")
    params:
        outdir = join(PROJECT_DIR, "01_processing/02_dereplicate/"),
        details_orp = join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_unpaired_details.txt")
    threads: 2
    resources:
        mem=32,
        time=24
    singularity: "docker://quay.io/biocontainers/seqkit:0.10.1--1"
    benchmark: join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_unpaired_time.txt")
    shell: """
        mkdir -p {params.outdir}
        seqkit rmdup --by-seq {input.orp} -D {params.details_orp} > {output.orp}
    """

################################################################################
rule sync:
    input:
        fwd = rules.dereplicate_fwd.output.fwd,
        rev = rules.dereplicate_rev.output.rev,
        orp = rules.dereplicate_orp.output.orp
    output:
        fwd = join(PROJECT_DIR, "01_processing/03_sync/{sample}_1.fq"),
        rev = join(PROJECT_DIR, "01_processing/03_sync/{sample}_2.fq"),
        orp = join(PROJECT_DIR, "01_processing/03_sync/{sample}_orphans.fq")
    resources:
        time = 12,
        mem = 128
    params: 
        outfile_paired_1 = rules.dereplicate_fwd.output.fwd + ".paired.fq",
        outfile_single_1 = rules.dereplicate_fwd.output.fwd + ".single.fq",
        outfile_paired_2 = rules.dereplicate_rev.output.rev + ".paired.fq",
        outfile_single_2 = rules.dereplicate_rev.output.rev + ".single.fq",
    singularity: "shub://bsiranosian/bens_1337_workflows:fastq-pair"
    benchmark: join(PROJECT_DIR, "01_processing/03_sync/{sample}_time.txt")
    shell: """
        fastq_pair {input.fwd} {input.rev}
        # move files to expected locations
        mv {params.outfile_paired_1} {output.fwd}
        mv {params.outfile_paired_2} {output.rev}
        cat {input.orp} {params.outfile_single_1} {params.outfile_single_2} > {output.orp}
        rm {params.outfile_single_1} {params.outfile_single_2}
        """

################################################################################
rule rm_host_reads:
    input:
        bwa_index = config['rm_host_reads']['host_genome'],
        fwd       = rules.sync.output.fwd,
        rev       = rules.sync.output.rev,
        orp       = rules.sync.output.orp
    output:
        unmapped_1 = join(PROJECT_DIR, "01_processing/04_host_align/{sample}_rmHost_1.fq"),
        unmapped_2 = join(PROJECT_DIR, "01_processing/04_host_align/{sample}_rmHost_2.fq"),
        unmapped_orp = join(PROJECT_DIR, "01_processing/04_host_align/{sample}_rmHost_unpaired.fq")
    threads: 4
    resources:
        mem=32,
        time=24
    singularity: "shub://bhattlab/bhattlab_workflows:align"
    benchmark: join(PROJECT_DIR, "01_processing/04_host_align/{sample}_time.txt")
    shell: """
        mkdir -p {PROJECT_DIR}/01_processing/04_host_align/
        # if an index needs to be built, use bwa index ref.fa
        # run on paired reads
        bwa mem -t {threads} {input.bwa_index} {input.fwd} {input.rev} | samtools view -b - | samtools fastq -t -T BX -f 4 -1 {output.unmapped_1} -2 {output.unmapped_2} -
        # run on unpaired reads
        bwa mem -t {threads} {input.bwa_index} {input.orp} | samtools view -bS - | samtools fastq -t -T BX -f 4 - > {output.unmapped_orp}
    """

################################################################################
rule rm_host_sync:
    input:
        # rep_fwd  = rules.sort_trimmed_fwd.output.fwd,
        fwd = rules.rm_host_reads.output.unmapped_1,
        rev = rules.rm_host_reads.output.unmapped_2,
        orp = rules.rm_host_reads.output.unmapped_orp
    output:
        fwd = join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq"),
        rev = join(PROJECT_DIR, "01_processing/05_sync/{sample}_2.fq"),
        orp = join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq")
    threads: 1  
    resources:
        time = 12,
        mem = 128
    params:
        outfile_paired_1 = rules.rm_host_reads.output.unmapped_1 + ".paired.fq",
        outfile_single_1 = rules.rm_host_reads.output.unmapped_1 + ".single.fq",
        outfile_paired_2 = rules.rm_host_reads.output.unmapped_2 + ".paired.fq",
        outfile_single_2 = rules.rm_host_reads.output.unmapped_2 + ".single.fq",
    singularity: "shub://bsiranosian/bens_1337_workflows:fastq-pair"
    benchmark: join(PROJECT_DIR, "01_processing/05_sync/{sample}_time.txt")
    shell: """
        fastq_pair {input.fwd} {input.rev}
        # move files to expected locations
        mv {params.outfile_paired_1} {output.fwd}
        mv {params.outfile_paired_2} {output.rev}
        cat {input.orp} {params.outfile_single_1} {params.outfile_single_2} > {output.orp}
        rm {params.outfile_single_1} {params.outfile_single_2}
        """

###############################################################################
rule zip:
    input: 
        fwd = join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq"),
        rev = join(PROJECT_DIR, "01_processing/05_sync/{sample}_2.fq"),
        orp = join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq"),
    output:
        fwd = join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq.gz"),
        rev = join(PROJECT_DIR, "01_processing/05_sync/{sample}_2.fq.gz"),
        orp = join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"),
    threads: 8
    resources:
        time = 1,
        mem = 32
    benchmark: join(PROJECT_DIR, "01_processing/05_sync/{sample}_zip_time.txt")
    shell: """
        pigz -p 8 {input.fwd}
        pigz -p 8 {input.rev}
        pigz -p 8 {input.orp}
        # rm {input}
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
        trimmed = expand(join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_" + READ_SUFFIX[0] + "_val_1.fq" + gz_ext), sample=SAMPLE_PREFIX),
        dedup = expand(join(PROJECT_DIR, "01_processing/03_sync/{sample}_orphans.fq"), sample=SAMPLE_PREFIX),
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
            outf.writelines('Sample\traw_reads\ttrimmed_reads\ttrimmed_frac\tdeduplicated_reads\tdeduplicated_frac\thost_removed_reads\thost_removed_frac\torphan_reads\torphan_frac\n')
            for sample in SAMPLE_PREFIX:
                raw_file = join(DATA_DIR, sample + "_" + READ_SUFFIX[0] + EXTENSION)
                trimmed_file = join(PROJECT_DIR, "01_processing/01_trimmed/" + sample + "_" + READ_SUFFIX[0] + "_val_1.fq" + gz_ext)
                dedup_file = join(PROJECT_DIR, "01_processing/03_sync/" + sample + "_1.fq")
                rmhost_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq")
                orphans_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_orphans.fq")

                raw_reads = int(file_len(raw_file) / 4)
                trimmed_reads = int(file_len(trimmed_file) / 4)
                dedup_reads = int(file_len(dedup_file) / 4)
                rmhost_reads = int(file_len(rmhost_file) / 4)
                orphans_reads = int(file_len(orphans_file) / 4)

                trimmed_frac = round(trimmed_reads / float(raw_reads), 3)
                dedup_frac = round(dedup_reads / float(raw_reads), 3)
                rmhost_frac = round(rmhost_reads / float(raw_reads), 3)
                orphans_frac = round(orphans_reads / float(raw_reads), 3)

                line = '\t'.join([sample, str(raw_reads),
                    str(trimmed_reads), str(trimmed_frac),
                    str(dedup_reads), str(dedup_frac),
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
    shell:
        "rm {PROJECT_DIR}/01_processing/01_trimmed/*.fq.gz / "
        "rm {PROJECT_DIR}/01_processing/02_dereplicate/*.fastq / "
        "rm {PROJECT_DIR}/01_processing/03_sync/*.fq / "
        "rm {PROJECT_DIR}/01_processing/04_host_align/*.fq / "
        "touch {output}"
'''