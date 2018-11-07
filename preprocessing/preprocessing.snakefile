import os,re

'''
Aim: A simple wrapper for metagenomics QC using paired end reads. To use this pipeline, edit parameters in the config.yaml, and specify the proper path to config file in the submission script.

This program runs under the assumption samples are named:
PREFIX_R1.fastq.gz and PREFIX_R2.fastq.gz.

This script will create the following folders:
PROJECT_DIR/qc/00_qc_reports/pre_fastqc
PROJECT_DIR/qc/00_qc_reports/post_fastqc
PROJECT_DIR/qc/01_trimmed
PROJECT_DIR/qc/02_dereplicate
PROJECT_DIR/qc/03_interleave
PROJECT_DIR/qc/04_host_align
'''

################################################################################
# specify project directories
DATA_DIR    = config["raw_reads_directory"]
PROJECT_DIR = config["output_directory"]

################################################################################
# get the names of the files in the directory
FILES = [f for f in os.listdir(DATA_DIR) if f.endswith(tuple(['fastq.gz', 'fq.gz']))]

SAMPLE_PREFIX = list(set([re.split('_1|_2|_R1|_R2', i)[0] for i in FILES]))
print(FILES)

################################################################################
# specify which rules do not need to be submitted to the cluster
localrules: sync, assembly_meta_file

rule all:
	input:
		expand(os.path.join(PROJECT_DIR,  "qc/00_qc_reports/pre_fastqc/{sample}_R{read}_fastqc.html"),
		      sample=SAMPLE_PREFIX, read=['1', '2']),
		expand(os.path.join(PROJECT_DIR, "qc/00_qc_reports/post_fastqc/{sample}_nodup_PE{read}_fastqc.html"),
		      sample=SAMPLE_PREFIX, read=['1','2']),
		expand(os.path.join(PROJECT_DIR, "qc/03_sync/{sample}_1.fastq"), sample=SAMPLE_PREFIX),
		#os.path.join(PROJECT_DIR, "qc/assembly_input.txt")
		#expand(os.path.join(PROJECT_DIR, "/qc/04_host_align/{sample}_{reference_name}_unmapped_{read}.fq"), sample=SAMPLE_PREFIX, reference_name=config['rm_host_reads']['host_pre'], read = ['1', '2'])

################################################################################
rule pre_fastqc:
	input:  os.path.join(DATA_DIR, "{sample}_R{read}.fastq.gz")
	output: os.path.join(PROJECT_DIR,  "qc/00_qc_reports/pre_fastqc/{sample}_R{read}_fastqc.html")
	threads: 1
	resources:
        	time = 1,
        	mem = 16
	shell: """
	   mkdir -p {PROJECT_DIR}/qc/00_qc_reports/pre_fastqc/
	   fastqc {input} --outdir {PROJECT_DIR}/qc/00_qc_reports/pre_fastqc/
	"""

################################################################################
rule trim_galore:
	input:
		fwd = os.path.join(DATA_DIR, "{sample}_R1.fastq.gz"),
		rev = os.path.join(DATA_DIR, "{sample}_R2.fastq.gz")
	output:
		fwd = os.path.join(PROJECT_DIR, "qc/01_trimmed/{sample}_val_1.fq.gz"),
		rev = os.path.join(PROJECT_DIR, "qc/01_trimmed/{sample}_val_2.fq.gz")
	threads: 4
	params:
		q_min   = config['trim_galore']['quality'],
		left    = config['trim_galore']['start_trim'],
		min_len = config['trim_galore']['min_read_length']
	resources:
        	time = 3,
        	mem = 24
	shell: """
		 mkdir -p {PROJECT_DIR}/qc/01_trimmed/
		 trim_galore --quality {params.q_min} \
			     --clip_R1 {params.left} --clip_R2 {params.left} \
			     --length {params.min_len} \
			     --output_dir {PROJECT_DIR}/qc/01_trimmed/ \
			     --paired {input.fwd} {input.rev}
	"""

################################################################################
rule dereplicate:
	input:
	 	fwd = rules.trim_galore.output.fwd,
		rev = rules.trim_galore.output.rev
	output:
	 	fwd = os.path.join(PROJECT_DIR, "qc/02_dereplicate/{sample}_nodup_PE1.fastq"),
		rev = os.path.join(PROJECT_DIR, "qc/02_dereplicate/{sample}_nodup_PE2.fastq")
	threads: 2
	resources:
        	time = 2,
        	mem = 16
	shell: """
		mkdir -p {PROJECT_DIR}/qc/02_dereplicate/
		seqkit rmdup -s {input.fwd} > {output.fwd}; seqkit rmdup -s {input.rev} > {output.rev}
	"""

################################################################################
rule post_fastqc:
	input:  rules.dereplicate.output
	output: os.path.join(PROJECT_DIR,  "qc/00_qc_reports/post_fastqc/{sample}_nodup_PE{read}_fastqc.html")
	threads: 1
	resources:
        	time = 1,
        	mem = 16
	shell: """
	   mkdir -p {PROJECT_DIR}/qc/00_qc_reports/post_fastqc/
	   fastqc {input} -f fastq --outdir {PROJECT_DIR}/qc/00_qc_reports/post_fastqc/
	 """

################################################################################
rule sync:
	input:
		rules.trim_galore.output[0],
		rules.dereplicate.output
	output:
		os.path.join(PROJECT_DIR, "qc/03_sync/{sample}_1.fastq"),
		os.path.join(PROJECT_DIR, "qc/03_sync/{sample}_2.fastq"),
		os.path.join(PROJECT_DIR, "qc/03_sync/{sample}_orphans.fastq")
	shell:
		"scripts/sync.sh {input} {output}"


################################################################################
rule assembly_meta_file:
	input:  rules.sync.output[0]
	output: os.path.join(PROJECT_DIR, "qc/assembly_input.txt")
	shell: """
		prefix=$(echo {input} | cut -d'_' -f1)
		files=$(find $(echo {input} | cut -d'_' -f1)* -maxdepth 0 -type f | tr '\n' ',' | sed 's/\(.*\),/\1 /' )
 		echo "$prefix	$files" >> {output}
	"""


################################################################################
rule align_host:
	input:
		fwd    = rules.sync.output[0],
		rev    = rules.sync.output[1],
		orphan = rules.sync.output[2],
		ref    = config['rm_host_reads']['host_genome']
	output: os.path.join(PROJECT_DIR, "qc/04_host_align/{sample}_{reference_name}.bam")
	threads: 8
	resources:
		time = 6,
		mem = 16
	shell: """
		mkdir -p {PROJECT_DIR}/qc/04_host_align/
		bowtie2 -x {input.ref} -1 {input.fwd} -2 {input.rev} -U {input.orphan} \
			--threads {threads} | samtools view -@ {threads} -bS - > {output}
	"""

################################################################################
rule rm_unmapped:
    input: rules.align_host.output
    output:
        unmapped_1 = os.path.join(PROJECT_DIR, "/qc/04_host_align/{sample}_{reference_name}_unmapped_1.fq"),
        unmapped_2 = os.path.join(PROJECT_DIR, "/qc/04_host_align/{sample}_{reference_name}_unmapped_2.fq"),
        unmapped_orphans1 = temp(os.path.join(PROJECT_DIR, "/qc/04_host_align/{sample}_{reference_name}_unmapped_orphans1.fq")),
        unmapped_orphans2 = temp(os.path.join(PROJECT_DIR, "/qc/04_host_align/{sample}_{reference_name}_unmapped_orphans2.fq")),
	unmapped_orphans = temp(os.path.join(PROJECT_DIR, "/qc/04_host_align/{sample}_{reference_name}_unmapped_orphans.fq"))
    threads: 4
    resources:
        time = 1,
        mem = 8
    shell: """
        samtools bam2fq -@ {threads} -f 4 -F 256 \
			-1 {output.unmapped_1} -2 {output.unmapped_2} \
			-s {output.unmapped_orphans} -0 {output.unmapped_orphans2} -n {input}
        cat {output.unmapped_orphans1} {output.unmapped_orphans2} > {output.unmapped_orphans}
    """
