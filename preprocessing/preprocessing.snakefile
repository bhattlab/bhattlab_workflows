import re,os,subprocess
from os.path import join, expanduser, abspath

################################################################################
# specify project directories
DATA_DIR    = config["raw_reads_directory"]
PROJECT_DIR = config["output_directory"]
READ_SUFFIX = config["read_specification"]
EXTENSION   = config["extension"]

# convert PROJECT_DIR to absolute path
if PROJECT_DIR[0] == '~':
	PROJECT_DIR = expanduser(PROJECT_DIR)
PROJECT_DIR = abspath(PROJECT_DIR)

# get file names
FILES = [f for f in os.listdir(DATA_DIR) if f.endswith(tuple(['fastq.gz', 'fq.gz']))]
SAMPLE_PREFIX = list(set([re.split('_1.f|_2.f|_R1|_R2|_PE1|_PE2', i)[0] for i in FILES]))

################################################################################
localrules: assembly_meta_file, pre_multiqc, post_multiqc
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
	input:  join(DATA_DIR, "{sample}_{read}") + EXTENSION
	output: join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{sample}_{read}_fastqc.html")
	params:
		outdir = join(PROJECT_DIR, "01_processing/00_qc_reports/pre_fastqc/")
	threads: 1
	resources:
			time = 1,
			mem = 32
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
	shell: """
		multiqc --force {params.indir} -o {params.outdir}
	"""

################################################################################
rule trim_galore:
	input:
		fwd = join(DATA_DIR, "{sample}_") + READ_SUFFIX[0] + EXTENSION,
		rev = join(DATA_DIR, "{sample}_") + READ_SUFFIX[1] + EXTENSION
	output:
		fwd = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_") + READ_SUFFIX[0] + "_val_1.fq.gz",
		rev = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_") + READ_SUFFIX[1] + "_val_2.fq.gz",
		orp = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_unpaired.fq.gz")
	threads: 4
	resources:
		mem=32,
		time=lambda wildcards, attempt: attempt * 6
	params:
		orp_fwd = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_") + READ_SUFFIX[0] + "_unpaired_1.fq.gz",
		orp_rev = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_") + READ_SUFFIX[1] + "_unpaired_2.fq.gz",
		q_min   = config['trim_galore']['quality'],
		left    = config['trim_galore']['start_trim'],
		min_len = config['trim_galore']['min_read_length']
	shell: """
		mkdir -p {PROJECT_DIR}/01_processing/01_trimmed/
		trim_galore --quality {params.q_min} \
			--clip_R1 {params.left} --clip_R2 {params.left} \
			--length {params.min_len} \
			--output_dir {PROJECT_DIR}/01_processing/01_trimmed/ \
			--paired {input.fwd} {input.rev} \
			--retain_unpaired
		# merge unpaired reads - delete intermediate files
		zcat -f {params.orp_fwd} {params.orp_rev} | gzip >> {output.orp}
		rm {params.orp_fwd} {params.orp_rev}
	"""

################################################################################
rule dereplicate:
	input:
		fwd = rules.trim_galore.output.fwd,
		rev = rules.trim_galore.output.rev,
		orp = rules.trim_galore.output.orp
	output:
		fwd = join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_nodup_PE1.fastq"),
		rev = join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_nodup_PE2.fastq"),
		orp = join(PROJECT_DIR, "01_processing/02_dereplicate/{sample}_nodup_unpaired.fastq")
	threads: 2
	resources:
		mem=32,
		time=24
	shell: """
		mkdir -p {PROJECT_DIR}/01_processing/02_dereplicate/
		seqkit rmdup --by-seq {input.fwd} > {output.fwd} \
			-D {PROJECT_DIR}/01_processing/02_dereplicate/{wildcards.sample}_1_details.txt
		seqkit rmdup --by-seq {input.rev} > {output.rev} \
			-D {PROJECT_DIR}/01_processing/02_dereplicate/{wildcards.sample}_2_details.txt
		seqkit rmdup --by-seq {input.orp} > {output.orp} \
			-D {PROJECT_DIR}/01_processing/02_dereplicate/{wildcards.sample}_unpaired_details.txt
	"""
################################################################################
rule sync:
	input:
		rep_fwd  = rules.trim_galore.output.fwd,
		fwd = rules.dereplicate.output.fwd,
		rev = rules.dereplicate.output.rev
	output:
		fwd = join(PROJECT_DIR, "01_processing/03_sync/{sample}_1.fq"),
		rev = join(PROJECT_DIR, "01_processing/03_sync/{sample}_2.fq"),
		orp = join(PROJECT_DIR, "01_processing/03_sync/{sample}_orphans.fq")
	params:
		scripts_folder = config["scripts_dir"]
	shell: """
		mkdir -p {PROJECT_DIR}/01_processing/03_sync/
		{params.scripts_folder}/sync.py {input.rep_fwd} {input.fwd} {input.rev} {output.fwd} {output.rev} {output.orp}
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
	shell: """
		mkdir -p {PROJECT_DIR}/01_processing/04_host_align/
		# if an index needs to be built, use bwa index ref.fa
		# run on paired reads
		bwa mem -t {threads} {input.bwa_index} {input.fwd} {input.rev} | samtools view -bS - | samtools bam2fq -f 4 -1 {output.unmapped_1} -2 {output.unmapped_2} -
		# run on unpaired reads
		bwa mem -t {threads} {input.bwa_index} {input.orp} | samtools view -bS - | samtools bam2fq -f 4 - > {output.unmapped_orp}
	"""


################################################################################
rule rm_host_sync:
	input:
		rep_fwd  = rules.trim_galore.output.fwd,
		fwd = rules.rm_host_reads.output.unmapped_1,
		rev = rules.rm_host_reads.output.unmapped_2,
		orp = rules.rm_host_reads.output.unmapped_orp
	output:
		fwd = join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq"),
		rev = join(PROJECT_DIR, "01_processing/05_sync/{sample}_2.fq"),
		orp = join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq")
	params:
		scripts_folder = config["scripts_dir"]
	shell: """
		mkdir -p {PROJECT_DIR}/01_processing/05_sync/
		{params.scripts_folder}/sync.py {input.rep_fwd} {input.fwd} {input.rev} {output.fwd} {output.rev} {output.orp}
		# concatinate the orphans reads after host removal with the orphan reads from syncing prior to host removal
		cat {input.orp} >> {output.orp}

	"""

################################################################################
rule zip:
	input: join(PROJECT_DIR, "01_processing/05_sync/{sample}_{read}.fq")
	output: join(PROJECT_DIR, "01_processing/05_sync/{sample}_{read}.fq.gz")
	params:
		scripts_folder = config["scripts_dir"]
	threads: 8
	resources:
		time = 1,
		mem = 8
	shell:
		"pigz -8 {input}"

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
	shell: """
		mkdir -p {params.outdir}
		fastqc {input} -f fastq --outdir {params.outdir}
	"""

rule post_multiqc:
	input: expand(join(PROJECT_DIR,  "01_processing/00_qc_reports/post_fastqc/{sample}_{read}_fastqc.html"), sample=SAMPLE_PREFIX, read=['1', '2', 'orphans']),
	output: join(PROJECT_DIR,  "01_processing/00_qc_reports/post_multiqc/multiqc_report.html")
	params:
		indir = join(PROJECT_DIR,  "01_processing/00_qc_reports/post_fastqc"),
		outdir = join(PROJECT_DIR,  "01_processing/00_qc_reports/post_multiqc/")
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
				join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq"),
				join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_2.fq"),
				join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_orphans.fq")])]
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
				join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq"),
				join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_2.fq")])]
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
		raw = expand(join(DATA_DIR, "{sample}_") + READ_SUFFIX[0] + EXTENSION, sample=SAMPLE_PREFIX),
		trimmed = expand(join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_") + READ_SUFFIX[0] + "_val_1.fq.gz", sample=SAMPLE_PREFIX),
		dedup = expand(join(PROJECT_DIR, "01_processing/03_sync/{sample}_orphans.fq"), sample=SAMPLE_PREFIX),
		rmhost = expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq"), sample=SAMPLE_PREFIX),
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
				raw_file = join(DATA_DIR, sample + "_") + READ_SUFFIX[0] + EXTENSION
				trimmed_file = join(PROJECT_DIR, "01_processing/01_trimmed/" + sample + "_") + READ_SUFFIX[0] + "_val_1.fq.gz"
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
	params:
		scripts_folder = config["scripts_dir"]
	script:
		"scripts/plot_readcounts.R"

################################################################################
rule cleanup:
	input: expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX)
	output: join(PROJECT_DIR, "cleaned")
	shell: """
		rm {PROJECT_DIR}/01_processing/01_trimmed/*.fq.gz
		rm {PROJECT_DIR}/01_processing/02_dereplicate/*.fastq
		rm {PROJECT_DIR}/01_processing/03_sync/*.fq
		rm {PROJECT_DIR}/01_processing/04_host_align/*.fq
		touch {output}
	""""
