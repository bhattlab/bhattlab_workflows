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
print(gz_ext)

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
			time = 1,
			mem = 32
	singularity: "docker://quay.io/biocontainers/fastqc:0.11.8--1"
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
	threads: 4
	resources:
		mem=32,
		time=lambda wildcards, attempt: attempt * 6
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

rule sort_trimmed_fwd:
	input:
		fwd = rules.trim_galore.output.fwd,
	output:
		fwd = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_" + READ_SUFFIX[0] + "_val_1_sorted.fq.gz"),
	resources:
		mem=64,
		time=lambda wildcards, attempt: attempt * 6
	shell: """
		export LC_ALL=C 
		zcat -f {input.fwd} | paste - - - - | sort -k1,1 | tr "\t" "\n" | gzip > {output.fwd} || test $? -eq 141	
		
		# remove old unsorted input
		# rm {input.fwd}
	"""

rule sort_trimmed_rev:
	input:
		rev = rules.trim_galore.output.rev,
	output:
		rev = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_" + READ_SUFFIX[1] + "_val_2_sorted.fq.gz"),
	resources:
		mem=64,
		time=lambda wildcards, attempt: attempt * 6
	shell: """
		export LC_ALL=C 
		zcat -f {input.rev} | paste - - - - | sort -k1,1 | tr "\t" "\n" | gzip > {output.rev} || test $? -eq 141	
		
		# remove old unsorted input
		# rm {input.rev}
	"""

rule sort_trimmed_orp:
	input:
		orp = rules.trim_galore.output.orp
	output:
		orp = join(PROJECT_DIR, "01_processing/01_trimmed/{sample}_unpaired_sorted.fq.gz")
	resources:
		mem=64,
		time=lambda wildcards, attempt: attempt * 6
	shell: """
		export LC_ALL=C 
		zcat -f {input.orp} | paste - - - - | sort -k1,1 | tr "\t" "\n" | gzip > {output.orp} || test $? -eq 141	
		
		# remove old unsorted input
		# rm {input.orp}
	"""


################################################################################
rule dereplicate_fwd:
	input:
		fwd = rules.sort_trimmed_fwd.output,
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
	shell: """
		mkdir -p {params.outdir}
		seqkit rmdup --by-seq {input.fwd} -D {params.details_fwd} > {output.fwd}
	"""

rule dereplicate_rev:
	input:
		rev = rules.sort_trimmed_rev.output,
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
	shell: """
		mkdir -p {params.outdir}
		seqkit rmdup --by-seq {input.rev} -D {params.details_rev} > {output.rev}
	"""

rule dereplicate_orp:
	input:
		orp = rules.sort_trimmed_orp.output
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
	shell: """
		mkdir -p {params.outdir}
		seqkit rmdup --by-seq {input.orp} -D {params.details_orp} > {output.orp}
	"""

################################################################################
rule sync:
	input:
		rep_fwd  = rules.sort_trimmed_fwd.output.fwd,
		fwd = rules.dereplicate_fwd.output.fwd,
		rev = rules.dereplicate_rev.output.rev,
		orp = rules.dereplicate_orp.output.orp
	output:
		fwd = join(PROJECT_DIR, "01_processing/03_sync/{sample}_1.fq"),
		rev = join(PROJECT_DIR, "01_processing/03_sync/{sample}_2.fq"),
		orp = join(PROJECT_DIR, "01_processing/03_sync/{sample}_orphans.fq")
	resources:
		time=12
	script:
		"scripts/sync.py"

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
		rep_fwd  = rules.sort_trimmed_fwd.output.fwd,
		fwd = rules.rm_host_reads.output.unmapped_1,
		rev = rules.rm_host_reads.output.unmapped_2,
		orp = rules.rm_host_reads.output.unmapped_orp
	output:
		fwd = join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq.gz"),
		rev = join(PROJECT_DIR, "01_processing/05_sync/{sample}_2.fq.gz"),
		orp = join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz")
	threads: 1
	resources:
		time = 12
	script:
		"scripts/sync.py"
	# shell:
	# 	" mkdir -p {PROJECT_DIR}/01_processing/05_sync/ " \
	# 	"{params.scripts_folder}/sync.py {input.rep_fwd} {input.fwd} {input.rev} " \
	# 	"{PROJECT_DIR}/01_processing/05_sync/{wildcards.sample}_1.fq " \
	# 	"{PROJECT_DIR}/01_processing/05_sync/{wildcards.sample}_2.fq {PROJECT_DIR}/01_processing/05_sync/{wildcards.sample}_orphans.fq " \
	# 	"# concatenate the orphans reads after host removal with the orphan reads from syncing prior to host removal " \
	# 	"cat {input.orp} >> {PROJECT_DIR}/01_processing/05_sync/{wildcards.sample}_orphans.fq " \
	# 	"pigz -8 {PROJECT_DIR}/01_processing/05_sync/{wildcards.sample}_1.fq " \
	# 	"pigz -8 {PROJECT_DIR}/01_processing/05_sync/{wildcards.sample}_2.fq " \
	# 	"pigz -8 {PROJECT_DIR}/01_processing/05_sync/{wildcards.sample}_orphans.fq"

################################################################################
# rule zip:
# 	input: join(PROJECT_DIR, "01_processing/05_sync/{sample}_{read}.fq")
# 	output: join(PROJECT_DIR, "01_processing/05_sync/{sample}_{read}.fq.gz")
# 	params:
# 		scripts_folder = config["scripts_dir"]
# 	threads: 8
# 	resources:
# 		time = 1,
# 		mem = 8
# 	shell:
# 		"pigz -8 {input}"

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
