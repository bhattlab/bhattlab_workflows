import re
from os.path import join, abspath
localrules: concat, plot_R_k21, plot_R_k31, plot_R_k51

preprocessing_dir = config["preprocessing_directory"]
outdir = config["output_directory"]
# convert outdir to absolute path
if outdir[0] == '~':
    outdir = expanduser(outdir)
outdir = abspath(outdir)

# coming from the output of the preprocessing pipeline
preprocessed_files=[f for f in os.listdir(join(preprocessing_dir)) if f.endswith('_1.fq.gz')]
sample_list = list(set([re.split('_1.fq.gz', i)[0] for i in preprocessed_files]))

rule all:
    input:
        expand(join(outdir, "01_concat_reads/{sample}_concat.fq"), sample=sample_list),
        expand(join(outdir, "02_trim_kmers/{sample}_concat.fq.abundtrim"), sample=sample_list),
        expand(join(outdir, "03_sourmash_signatures/{sample}.sig"), sample=sample_list),
        join(outdir, "04_sourmash_compare/compare_k21.csv"),
        join(outdir, "04_sourmash_compare/compare_k31.csv"),
        join(outdir, "04_sourmash_compare/compare_k51.csv"),
        join(outdir, "04_sourmash_compare/compare_k21_heatmap_complete.pdf"),
        join(outdir, "04_sourmash_compare/compare_k31_heatmap_complete.pdf"),
        join(outdir, "04_sourmash_compare/compare_k51_heatmap_complete.pdf"),

rule concat:
    input:
        fwd = join(preprocessing_dir, "{sample}_1.fq.gz"),
        rev = join(preprocessing_dir, "{sample}_2.fq.gz"),
        orp = join(preprocessing_dir, "{sample}_orphans.fq.gz"),
    output:
        allreads = join(outdir, "01_concat_reads/{sample}_concat.fq")
    shell: """
        zcat -f {input.fwd} {input.rev} {input.orp} > {output.allreads}
    """

rule trim_low_abund:
    input:
        rules.concat.output
    output:
        join(outdir, "02_trim_kmers/{sample}_concat.fq.abundtrim")
    params:
        outdir=join(outdir, "02_trim_kmers")
    threads:1
    resources:
        time=4,
        mem=40
    singularity: "docker://quay.io/biocontainers/khmer:3.0.0a3--py37hf484d3e_0"
    shell: """
        cd {params.outdir}
        trim-low-abund.py -C 3 -Z 18 -V -M 32e9 {input}
    """

rule compute:
    input:
        rules.trim_low_abund.output
    output:
        join(outdir, "03_sourmash_signatures/{sample}.sig")
    threads:1
    resources: 
        time=lambda wildcards, attempt: 12 * attempt,
    singularity: "docker://quay.io/biocontainers/sourmash:4.2.2--hdfd78af_0"
    shell: """
        sourmash compute --scaled 10000 \
        {input} -o {output} -k 21,31,51
    """

rule index:
    input:
        expand(join(outdir, "03_sourmash_signatures/{sample}.sig"), sample=sample_list)
    output:
        k21=join(outdir, "03_sourmash_signatures/{sample}_sig_db_k21.sbt.json"),
        k31=join(outdir, "03_sourmash_signatures/{sample}_sig_db_k31.sbt.json"),
        k51=join(outdir, "03_sourmash_signatures/{sample}_sig_db_k51.sbt.json")
    params:
        outfolder_k21="k21",
        outfolder_k31="k31",
        outfolder_k51="k51",
        filestr="*.sig",
        outdir=join(outdir, "03_sourmash_signatures")
    singularity: "docker://quay.io/biocontainers/sourmash:4.2.2--hdfd78af_0"
    shell: """
        cd {params.outdir}
        sourmash index {params.outfolder_k21} {params.filestr} -k21
        sourmash index {params.outfolder_k31} {params.filestr} -k31
        sourmash index {params.outfolder_k51} {params.filestr} -k51
    """

rule compare_k21:
    input:
        expand(join(outdir, "03_sourmash_signatures/{sample}.sig"), sample=sample_list)
    output:
        k21=join(outdir, "04_sourmash_compare/compare_k21.csv")
    params:
        filestr="*.sig",
        outfile=join(outdir, "04_sourmash_compare/k21"),
        sigdir=join(outdir, "03_sourmash_signatures")
    threads: 4
    resources:
        time=6,
        mem=256
    singularity: "docker://quay.io/biocontainers/sourmash:4.2.2--hdfd78af_0"
    shell: """
        cd {params.sigdir}
        sourmash compare {params.filestr} -o {params.outfile} --csv {output.k21} -k 21 -p {threads}
    """

rule compare_k31:
    input:
        expand(join(outdir, "03_sourmash_signatures/{sample}.sig"), sample=sample_list)
    output:
        k31=join(outdir, "04_sourmash_compare/compare_k31.csv")
    params:
        filestr="*.sig",
        outfile=join(outdir, "04_sourmash_compare/k31"),
        sigdir=join(outdir, "03_sourmash_signatures")
    threads: 4
    resources:
        time=6,
        mem=256
    singularity: "docker://quay.io/biocontainers/sourmash:4.2.2--hdfd78af_0"
    shell: """
        cd {params.sigdir}
        sourmash compare {params.filestr} -o {params.outfile} --csv {output.k31} -k 31 -p {threads}
    """

rule compare_k51:
    input:
        expand(join(outdir, "03_sourmash_signatures/{sample}.sig"), sample=sample_list)
    output:
        k51=join(outdir, "04_sourmash_compare/compare_k51.csv")
    params:
        filestr="*.sig",
        outfile=join(outdir, "04_sourmash_compare/k51"),
        sigdir=join(outdir, "03_sourmash_signatures")
    threads: 4
    resources:
        time=6,
        mem=256
    singularity: "docker://quay.io/biocontainers/sourmash:4.2.2--hdfd78af_0"
    shell: """
        cd {params.sigdir}
        sourmash compare {params.filestr} -o {params.outfile} --csv {output.k51} -k 51 -p {threads}
    """

rule plot_R_k21:
    input:
        join(outdir, "04_sourmash_compare/compare_k21.csv")
    output:
        join(outdir, "04_sourmash_compare/compare_k21_heatmap_complete.pdf")
    params:
        outdir = join(outdir, "04_sourmash_compare")
    singularity: "shub://bsiranosian/bens_1337_workflows:plotting"
    script: "scripts/heatmaps.R"

rule plot_R_k31:
    input:
        join(outdir, "04_sourmash_compare/compare_k31.csv")
    output:
        join(outdir, "04_sourmash_compare/compare_k31_heatmap_complete.pdf")
    params:
        outdir = join(outdir, "04_sourmash_compare")
    singularity: "shub://bsiranosian/bens_1337_workflows:plotting"
    script: "scripts/heatmaps.R"

rule plot_R_k51:
    input:
        join(outdir, "04_sourmash_compare/compare_k51.csv")
    output:
        join(outdir, "04_sourmash_compare/compare_k51_heatmap_complete.pdf")
    params:
        outdir = join(outdir, "04_sourmash_compare")
    singularity: "shub://bsiranosian/bens_1337_workflows:plotting"
    script: "scripts/heatmaps.R"

