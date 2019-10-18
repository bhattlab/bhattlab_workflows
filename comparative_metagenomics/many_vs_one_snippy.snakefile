# many metagenomic samples vs a single reference
# to search for SNPs, multiallelic sites

from os.path import join, abspath, expanduser

localrules: compile_snippy, count_depth, count_filtered_snps, filter_multiallelic_snps, filter_vcf, decompose_vcf

outdir = config['outdir']
# convert outdir to absolute path
if outdir[0] == '~':
    outdir = expanduser(outdir)
outdir = abspath(outdir)


def get_sample_reads(sample_file):
    sample_reads = {}
    paired_end = ''
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            sample = s[0]
            # paired end specified
            if (len(s)==3):
                reads = [s[1],s[2]]
                if paired_end != '' and not paired_end:
                    sys.exit('All samples must be paired or single ended.')
                paired_end = True
            # single end specified
            elif len(s)==2:
                sys.exit('Currently not supported single end reads')
                reads=s[1]
                if paired_end != '' and paired_end:
                    sys.exit('All samples must be paired or single ended.')
                paired_end = False
            if sample in sample_reads:
                raise ValueError("Non-unique sample encountered!")
            sample_reads[sample] = reads
    return (sample_reads, paired_end)

# read in sample info and reads from the sample_file
sample_reads, paired_end = get_sample_reads(config['sample_file'])
if paired_end:
    paired_string = '--paired'
else:
    paired_string = ''
sample_names = sample_reads.keys()


# get reference names
reference_name = config['reference_name']
reference_file = config['reference_file']

comparison_list = [(s, reference_name) for s in sample_names]
comparison_names = [s + "__" + reference_name for s in sample_names]
comparison_files = {s + "__" + reference_name:(sample_reads[s], reference_file) for s,r in comparison_list}
# print(comparison_names)
# print(comparison_files)

rule all:
    input:
        expand(join(outdir, "01_snippy/{name}/snps.vcf"), name=comparison_names),
        expand(join(outdir, "01_snippy/{name}/snps_filtered.vcf"), name=comparison_names),
        expand(join(outdir, "01_snippy/{name}/snps_decomposed_filtered.vcf"), name=comparison_names),
        expand(join(outdir, "01_snippy/{name}/multiallelic_snps_decomposed_filtered.vcf"), name=comparison_names),
        expand(join(outdir, "01_snippy/{name}/multiallelic_0.1_snps_decomposed_filtered.vcf"), name=comparison_names),
        expand(join(outdir, "01_snippy/{name}/snps.depth"), name=comparison_names),
        expand(join(outdir, "01_snippy/{name}/chromosome_min_depth_positions.txt"), name=comparison_names),
        expand(join(outdir, "01_snippy/{name}/snps_filtered_count_by_type.txt"), name=comparison_names),
        # expand(join(outdir, "top10_reference_SNP_total.txt"), name=comparison_names),

########################################################################
##### SNIPPY PIPELINE ##################################################
########################################################################
rule run_snippy:
    input:
        r1 = lambda wildcards: comparison_files[wildcards.name][0][0],
        r2 = lambda wildcards: comparison_files[wildcards.name][0][1],
        ref = reference_file
    output:
        join(outdir, '01_snippy/{name}/snps.vcf'),
        join(outdir, '01_snippy/{name}/snps.raw.vcf'),
        join(outdir, '01_snippy/{name}/snps.bam')
    params:
        outdir = join(outdir, '01_snippy/{name}')
    threads: 8
    singularity:
        "shub://bsiranosian/bens_1337_workflows:snippy"
    shell: """
        snippy --force --cpus {threads} --outdir {params.outdir} --ref {input.ref} --R1 {input.r1} --R2 {input.r2}
    """

rule samtools_depth:
    input:
        bam = join(outdir, '01_snippy/{name}/snps.bam')
    output:
        depth = join(outdir, '01_snippy/{name}/snps.depth'),
    singularity:
        "shub://bsiranosian/bens_1337_workflows:snippy"
    shell: """
        samtools depth -a {input.bam} > {output.depth}
    """

# compute depth stats for each chromosome - includes plasmids here
rule count_depth:
    input:
        join(outdir, '01_snippy/{name}/snps.depth')
    output:
        ten = join(outdir, '01_snippy/{name}/chromosome_min_depth_positions.txt'),
        fifty = join(outdir, '01_snippy/{name}/chromosome_min_depth_positions_50.txt'),
    params:
        min_depth_1 = 10,
        min_depth_2 = 50
    shell: """
        awk 'BEGIN {{OFS="\t"}} {{b[$1] ++}} $3 >= {params.min_depth_1} {{a[$1] ++}} END {{for (i in a) print i, b[i], a[i]}}' {input} > {output.ten}
        awk 'BEGIN {{OFS="\t"}} {{b[$1] ++}} $3 >= {params.min_depth_2} {{a[$1] ++}} END {{for (i in a) print i, b[i], a[i]}}' {input} > {output.fifty}
    """

rule filter_vcf:
    input:
        vcf = join(outdir, '01_snippy/{name}/snps.raw.vcf'),
        bam = join(outdir, '01_snippy/{name}/snps.bam')
    output:
        vcf_filt = join(outdir, '01_snippy/{name}/snps_filtered.vcf'),
        tab_filt = join(outdir, '01_snippy/{name}/snps_filtered.tab'),
        vcf_filt_fixed = join(outdir, '01_snippy/{name}/snps_fixed_filtered.vcf'),
        tab_filt_fixed = join(outdir, '01_snippy/{name}/snps_fixed_filtered.tab'),
    params:
        filter_string = 'QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0',
        filter_string_fixed = 'FMT/GT="1/1" && QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0',
        reference = join(outdir, '01_snippy/{name}/ref.fa'),
    singularity:
        "shub://bsiranosian/bens_1337_workflows:snippy"
    shell: """
        bcftools norm -f {params.reference} {input.vcf} | bcftools view --include \'{params.filter_string}\' > {output.vcf_filt}

        echo -e "SAMPLE\tCHROM\tPOS\tREF\tALT\tDP\tRO\tAO\tTYPE" > {output.tab}
        bcftools query  -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t%DP\t%RO\t%AO\t%TYPE\n' {output.vcf_filt} >> {output.tab}

        bcftools norm -f {params.reference} {input.vcf} | bcftools view --include \'{params.filter_string}\' > {output.vcf_filt}

        
    """

rule decompose_vcf:
    input:
        vcf = join(outdir, '01_snippy/{name}/snps.raw.vcf'),
        bam = join(outdir, '01_snippy/{name}/snps.bam')
    output:
        vcf_raw = join(outdir, '01_snippy/{name}/snps_decomposed_raw.vcf'),
        vcf_filt = join(outdir, '01_snippy/{name}/snps_decomposed_filtered.vcf'),
        tab = join(outdir, '01_snippy/{name}/snps_decomposed_filtered.tab'),
    params:
        filter_string = 'QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0',
        reference = join(outdir, '01_snippy/{name}/ref.fa'),
    singularity:
        "shub://bsiranosian/bens_1337_workflows:snippy"
    shell: """
        vt decompose -s {input.vcf} | vt decompose_blocksub -a - | bcftools norm -f {params.reference} -m -any > {output.vcf_raw}
        bcftools view --include \'{params.filter_string}\' {output.vcf_raw} > {output.vcf_filt}

        echo -e "SAMPLE\tCHROM\tPOS\tREF\tALT\tDP\tRO\tAO\tTYPE" > {output.tab}
        bcftools query  -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t%DP\t%RO\t%AO\t%TYPE\n' {output.vcf_filt} >> {output.tab}
    """

rule filter_multiallelic_snps:
    input:
        vcf_filt = join(outdir, '01_snippy/{name}/snps_decomposed_filtered.vcf'),
    output: 
        vcf_multiallelic = join(outdir, '01_snippy/{name}/multiallelic_snps_decomposed_filtered.vcf'),
        tab_multiallelic = join(outdir, '01_snippy/{name}/multiallelic_snps_decomposed_filtered.tab'),
        vcf_multiallelic_1 = join(outdir, '01_snippy/{name}/multiallelic_0.1_snps_decomposed_filtered.vcf'),
        tab_multiallelic_1 = join(outdir, '01_snippy/{name}/multiallelic_0.1_snps_decomposed_filtered.tab'),
    params:
        filter_string = 'QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0 && TYPE="SNP" && FMT/RO>=5 && FMT/AO>=5',
        filter_string_1 = 'QUAL>=100 && FMT/DP>=50 && (FMT/AO)/(FMT/DP)>=0.1 && (FMT/RO)/(FMT/DP)>=0.1 && TYPE="SNP"',
    singularity:
        "shub://bsiranosian/bens_1337_workflows:snippy"
    shell: """
        bcftools view --include \'{params.filter_string}\' {input.vcf_filt} > {output.vcf_multiallelic}
        echo -e "SAMPLE\tCHROM\tPOS\tREF\tALT\tDP\tRO\tAO\tTYPE" > {output.tab_multiallelic}
        bcftools query  -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t%DP\t%RO\t%AO\t%TYPE\n' {output.vcf_multiallelic} >> {output.tab_multiallelic}
        # and do the same for sites of AF at least 0.1
        bcftools view --include \'{params.filter_string_1}\' {input.vcf_filt} > {output.vcf_multiallelic_1}
        echo -e "SAMPLE\tCHROM\tPOS\tREF\tALT\tDP\tRO\tAO\tTYPE" > {output.tab_multiallelic_1}
        bcftools query  -f '[%SAMPLE]\t%CHROM\t%POS\t%REF\t%ALT\t%DP\t%RO\t%AO\t%TYPE\n' {output.vcf_multiallelic_1} >> {output.tab_multiallelic_1}

    """

rule count_filtered_snps:
    input:
        filtered = join(outdir, '01_snippy/{name}/snps_filtered.tab'),
        decomposed_filtered = join(outdir, '01_snippy/{name}/snps_decomposed_filtered.tab'),
        multiallelic_decomposed_filtered = join(outdir, '01_snippy/{name}/multiallelic_snps_decomposed_filtered.tab'),
        multiallelic_decomposed_filtered_1 = join(outdir, '01_snippy/{name}/multiallelic_0.1_snps_decomposed_filtered.tab'),
    output: 
        filtered = join(outdir, '01_snippy/{name}/snps_filtered_count_by_type.txt'),
        decomposed_filtered = join(outdir, '01_snippy/{name}/snps_decomposed_filtered_count_by_type.txt')
    shell: """
        tail -n +2 {input.filtered} | cut -f 2,9 | awk 'BEGIN {{FS="||"; OFS="\t"}} {{a[$1] ++}} END {{for (i in a) print i, a[i]}}' | sort -k1,1 -k2,2 > {output.filtered}
        tail -n +2 {input.decomposed_filtered} | cut -f 2,9 | awk 'BEGIN {{FS="||"; OFS="\t"}} {{a[$1] ++}} END {{for (i in a) print i, a[i]}}' | sort -k1,1 -k2,2 > {output.decomposed_filtered}
        echo -e $(tail -n +2 {input.multiallelic_decomposed_filtered} | head -n 1 | cut -f 2)"\tSNP;multiallelic\t"$(tail -n +2 {input.multiallelic_decomposed_filtered} | wc -l) >> {output.decomposed_filtered}
        echo -e $(tail -n +2 {input.multiallelic_decomposed_filtered_1} | head -n 1 | cut -f 2)"\tSNP;multiallelic_0.1\t"$(tail -n +2 {input.multiallelic_decomposed_filtered_1} | wc -l) >> {output.decomposed_filtered}
    """
# compile a table with counts for each sample
# 

# not implemented here yet
rule compile_snippy:
    input:
        expand(join(outdir, "01_snippy/{name}/snps_filtered_count_by_type.txt"), name=comparison_names),
        expand(join(outdir, "01_snippy/{name}/chromosome_min_depth_positions.txt"), name=comparison_names),
    output:
        join(outdir, 'top10_reference_SNP_total.txt')
    singularity:
        "shub://bsiranosian/bens_1337_workflows:snippy"
    params: 
        comparison_names = comparison_names, 
        reference_names = reference_name,
        basedir = join(outdir, "01_snippy"),
        outdir = outdir, 
        sample_name = sample_names
    script: "scripts/snippy_one_vs_many_plots.R"
