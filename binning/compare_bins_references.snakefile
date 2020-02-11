from os.path import join, abspath, expanduser, exists
import random
import string
import re
# compare bins against reference collection, do dotplots and etc
# takes in the end result of binnig 
# and edits the output tables?

# function from other binning pipeline to get assemblies and reads
# don't need unless we're doing snp comparisons
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

# for now all we need is sample names
def get_samples(sample_file):
    sample_list = []
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            sample = s[0]
            if sample in sample_list:
                print(sample)
                raise ValueError("Non-unique sample encountered!")
            sample_list.append(sample)
    return sample_list
   

# Read in sample and outdir from config file
sample_file = config['sample_file']
outdir = config['outdir_base']
# sample_reads, sample_assemblies = get_sample_assemblies_reads(sample_file)
# sample_list = list(sample_reads.keys())
sample_list = get_samples(sample_file)

print('##################################################################')
print(' SAMPLE LIST ')
print(sample_list)
print('##################################################################')

# convert outdir to absolute path
if outdir[0] == '~':
    outdir = expanduser(outdir)
outdir = abspath(outdir)

# keep top mash references
keep_mash_matches = 50
keep_fastani_matches = 10
# paramaters for fastani
minFrag = 5
fragLen = 2000
# settings of databases
genbank_sketch_file = "/labs/asbhatt/bsiranos/databases/genbank_nov2019/genbank_mash.msh"
# genbank_asm_file = "/labs/asbhatt/bsiranos/databases/genbank_nov2019/assembly_summary_chromosome_latest_combined.txt"
genbank_names_file = "/labs/asbhatt/bsiranos/databases/genbank_nov2019/fixed_organism_names.txt"
mag_sketch_file = "/labs/asbhatt/data/MAG_databases/mag_mash.msh"
# mag_asm_file = "/labs/asbhatt/bsiranos/databases/genbank_nov2019/assembly_summary_chromosome_latest_combined.txt"
mag_names_file = "/labs/asbhatt/data/MAG_databases/mag_mash_names.txt"
db_choices = ['genbank', 'mags']
sketch_files = {'genbank': genbank_sketch_file,
                'mags': mag_sketch_file}
names_files = {'genbank': genbank_names_file,
               'mags': mag_names_file}

# set binner used from the name of the folder
if (exists(join(outdir, sample_list[0], "DAS_tool_bins"))):
    binner_used = "DAS_tool"
else: 
    binner_used = "metabat"
bin_location_dict = {"DAS_tool": "DAS_tool_bins", 
                "metabat": "bins"}
bin_location = bin_location_dict[binner_used]

def get_bins(sample):
    print(sample)
    # remove unbinned
    bins = [b for b in glob_wildcards(join(outdir, sample, bin_location, "{bin}.fa")).bin if not(re.match('.*unbinned.*', b))]
    return bins

sample_bins = {s:get_bins(s) for s in sample_list}
print(sample_bins)

def aggregate_fastani(wildcards):
    return expand(join(outdir, "{sample}/reference_comparison/fastani_{choice}/{bin}_best.txt"),
        sample=wildcards.sample, 
        choice=wildcards.choice,
        bin=sample_bins[wildcards.sample])

def aggregate_nucmer(wildcards):
    return expand(join(outdir, "{sample}/reference_comparison/nucmer_{choice}/{bin}/report_stats.txt"),
        sample=wildcards.sample, 
        choice=wildcards.choice,
        bin=sample_bins[wildcards.sample])


rule all:
    input:
        # expand(join(outdir, "{sample}/reference_comparison/mash_{choice}/{bin}_mash_dist_best.tsv"), 
            # sample=sample_list, choice=db_choices, bin=lambda wildcards: sample_bins[wildcards.sample]),
        expand(join(outdir, "{sample}/reference_comparison/{choice}_done.txt"), 
            sample=sample_list, choice=db_choices),
        expand(join(outdir, "{sample}/reference_comparison/nucmer_top_{choice}.txt"), 
            sample=sample_list, choice=db_choices),


# compare to all reference sketch with MASH
rule compare_to_references_mash:
    input:
        bin = join(outdir, "{sample}", bin_location, "{bin}.fa"),
        sketch_file = lambda wildcards: sketch_files[wildcards.choice],
        names_file = lambda wildcards: names_files[wildcards.choice],
    output:
        best_dist = join(outdir, "{sample}/reference_comparison/mash_{choice}/{bin}_mash_dist_best.tsv"),
        best_files = join(outdir, "{sample}/reference_comparison/mash_{choice}/{bin}_mash_dist_best_files.tsv"),
        best_references = join(outdir, "{sample}/reference_comparison/mash_{choice}/{bin}_mash_dist_best_reference_names.tsv"),
    params:
        keep_mash_matches = keep_mash_matches
    threads: 1
    shell: """
        set +o pipefail
        mash dist -p {threads} {input.sketch_file} {input.bin} |sort -k 3,3 | head -n {params.keep_mash_matches} > {output.best_dist}
        # get matches from original reference list
        grep -f <(cut -f 1 {output.best_dist}) {input.names_file} > {output.best_references}
        cut -f 1 {output.best_dist} > {output.best_files}
    """

# then fastani with these genomes
rule compare_to_references_fastani:
    input:
        best_files = rules.compare_to_references_mash.output.best_files,
        bin = join(outdir, "{sample}", bin_location, "{bin}.fa"),
    output:
        result = join(outdir, "{sample}/reference_comparison/fastani_{choice}/{bin}.txt"),
    params:
        minFrag = minFrag,
        fragLen = fragLen,
    threads: 4
    resources:
        mem = 128,
        time = lambda wildcards, attempt: 2 ** attempt
    shell: """
        fastani --rl {input.best_files} -q {input.bin} -o {output} --minFrag {params.minFrag} --fragLen {params.fragLen} -t {threads}
    """

rule improve_fastani:
    input:
        result = rules.compare_to_references_fastani.output.result,
        reference_names_file = lambda wildcards: names_files[wildcards.choice],
    output:
        top_result = join(outdir, '{sample}/reference_comparison/fastani_{choice}/{bin}_best.txt'),
    params:
        keep_min_ani =80,
        keep_fastani_matches = keep_fastani_matches,
    threads: 1
    shell: """
        set +o pipefail 
        paste <(sort -k 3,3 -hr {input.result} | head -n {params.keep_fastani_matches} | awk '$3 >= {params.keep_min_ani} {{print}}') \
         <(sort -k 3,3 -hr {input.result} | head -n {params.keep_fastani_matches} | awk '$3 >= {params.keep_min_ani} {{print}}' | cut -f 2 | grep -f - {input.reference_names_file} | cut -f 2) \
          > {output.top_result}
    """


# nucmer dotplots and reports to best references
rule nucmer_from_fastani:
    input:
        best_references = rules.improve_fastani.output.top_result,
        bin = join(outdir, "{sample}", bin_location, "{bin}.fa"),
    output:
        join(outdir, "{sample}/reference_comparison/nucmer_{choice}/{bin}/report_stats.txt")
    params:
        outdir_nucmer = join(outdir, "{sample}/reference_comparison/nucmer_{choice}/{bin}"),
        outdir_reports = join(outdir, '{sample}/reference_comparison/nucmer_{choice}/{bin}/reports'),
        outdir_dotplots = join(outdir, '{sample}/reference_comparison/nucmer_{choice}/{bin}/dotplots'),
        outdir_dotplots_snps = join(outdir, '{sample}/reference_comparison/nucmer_{choice}/{bin}/dotplot_snps'),
        outdir_other = join(outdir, '{sample}/reference_comparison/nucmer_{choice}/{bin}/other'),
        tmpdir = join('/tmp/{bin}', ''.join(random.choice(string.ascii_letters) for i in range(6)))
    shell: """
        # only do something if theres lines in the input file
        if [[ $(wc -l <{input.best_references}) -ge 1 ]]; then
            mkdir -p {params.outdir_nucmer}
            mkdir -p {params.outdir_reports}
            mkdir -p {params.outdir_dotplots}
            mkdir -p {params.outdir_dotplots_snps}
            mkdir -p {params.outdir_other}
            cd {params.outdir_nucmer}
            while read line; do
                file=$(echo $line | cut -f 2 -d " ")
                filebase=$(basename "$file")
                name=$(echo $line | cut -f 6 -d " " | cut -d "|" -f 1)
                # echo "$line"
                echo starting: "$name"
                echo file: "$file"
                # if gzipped, unzip to tmp
                if [[ $file =~ \.gz$ ]]; then
                    # echo "UNZIPPING!"
                    mkdir -p {params.tmpdir}
                    gunzip -c $file > {params.tmpdir}/"${{filebase%.*}}"
                    file={params.tmpdir}/"${{filebase%.*}}"
                fi
                nucmer --prefix="$name" {input.bin} "$file"
                show-coords -rcl "$name".delta > "$name".coords
                # ensure we have data in the alignments before moving along with the pipeline
                if [[ $(wc -l <"$name".delta) -ge 3 ]]; then
                    mummerplot -f -l -png --large -p "$name"_dotplot "$name".delta
                    mummerplot -f -l -png --large -S -p "$name"_dotplot_snps "$name".delta
                    # snp detection
                    show-snps -Clr "$name".delta > "$name".snps
                    # reports
                    dnadiff -d "$name".delta -p "$name"
                else
                    # make some fakes
                    touch "$name"_dotplot.png
                    touch "$name"_dotplot_snps.png
                    echo -e "AvgIdentity 0\nTotalLength 0\n" > "$name".report
                fi
                # move output files around
                mv "$name".report {params.outdir_reports}
                mv "$name"_dotplot.png {params.outdir_dotplots}
                mv "$name"_dotplot_snps.png {params.outdir_dotplots_snps}
                rm "$name"* 
            done < {input.best_references} 

            for f in {params.outdir_reports}/*.report; do
                idty=$(grep "AvgIdentity" "$f" | head -n 1 | awk '$1=$1' | cut -f2 -d " ")
                ttl=$(grep "TotalLength" "$f" | head -n 1 | awk '$1=$1' | cut -f2 -d " ")
                comp_name1=$(basename $f)
                comp_name=$(echo $comp_name1 | sed "s/.report//g")
                printf '%s\t%s\t%s\n' "$comp_name" "$idty" "$ttl" >> {output}.tmp
            done
            sort -k2,2 -k3,3 -nr {output}.tmp > {output}
            rm {output}.tmp
        else
            echo "Skipping nucmer."
            touch {output}
        fi
    """

rule done_comparison:
    input:
        aggregate_fastani,
        aggregate_nucmer
    output:
        join(outdir, "{sample}/reference_comparison/{choice}_done.txt"), 
    shell: """
        touch {output}
    """


# append the top results to the final binning file 
rule top_nucmer:
    input:
        aggregate_nucmer
    output:
        join(outdir, "{sample}/reference_comparison/nucmer_top_{choice}.txt"), 
    params:
        nucmer_dir = join(outdir, "{sample}/reference_comparison/nucmer_{choice}"),
    shell: """
        cd {params.nucmer_dir}
        echo -e "bin\tgenome\tANI\talignment_length" > {output}
        for i in *; do 
            paste <(echo $i) <(sort -k2,2 -k3,3 -rn "$i"/report_stats.txt | head -n 1) >> {output}
        done
    """


# concatenate these all into one for all samples
rule concat_top_nucmer:
    input:
        expand(join(outdir, "{sample}/reference_comparison/nucmer_top_{choice}.txt"), sample=sample_list)
    output:
        join(outdir, "allsample_nucmer_top_{choice}.txt"), 
        outdir=outdir,
        choice=lambda wildcards: wildcards.choice,
    shell: """
        cd {params.outdir}
        echo -e "sample\tbin\tgenome\tANI\talignment_length" > {output}
        for i in $(ls -d */); do
            s=$(basename "$i")
            echo "$s"
            comp_file="$s"/reference_comparison/nucmer_top_{params.choice}.txt
            nl=$(tail -n +2 "$comp_file" | wc -l)
            # pase column of sample name and report info together
            # fix lines that end in tab (that have no alignment data) 
            # by adding in NA cols
            paste <(yes "$s" | head -n "$nl") <(tail -n +2 "$comp_file") | sed "s/\t$/\tNA\tNA\tNA/g" >> {output}
        done

    """

# run this last rule separately in case the code above wasnt ported to snakemake well
'''
out=allsample_nucmer_top_genbank.txt
echo -e "sample\tbin\tgenome\tANI\talignment_length" > "$out"
for i in $(ls -d */); do
    s=$(basename "$i")
    echo "$s"
    comp_file="$s"/reference_comparison/nucmer_top_genbank.txt
    nl=$(tail -n +2 "$comp_file" | wc -l)
    # pase column of sample name and report info together
    # fix lines that end in tab (that have no alignment data) 
    # by adding in NA cols
    paste <(yes "$s" | head -n "$nl") <(tail -n +2 "$comp_file") | sed "s/\t$/\tNA\tNA\tNA/g" >> "$out"
done

out=allsample_nucmer_top_mags.txt
echo -e "sample\tbin\tgenome\tANI\talignment_length" > "$out"
for i in $(ls -d */); do
    s=$(basename "$i")
    echo "$s"
    comp_file="$s"/reference_comparison/nucmer_top_mags.txt
    nl=$(tail -n +2 "$comp_file" | wc -l)
    # pase column of sample name and report info together
    # fix lines that end in tab (that have no alignment data) 
    # by adding in NA cols
    paste <(yes "$s" | head -n "$nl") <(tail -n +2 "$comp_file") | sed "s/\t$/\tNA\tNA\tNA/g" >> "$out"
done

'''