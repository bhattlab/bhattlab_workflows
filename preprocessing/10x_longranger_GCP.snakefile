from os.path import join
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()


GS_PREFIX = "gbsc-gcp-lab-bhatt_user-bsiranos"
samples, *_ = GS.glob_wildcards(GS_PREFIX + '/raw_data_renamed/{sample}_S1_L001_R1_001.fastq.gz')
print(samples)

rule all:
    input:
        expand('barcoded_fastq_deinterleaved/{sample}_1.fq.gz', sample=samples)

rule longranger:
    input: 
        r1 = 'raw_data_renamed/{sample}_S1_L001_R1_001.fastq.gz',
        r2 = 'raw_data_renamed/{sample}_S1_L001_R2_001.fastq.gz'
    output: 'barcoded_fastq/{sample}_barcoded.fastq.gz'
    singularity: "docker://biocontainers/longranger:v2.2.2_cv2"
    threads: 15
    resources:
        mem=30,
        time=12
    params:
        fq_dir = join(GS_PREFIX, 'raw_data_renamed'),
        outdir = join(GS_PREFIX, '{sample}'),
    shell: """
        longranger basic --fastqs {params.fq_dir} --id {wildcards.sample} \
            --sample {wildcards.sample} --disable-ui --localcores={threads}
        mv {wildcards.sample}/outs/barcoded.fastq.gz {output}
    """

rule deinterleave:
    input:
        rules.longranger.output
    output:
        r1 = 'barcoded_fastq_deinterleaved/{sample}_1.fq.gz',
        r2 = 'barcoded_fastq_deinterleaved/{sample}_2.fq.gz'
    conda: "envs/pigz.yaml"
    threads: 7
    resources: 
        mem=8,
        time=12
    shell: """
        # code inspired by https://gist.github.com/3521724
        zcat {input} | paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" |
            pigz --best --processes {threads} > {output.r1}) | \
            cut -f 5-8 | tr "\t" "\n" | pigz --best --processes {threads} > {output.r2}
    """

'''
To start kubernetes cluser
# need to do something to enable autoscaling from here
# and maybe increase ephemeral storage for these jobs that use lots of tmp
export CLUSTER_NAME="snakemake-cluster-big3"
export ZONE="us-west1-b"
gcloud container clusters create $CLUSTER_NAME \
    --zone=$ZONE --num-nodes=50 \
    --machine-type="n1-standard-8" \
    --scopes storage-rw \
    --image-type=UBUNTU \
    --disk-size=500GB \
    --enable-autoscaling \
    --max-nodes=50 \
    --min-nodes=0
gcloud container node-pools create pool2 \
    --cluster $CLUSTER_NAME \
    --zone=$ZONE --num-nodes=17 \
    --machine-type="n1-standard-16" \
    --scopes storage-rw \
    --image-type=UBUNTU \
    --disk-size=500GB \
    --enable-autoscaling \
    --max-nodes=88 \
    --min-nodes=0

gcloud container clusters get-credentials --zone=$ZONE $CLUSTER_NAME
/labs/asbhatt/bsiranos/miniconda3/envs/mgwf/bin/snakemake -s preprocessing/10x_longranger.snakefile --default-remote-provider GS --default-remote-prefix gbsc-gcp-lab-bhatt_user-bsiranos  --use-singularity --kubernetes -j 99999 --use-conda

# THEN to delete...
gcloud container clusters delete --zone $ZONE $CLUSTER_NAME

'''