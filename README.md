# preprocessing

agreed steps:

1) qc - fastqc

2) trim - trim galore: all adapters, min read 60

3) deduplication - seakit: fwd and rev separately

4) sync - eli's script

5) post qc

6) multiqc

7) separate reference reads and give ref alignment stats


# assembly

spades

megahit

quast

# classification

kraken, barplots

metaphlan, heatmaps
