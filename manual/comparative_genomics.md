
## Comparative (meta)genomics pipelines
This is a collection of tools that are useful for comparing metagenomic datasets to each other and to references. They were also used in the Transmission of crAssphage paper (reference once it's available).

 1. Compare many metagenomic datasets (reads) to a single reference for SNP calling, multiallelic site identification
 2. Pairwise comparison of metagenomic assemblies of an microbe
 3. Compare a single metagenomic sample to a collection of references to identify what strain it's closest to

### SNP calling on many metagenomic datasets
`many_vs_one_snippy.snakefile`

This takes as input a set of sequencing reads and a reference genome. SNPs are called against the reference using [snippy](https://github.com/tseemann/snippy/). Variants are filtered for high quality. To compare many samples against each other, variants are normalized and decomposed using [vt](https://genome.sph.umich.edu/wiki/Vt). 

A heatmap of pairwise SNP similarity between samples is generated at the end. 
### Pairwise comparison of metagenomic assemblies of a single microbe
`compare_assembled_contigs.snakefile	`

This pipeline takes as input a set of metagenomic assemblies, filters for contigs >500bp in length, and aligns contigs to a reference genome. Contigs that align are then compared against each other, pairwise, using [nucmer](http://mummer.sourceforge.net/). 
<!--stackedit_data:
eyJoaXN0b3J5IjpbNDE3NTg4NzY0XX0=
-->