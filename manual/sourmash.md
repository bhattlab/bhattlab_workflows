## Sourmash
Workflow for kmer-based comparison of sequencing reads. For the details on what sourmash and MinHash are, see the [sourmash docs](https://sourmash.readthedocs.io/en/latest/). Install and activate the `sourmash.yaml` conda environment in the `envs` folder.

There are two inputs required in the config file: the final output directory of the preprocessing pipeline (something like `preprocessing/01_processing/05_sync`) and the output directory. By default, kmer comparisons are calculated at k=21, k=31 and k=51. Output files include:
- Matrices of pairwise jaccard distances between all samples (`04_sourmash_compare/compare_k21.csv`)
- Clustered heatmaps of pairwise distances (`04_sourmash_compare/compare_k21_heatmap_ward.D2.pdf`)

Steps of the workflow include:
- Concatenate all reads into a single file
- Trim lowly abundant kmers
- Build MinHash sketches for each file
- Compare each pair of signatures to get pairwise Jaccard distances
- Downstream processing and plotting
<!--stackedit_data:
eyJoaXN0b3J5IjpbNDI1NTY4NTEzXX0=
-->