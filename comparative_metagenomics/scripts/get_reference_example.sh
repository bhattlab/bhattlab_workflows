cd /labs/asbhatt/bsiranos/transmit_crass/compare_aligned_assembled_reads/references/enterococcus_faecalis_genomes
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Enterococcus_faecalis/assembly_summary.txt
cut -f 20 assembly_summary.txt | tail -n +3  | xargs -n1 -P 4 -I foo wget foo/*[^m]_genomic.fna.gz
ls *.fna.gz | xargs -n1 -P4 gunzip 
cat *.fna > ../enterococcus_faecalis_composite.fasta