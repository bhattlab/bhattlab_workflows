import sys

import pandas as pd

# example comparisons with rounding.

closestR = 'original/rare-insertions/closest_genes/1280--Staphylococcus_aureus.saur_GCF_000237265.closest_genes.tsv'
closestPy = 'rare_insertion_analysis/outputs/closest_genes/1280--Staphylococcus_aureus.closest_genes.tsv'

py = pd.read_csv(closestPy, delimiter='\t')
R = pd.read_csv(closestR, delimiter='\t')

py.iloc[:,[6]] = py.iloc[:,[6]].round(10)
R.iloc[:,[6]] = R.iloc[:,[6]].round(10)
print(py)

print(R)
print(py.compare(R))

sys.exit(0)

tsvPy = "./outputs/all_top_windows_summarized.tsv"
tsvR = "/Users/cnicolaou/Dropbox/github.com/bhattlab/bhattlab_workflows/mge-hotspots/14.GENOME_WIDE_INSERTION_ENRICHMENT/OUTPUT/all_top_windows_summarized.tsv"
py =  pd.read_csv(tsvPy, delimiter='\t')
R = pd.read_csv(tsvR, delimiter='\t')
py['pvalue'] = py['pvalue'].round(10)
R['pvalue'] = R['pvalue'].round(10)
print(py.compare(R))

bedPy = "./outputs/all_top_windows_summarized.bed"
bedR  = "/Users/cnicolaou/Dropbox/github.com/bhattlab/bhattlab_workflows/mge-hotspots/14.GENOME_WIDE_INSERTION_ENRICHMENT/OUTPUT/all_top_windows_summarized.bed"

py = pd.read_csv(bedPy, delimiter='\t')
py.columns = ['contig', 'start', 'end', 'pvalue']
R = pd.read_csv(bedR, delimiter='\t')
R.columns = ['contig', 'start', 'end', 'pvalue']
py['pvalue'] = py['pvalue'].round(10)
R['pvalue'] = R['pvalue'].round(10)
print(py.compare(R))

