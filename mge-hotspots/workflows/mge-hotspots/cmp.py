
import pandas as pd

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
