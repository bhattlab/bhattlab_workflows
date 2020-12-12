import os

def genome_filename(dir, genome):
    return dir +"/"+genome+".fna"

class Outputs:
    '''
    Outputs manages the directory and file hierarchy for the outputs
    of the series of scripts that calculate insertion hot spots.
    '''
    def __init__(self, prefix):
        self.__prefix = prefix
        for dir in (
            "insert_beds",
            "genome_lengths",
            "sliding_windows",
            "window_counts",
            "window_significance_results",
            "window_significance_plots",
            ):
            os.makedirs(self.__prefix+"/"+dir, exist_ok=True)

    @property
    def insertions_dir(self):
        return self.__prefix+"/insert_beds/"

    @property
    def lengths_dir(self):
        return self.__prefix+"/genome_lengths/"

    @property
    def lengths_tsv(self):
        return self.__prefix+"/genome_lengths/all_genomes.tsv"

    @property
    def windows_dir(self):
        return self.__prefix+"/sliding_windows/"

    @property
    def counts_dir(self):
        return self.__prefix+"/window_counts/"

    @property
    def significance_results_dir(self):
        return self.__prefix+"/window_significance_results/"

    @property
    def significance_plots_dir(self):
        return self.__prefix+"/window_significance_plots/"

    def insertions_filename(self, species, suffix):
        return self.insertions_dir+species+".insertions."+suffix+".bed"
    
    def length_filename(self, genome):
        return self.lengths_dir+genome+".genomeFile.txt"
    
    def windows_filename(self, genome):
        return self.windows_dir+genome+".sliding_windows.bed"

    def counts_filename(self, species):
        return self.counts_dir+species+".window_counts.tsv"

    def significance_results_filename(self, species, genome):
        return self.significance_results_dir+species+"."+genome+".window_signif.tsv"

    def significance_plots_filename(self, species, genome):
        return self.significance_plots_dir+species+"."+genome+".window_signif.pdf"
