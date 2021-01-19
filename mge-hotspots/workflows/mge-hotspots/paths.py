import glob
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
            "prokka_annot",
            "prokka_genomes",
            "insert_beds",
            "genome_lengths",
            "sliding_windows",
            "window_counts",
            "window_significance_results",
            "window_significance_plots",
            "closest_genes",
            ):
            os.makedirs(self.__prefix+"/"+dir, exist_ok=True)

    @property
    def prokka_genomes_dir(self):
        return self.__prefix+"/prokka_genomes/"

    @property
    def prokka_annotations_dir(self):
        return self.__prefix+"/prokka_annot/"

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

    @property
    def closest_genes_dir(self):
        return self.__prefix+"/closest_genes/"

    def prokka_genome_dir(self, species):
       return self.prokka_genomes_dir+species

    def prokka_gff_files(self, species):
       return glob.glob(self.prokka_genomes_dir+species+"/*.gff")

    def prokka_bed_file(self, species):
       return self.prokka_annotations_dir+species+".coding_seqs.bed"

    def prokka_name_conversions(self):
       return self.__prefix+"/all_prokka_genomes_name_conversion.tsv"

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

    def top_windows_detail(self):
        return self.__prefix+"/all_top_windows_detailed.tsv"

    def top_windows_summary(self):
        return self.__prefix+"/all_top_windows_summarized.tsv"

    def top_windows_summary_bedfile(self):
        return self.__prefix+"/all_top_windows_summarized.bed"

    def top_windows_summary_merged_bedfile(self):
        return self.__prefix+"/all_top_windows_summarized.merged.bed"

    def top_windows_plot(self):
        return self.__prefix+"/genomewide_insertion_hotspots.pdf"

    def closest_genes_file(self, species):
        return self.closest_genes_dir+species+".closest_genes.tsv"

    def intergenic_hotspots_file(self, species):
        return self.closest_genes_dir+species+".intergenic_insertions.tsv"

    def all_hotspots_with_closest_genes(self):
        return self.__prefix+"/all_hotspots_with_closest_genes.tsv"

    def intergenic_hotspots_with_closest_genes(self):
        return self.__prefix+"/intergenic_hotspots_with_closest_genes.tsv"
