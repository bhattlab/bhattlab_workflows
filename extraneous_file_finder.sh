Folders you can most likely delete, which contain extra files and snakemake logs

find . -print0 > find_tmp.txt

echo "# finding likely extra directories" > extra_files.txt
echo "# .snakemake folders" >> extra_files.txt
# .snakemake directories
grep "\.snakemake$" find_tmp.txt  >> extra_files.txt

# Spades assembly temp files and graphs
echo "# extra spades assembly directories" >> extra_files.txt
grep "spades.*K[0-9][0-9]$" find_tmp.txt  >> extra_files.txt
grep "spades.*misc$" find_tmp.txt  >> extra_files.txt
grep "spades.*corrected$" find_tmp.txt  >> extra_files.txt
grep "spades.*tmp$" find_tmp.txt  >> extra_files.txt

# preprocessing extra files
echo "# extra preprocessing files" >> extra_files.txt
grep "preprocessing.*01_trimmed.*.fq.gz$" find_tmp.txt >> extra_files.txt
grep "preprocessing.*02_dereplicate.*fastq$"  find_tmp.txt >> extra_files.txt
grep "preprocessing.*03_sync.*fq$"  find_tmp.txt >> extra_files.txt
grep "preprocessing.*04_.*fq$"  find_tmp.txt >> extra_files.txt

