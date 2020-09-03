## Antibiotic Resistance Gene (ARG) detection
#### Using the RGI tool from CARD

These pipelines allow you to use the ARG detection tools from CARD easily. First, familiarize yourself with the tools [here](https://github.com/arpcard/rgi). There are two detection tools available: one that works on assembles genomes and one that works on read data. 

### Installation
Follow the instructions at the RGI github page linked above for installation and dependencies. Recommend creating a conda environment for this and all the dependencies. 

### Database configuration
You need to set up rgi to use a database of antibiotic resistance genes. Either download a new database and set it up manually (instructions are at the github), or use the database in Ben's folder. Instructions for configuring the database using Ben's folder: 
```
db_folder="/oak/stanford/scg/lab_asbhatt/bsiranos/databases/CARD"
db_version="3.1.0"
rgi load --card_json "$db_folder"/"$db_version"/card.json 
rgi load -i "$db_folder"/"$db_version"/card.json --card_annotation "$db_folder"/"$db_version"/card_database_v"$db_version".fasta 
rgi load --wildcard_annotation "$db_folder"/"$db_version"/wildcard/wildcard_database_v"$db_version".fasta --wildcard_index "$db_folder"/"$db_version"/wildcard/index-for-model-sequences.txt --card_annotation "$db_folder"/"$db_version"/card_database_v"$db_version".fasta 
rgi load --kmer_database "$db_folder"/"$db_version"/wildcard/61_kmer_db.json --amr_kmers "$db_folder"/"$db_version"/wildcard/all_amr_61mers.txt --kmer_size 61 --debug > kmer_load.61.log 2>&1
# check that this correctly reports the version above
rgi database -v
```

### RGI main 
This is the main RGI tool, which identifies ARGs found in an assembly. Use `arg_detection/rgi.snakefile` for this. See the corresponding config file for the necessary inputs, just a sample table and output directory. 

RGI identifies ARGs on contigs, but this doesn't tell anything about the organism that contains the ARG. If you have binning data, you can use script in `arg_detection/process_rgi_add_bins.R` to incorporate binning info into the RGI calls. This script is not polished and requires editing to use with your data, contact Ben for help. 

### RGI bwt
This is an experimental pipeline from CARD that maps reads against the database of ARGs. Use `arg_detection/rgi_bwt.snakefile` and `arg_detection/config_rgi_bwt.yaml`. This pipeline uses a table of samples to forward and reverse reads as input instead of assemblies. 
