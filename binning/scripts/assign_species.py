#!/usr/bin/env/python
import sys
import operator
import os
import requests
import pandas as pd


krakf = snakemake.input["krak"]
binfolder = snakemake.params["binfolder"]
outf = snakemake.output[0]
custom_taxonomy = snakemake.params["custom_taxonomy"]

# TODO: work with tree taxonomy stored locally instead
# of requesting to ncbi for each...
# requests module necessary to get NCBI taxonomy information
def get_taxonomy_json(taxids):
    # print(taxids)
    url_base = 'http://taxonomy.jgi-psf.org/tax/id/ancestor/'
    url = url_base + ','.join(taxids)
    r = requests.get(url)
    try: 
        json=r.json()
    except:
        print("Failing taxonomy lookup on taxid:" + str(taxids))
        json=''
    return json


binfaifs = [os.path.join(binfolder, f) for f in os.listdir(binfolder) if os.path.splitext(f)[-1] == '.fai']
winning_margin = 66 #int(snakemake.input[1])

tig_species = {}
called_species = []

# read kraken file
krak = pd.read_csv(krakf, delimiter='\t', header=None, usecols=[1,2,3], index_col=0)
krak.columns = ['taxid', 'size']

df_columns = ['Bin', 'Size.Mb','lca_species', 'lca_level', 'lca_fraction',
              'best_species', 'best_level', 'best_fraction']
out_df = pd.DataFrame(columns=df_columns)
for binfaif in binfaifs:
    species_votes = {}

    bin_size = 0
    fai_df = pd.read_csv(binfaif, delimiter='\t', header=None)
    for index, row in fai_df.iterrows():
        contig = row[0]
        # number of votes = contig length
        votes = int(row[1])
        try:
            taxid = krak.loc[contig, 'taxid']
        except:
            taxid=0
        species_votes.setdefault(taxid, 0)
        species_votes[taxid] = species_votes[taxid] + votes
        try:
            bin_size += krak.loc[contig,'size']
        except:
            bin_size +=0
    species_votes_sorted = sorted(species_votes.items(), key=operator.itemgetter(1))
    total_votes = sum(species_votes.values())
    vote_thresh = total_votes * winning_margin / 100
    taxid_list = []
    accum_votes = 0
    size_mb = round(bin_size /1000000, 2)

    # take LCA of species that collectively make up 2/3 of the bin
    while accum_votes <= vote_thresh:
        species_tuple = species_votes_sorted.pop()
        taxid_list.append(str(species_tuple[0]))
        accum_votes += species_tuple[1]

    species_votes_sorted = sorted(species_votes.items(), key=operator.itemgetter(1))
    best_fraction = species_votes_sorted[-1][1] / total_votes
    best_taxid = species_votes_sorted[-1][0]
    # ensure best taxid is not zero, aka unclassified
    if best_taxid != 0:
        if not custom_taxonomy:
            best_json = get_taxonomy_json([str(best_taxid)])
            if 'name' in best_json and 'level' in best_json:
                best_species = best_json['name']
                best_level = best_json['level']
            else:
                best_species = "CLASSIFICATION ERROR"
                best_level = "CLASSIFICATION ERROR"
        # if working with custom taxonomy
        else:
            best_species = best_taxid
            best_level = 'Custom'
    else:
        best_species = 'Unclassified'
        best_level = 'Unclassified'

    # LCA calculations
    lca_fraction = accum_votes / total_votes
    # LCA name from ncbi request
    if not custom_taxonomy:
        lca_json = get_taxonomy_json(taxid_list)
        if 'name' in lca_json and 'level' in lca_json:
            lca_species = lca_json['name']
            lca_level = lca_json['level']
        else:
            lca_species = "CLASSIFICATION ERROR"
            lca_level = "CLASSIFICATION ERROR"
    # if custom taxonomy
    else:
        # best we can do now is report all of them 
        lca_species = ", ".join([str(a) for a in taxid_list])
        lca_level = "Custom"
        
    bin = binfaif.split("/")[-1].replace('.fai', '')
    
    fai_df = pd.DataFrame([bin, size_mb, lca_species, lca_level,
                           round(lca_fraction,3), best_species,
                           best_level, round(best_fraction,3)]).transpose()
    fai_df.columns = df_columns
    out_df = out_df.append(fai_df)
out_df.to_csv(outf, sep='\t', index=False)
