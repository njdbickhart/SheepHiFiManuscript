#!/usr/bin/env python
# Originally from Magpy. Altered to fit this pipeline

import sys

from ete3 import NCBITaxa

# get NCBI taxonomy object
ncbi = NCBITaxa()

# open the file
checkm_file = snakemake.input[0]
outfile = snakemake.output[0]

# skip three lines
for i in range(3):
    checkm_file.readline()

# print titles for the output
titles = ["bin_id",
		"marker_lineage",
		"ml_uid",
		"n_genomes",
		"n_markers",
		"n_marker_sets",
		"n_0",
		"n_1",
		"n_2",
		"n_3",
		"n_4",
		"n_5_plus",
		"completeness",
		"contamination",
		"strain_het",
		"superkingdom",
		"kingdom",
		"phylum",
		"class",
		"order",
		"family",
		"genus"]



# iterate over file
with open(checkm_file, 'r') as input, open(outfile, 'w') as out:
    out.write('\t'.join(map(str,titles)) + '\n')
    for row in input:

    	# split on whitespace
    	arr = row.split()

    	# only consider data lines
    	if (len(arr) > 1):

    		# get taxonomy free of the k__ bit
    		if (arr[1]=="root"):
    			tax = "root"
    		else:
    			tax = arr[1].split("__")[1]

    		# map taxid and tax names
    		name2taxid = ncbi.get_name_translator([tax])

    		# empty variables unless we change them
            lineages = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"]
    		values = {x : "" for x in lineages}

    		# check we got what we asked for
    		if tax in name2taxid.keys():

    			# we want the taxonomy ID
    			taxid = name2taxid[tax]

    			# get entire lineage from this tax id
    			lineage = ncbi.get_lineage(taxid[0])

    			# get all names for that lineage
    			names = ncbi.get_taxid_translator(lineage)

    			# iterate up the lineage mapping names
    			# to each of our variables
    			for l in lineage:
    				rank = ncbi.get_rank([l])
                    if rank[l] in values:
                        values[rank[l]] = names[l]

            ordered = [values[x] for x in lineages]
    		# print it all out
    		out.write('\t'.join(map(str,arr)) + '\t')
    		out.write("\t".join(ordered) + '\n')
