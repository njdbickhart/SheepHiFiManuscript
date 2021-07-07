#!/usr/bin/env python
# Originally from Magpy. Altered to fit this pipeline

import sys

from ete3 import NCBITaxa

# get NCBI taxonomu object
ncbi = NCBITaxa()

if len(sys.argv) == 1:
	print("Please provide a filename")
	sys.exit()


# open the file
checkm_file = snakemake.input[0]
outfile = snakemake.output[0]

# skip three lines
row1 = checkm_file.readline()

# print titles for the output
titles = ["name",
		"nprots",
		"nhits",
		"nfull",
		"genus",
		"ngenus",
		"species",
		"nspecies",
		"avgpid",
		"Superkingdom",
		"kingdom",
		"phylum",
		"class",
		"order",
		"family",
		"genus"]

# iterate over file
with open(checkm_file, 'r') as input, open(outfile, 'w') as out:
    out.write('\t'.join(map(str,titles)) + '\n')
    for row in checkm_file:
    	# split on whitespace
    	arr = row.rstrip('\n\r').split('\t')
    	# only consider data lines
    	if (len(arr) > 1):
    		tax = arr[4]
    		# map taxid and tax name
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
