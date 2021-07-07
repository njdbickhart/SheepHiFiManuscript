import os
from collections import defaultdict
import subprocess

os.makedirs(snakemake.ouput["dir"])
bins = defaultdict(list)
with open(snakemake.input["cluster"], 'r') as clust:
    for l in clust:
        s = l.rstrip().split()
        bins[s[1]].append(s[0])
with open(snakemake.output["counts"], 'w') as out:
    for bins, ctgs in bins.items():
        outfile = snakemake.output["dir"] + f'/{bins}.fa'
        cmd = f'samtools faidx {snakemake.params["reference"]} ' + ' '.join(ctgs) + f' > {outfile}'
        print(cmd)
        subprocess.Popen(cmd, shell=True)
        out.write(f'{bins}\t{len(ctgs)}\n')
