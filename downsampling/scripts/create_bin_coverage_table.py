import os
import sys
import numpy as np
from collections import defaultdict

usage = f'{sys.argv[0]} <assembly> <input bin scaffold file> <input blobtools cov file> <input dastool file> <output file>\n'


class bin:

    def __init__(self, name, asm):
        self.name = name
        self.asm = asm
        self.comp = 0.0
        self.cont = 0.0
        self.cov = 0.0

    def getPrintOut(self):
        return f'{self.name}\t{self.cov}\t{self.comp}\t{self.cont}\t{self.asm}\n'

if len(sys.argv) < 6:
    print(usage)
    sys.exit()

# Create contig to bin tool and bin dictionary
ctgToBin = dict()
binDict = dict()
with open(sys.argv[2], 'r') as binfile:
    for l in binfile:
        s = l.rstrip().split()
        ctgToBin[s[0]] = s[1]
        if s[1] not in binDict:
            binDict[s[1]] = bin(s[1], sys.argv[1])

# Run through blobtools file
binctgcovs = defaultdict(list)
with open(sys.argv[3], 'r') as blob:
    covidx = 0
    for l in blob:
        if l.startswith('##'):
            continue
        if l.startswith('#'):
            s = l.rstrip().split('\t')
            for i, x in enumerate(s):
                if x == "cov_sum" or x == "cov0":
                    covidx = i
                    print(f'Covidx = {covidx}')
        s = l.rstrip().split('\t')
        if s[0] in ctgToBin:
            binctgcovs[ctgToBin[s[0]]].append(float(s[covidx]))

print(f'Identified {len(binctgcovs)} bins with values out of {len(binDict)} bins')

for b, cs in binctgcovs.items():
    avg = np.mean(cs)
    binDict[b].cov = avg

# Run through dastool file
with open(sys.argv[4], 'r') as das:
    das.readline()
    for l in das:
        s = l.rstrip().split()
        if s[0] in binDict:
            binDict[s[0]].comp = float(s[-2])
            binDict[s[0]].cont = float(s[-1])

# Write it all out
with open(sys.argv[5], 'w') as out:
    out.write("Contig\tCoverage\tCompleteness\tContamination\tAssembly\n")
    for b, c in binDict.items():
        out.write(c.getPrintOut())
