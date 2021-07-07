#!/usr/bin/env python3
# This is a script designed to create a smaller subtable for blobtools plotting

import sys

usage = f'{sys.argv[0]} <input file> <output file>\n'
if len(sys.argv) < 3:
    print("Error in input! Please enter required arguments!")
    print(usage)
    sys.exit(-1)

with open(sys.argv[1], 'r') as input, open(sys.argv[2], 'w') as output:
    output.write("LEN\tGC\tKING\n")
    supkingcol = 17
    for l in input:
        if l.startswith('#'):
            s = l.rstrip().split()
            if len(s) < 2:
                continue
            if s[1] == "name":
                for i in range(len(s)):
                    if s[i].startswith("superkingdom"):
                        supkingcol = i - 1
                        break
            continue
        s = l.rstrip().split("\t")
        output.write(f'{s[1]}\t{s[2]}\t{s[supkingcol]}\n')

print(f'Succsessfully written to {sys.argv[2]}!')
