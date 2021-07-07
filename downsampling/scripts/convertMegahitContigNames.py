# Megahit doesn't follow the fasta format so we must do this!
import sys
import re

usage = "python3 script.py <input megahit fasta> <output reformated fasta>"

if len(sys.argv) != 3:
    print(usage)
    sys.exit(-1)

with open(sys.argv[1], 'r') as input, open(sys.argv[2], 'w') as out:
    seq = ''
    previous = ''
    for l in input:
        s = l.rstrip().split()
        if s[0].startswith('>'):
            if len(seq) > 1:
                start = 0
                llen = len(seq)
                out.write(previous + '\n')
                while llen - start > 60:
                    out.write(seq[start:start+60] + '\n')
                    start += 60
                out.write(seq[start:] + '\n')
            previous = s[0]
            seq = ''
        else:
            seq += s[0]

print("Finished with reformating fasta!")
