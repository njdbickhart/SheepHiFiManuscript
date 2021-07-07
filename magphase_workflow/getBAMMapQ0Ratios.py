# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 14:51:53 2021

@author: derek.bickhart-adm
"""

import argparse
import pysam
from collections import defaultdict

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "Window analysis to calculate ratio of mapq0 reads in bam files"
            )
    parser.add_argument('-b', '--bam', 
                        help="Input bam file. Must be coord sorted and indexed.",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file name. Output format is bed format",
                        required=True, type=str,
                        )
    parser.add_argument('-w', '--windowlength',
                        help="Window length for contig windows [5000]",
                        default = 5000, type=int,
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    with pysam.AlignmentFile(args.bam, 'rb') as bamfile:
        # create windows
        references = bamfile.references
        lengths = bamfile.lengths 
        winlist = defaultdict(list)
        for r, l in zip(references, lengths):
            for i in range(0, l, args.windowlength):
                winlist[r].append(window(r, i, i + args.windowlength))
        
        
        for c, w in winlist.items():
            for i, win in enumerate(w):
                count = 0
                mapq0 = 0
                for s in bamfile.fetch(c, win.start, win.end):
                    if s.is_secondary:
                        continue
                    count += 1
                    if s.mapping_quality == 0:
                        mapq0 += 1
                winlist[c][i].addCount(mapq0, count)
                
    # Now print it all out
    with open(args.output, 'w') as out:
        for c, w in winlist.items():
            for win in w:
                out.write(win.getBed())

class window:
    
    def __init__(self, contig, start, end):
        self.contig = contig
        self.start = start 
        self.end = end
        self.mapq0 = 0
        self.totReadCount = 0
    
    def addCount(self, mapq0, totReads):
        self.mapq0 = mapq0
        self.totReadCount = totReads
    
    def getRatio(self):
        return self.mapq0 / self.totReadCount
    
    def hasMapq0(self):
        return self.mapq0 != 0
    
    def getBed(self):
        return f'{self.contig}\t{self.start}\t{self.end}\t{self.hasMapq0()}\t{self.getRatio()}\t{self.totReadCount}\t{self.mapq0}\n'

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
