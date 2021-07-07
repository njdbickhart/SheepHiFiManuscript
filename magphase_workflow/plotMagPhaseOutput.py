# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 13:14:28 2020
@author: derek.bickhart-adm
"""

import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from matplotlib.collections import BrokenBarHCollection
from matplotlib import cm
from itertools import cycle
from collections import defaultdict
import argparse
import pandas
import numpy as np
import pysam
import seaborn as sns

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A tool to plot bin and contig level read depth differences in strain assignment"
            )
    parser.add_argument('-f', '--fai', 
                        help="Input reference fasta index file for the bin",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename. Output files are {output}.wins and {output}.pdf",
                        required=True, type=str,
                        )
    parser.add_argument('-b', '--bam', 
                        help="Input CCS read depth bam file",
                        required=True, type=str
                        )
    parser.add_argument('-u', '--human', 
                        help="Input human-readable position variant call file",
                        required=True, type=str
                        )
    parser.add_argument('-i', '--binsize',
                        help="Bin size in bases [5000 bp]",
                        type = int, default=5000
                        )
    parser.add_argument('-e', '--breaks',
                        help="Draw windows in the following bed regions [Optional]",
                        type = str, default="NO"
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    # Get the contig length list
    ctglens = dict()
    with open(args.fai, 'r') as fai:
        for l in fai:
            s = l.rstrip().split()
            ctglens[s[0]] = int(s[1])
            
    # Create windows 
    winlist = defaultdict(list)
    # offset bp to add for stitching contigs together in one line
    ctgoffset = dict()
    breaks = list()
    lastbp = 0
    for c in ctglens:
        ctgoffset[c] = lastbp + 100
        for i in range(0, ctglens[c], args.binsize):
            winlist[c].append(window(c, i, i + args.binsize))
        lastbp += ctglens[c]
        breaks.append(lastbp)
        
    # read each sam region and count the reads
    with pysam.AlignmentFile(args.bam, 'rb') as bamfile:
        for c, w in winlist.items():
            for i, win in enumerate(w):
                count = 0
                for s in bamfile.fetch(c, win.start, win.end):
                    if s.is_secondary:
                        continue
                    count += 1
                winlist = updateWin(winlist, c, i, count)
                
    # Now, read in the human readable text file and process that 
    hapset = set()
    with open(args.human, 'r') as human:
        human.readline()
        for l in human:
            s = l.rstrip().split()
            if '?' in s[0]:
                continue
            # determine where the contig start falls
            for i, win in enumerate(winlist[s[2]]):
                if int(s[3]) < win.end and int(s[3]) >= win.start:
                    if s[5] == 'REF':
                        s[0] = '0'
                    winlist = updateWin(winlist, s[2], i, int(s[6]), s[0])
                    print(f'Updating window: {s[2]} {win.start} {win.end} to {s[6]} for Hap {s[0]}')
                    #hapset.add(s[4])
                    
    # OK, data is in! Let's try plotting
    raw = defaultdict(list)
    bars = list()
    for c, w in winlist.items():
        bars.append([ctgoffset[c], ctglens[c]])
        for win in w:
            for i, h in enumerate(sorted(win.count.keys())):
                raw["contig"].append(c)
                raw["start"].append(win.start + ctgoffset[c])
                raw["end"].append(win.end + ctgoffset[c])
                raw["hap"].append(i)
                raw["count"].append(win.getCount(h))
                
    df = pandas.DataFrame(raw)
    df.to_csv(args.output + '.wins', sep='\t', header=True)
    
    df['pos'] = df["contig"] + "_" + df["start"].astype(str)
    
    
    fig, ax = plt.subplots()
    #ax = df[['pos', 'hap', 'count']].plot.area(x='pos', y='count', hue='hap', ax=ax, colormap='viridis')
    #sns.barplot(data=df[['pos', 'hap','count']], x='pos', y='count', hue='hap', ax=ax, palette='dark')
    tdata = df[['start','hap','count']].pivot_table(index='start', columns='hap', values='count', fill_value=0, aggfunc='sum').unstack().to_frame().rename(columns={0:'Count'})
    print(tdata.head())
    #sns.lineplot(data=df[['start','hap','count']], x='start', y='count', hue='hap', ax=ax, palette='dark')
    sns.lineplot(data=tdata, x='start', y='Count', hue='hap', ax=ax, palette='dark')
    plt.xticks(rotation=45)
    
    if args.breaks != 'NO':
        maxy = tdata['Count'].max()
        with open(args.breaks, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                ctgstart = int(s[1]) + ctgoffset[s[0]]
                ctgend = int(s[2]) + ctgoffset[s[0]]
                midpoint = (ctgend - ctgstart) / 2 + ctgstart
                ax.axvline(ctgstart, ls='-', c='lightgrey', zorder=-1)
                ax.axvline(ctgend, ls='--', c='lightgrey', zorder=-1)
                ax.text(midpoint, maxy, s[3], rotation='vertical', verticalalignment='top', size=4.0)
                print(f'Drawing line for {s[0]} {s[1]} {s[2]} {s[3]} {midpoint} {maxy}')

    if len(breaks) > 1:
        for b in breaks:
            ax.axvline(b, ls='-', c='k', zorder=-1)
    
    #ax.add_collection(BrokenBarHCollection(bars, [-1, 1], facecolors=plt.get_cmap('tab20')))
    #ax.axis('tight')
    plt.savefig(args.output + '.pdf')
    
    

            
            
    
def updateWin(winlist, contig, winidx, count, haplotype = '0'):
    winlist[contig][winidx].addCount(haplotype, count)
    return winlist
    

class window:
    
    def __init__(self, contig, start, end):
        self.contig = contig
        self.start = start 
        self.end = end
        self.haps = dict()
        self.hapcount = 0
        self.count = dict()
    
    def addCount(self, hap, count):
        if hap != '0':
            if hap not in self.haps:
                self.hapcount += 1
                self.haps[hap] = self.hapcount
            hap = str(self.haps[hap])
        self.count[hap] = count
    
    def getCount(self, hap):
        if hap in self.count:
            return self.count[hap]
        else:
            return 0
    

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
