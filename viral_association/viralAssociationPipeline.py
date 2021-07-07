# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 15:19:15 2019

@author: Derek M Bickhart
Requirements: Minimap2 and Samtools on your PATH (or entered as arguments at the parser)
Python Version 3.6+
"""

import argparse
import re
import subprocess
import sys
from collections import defaultdict
import numpy as np
#import scipy.cluster.vq as sp
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rundef = re.compile(r'-undef')

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Process long-read alignments and Hi-C links to identify likely viral-host associations in metagenomic assembly data"
            )
    parser.add_argument('-l', '--long_read', 
                        help="Input long-read fastq file",
                        type=str, default = "None"
                        )
    parser.add_argument('-a', '--assembly', 
                        help="Fasta of the assembled contigs (in a single file) for analysis",
                        type=str, required=True
                        )
    parser.add_argument('-g', '--viral_contigs',
                        help="Fasta file of the separated contigs of viral origin (a subset of the full assembly)",
                        type=str, default = "None"
                        )
    parser.add_argument('-b', '--blob_tools',
                        help="Input blob tools or taxonomic data table",
                        type=str, required=True
                        )
    parser.add_argument('-i', '--hic_links',
                        help="Input sam/bam file with alignments of Hi-C reads",
                        type=str, default = "None"
                        )
    parser.add_argument('-v', '--viruses',
                        help="Tab delimited list of contigs containing viral sequence and their lengths",
                        type=str, required=True
                        )
    parser.add_argument('-c', '--link_thresh',
                        help="Filter for the number of stdevs of Hi-C links above the average count to be used in viral Hi-C association [2.5]",
                        type=float, default=2.5
                        )
    parser.add_argument('-e', '--overhang',
                        help="Filter for long-read overhang off of viral contigs [150]",
                        type=int, default=150
                        )
    parser.add_argument('-n', '--noplot',
                        help="[optional flag] Disables plotting of data",
                        action='store_true'
                        )
    parser.add_argument('-o', '--output',
                        help="Output basename",
                        type=str, required=True
                        )
    parser.add_argument('-m', '--minimap',
                        help="Path to the minimap executable",
                        type=str, default = "minimap"
                        )
    parser.add_argument('-s', '--samtools',
                        help="Path to the samtools executable",
                        type=str, default = "samtools"
                        )
    
    return parser.parse_args()

def generateViralSubset(samtools : str, assembly : str, vctgs : str, outfile = "temp_viral.fa"):
    vlist = list()
    with open(vctgs, 'r') as input:
        for l in input:
            l = l.rstrip()
            s = l.split()
            vlist.append(s[0])
            
    cmd = [samtools, 'faidx', assembly] + vlist
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, encoding='utf8') 
    with open(outfile, 'w') as out:
        while True:
            l = proc.stdout.readline()
            if l == '' and proc.poll is not None:
                break
            out.write(l)
    print("Created viral subset fasta: " + outfile)

def main(args):
    # Boolean counters to confirm that at least one data type is entered!
    reads = True if args.long_read != "None" else False
    hic = True if args.hic_links != "None" else False
    
    if not reads and not hic:
        print("ERROR: You must enter at least one data type! Either long-reads or hi-C links!")
        print(args.print_help())
        sys.exit(-1)
        
    # Create workhorse class and begin processing the file
    workhorse = viralComparison(args.viruses)
    
    # Now the script bifurcates. Process each data type in order
    if reads:
        vctgfile = args.viral_contigs
        if vctgfile == "None":
            print("Generating a subset list of viral contigs in file: temp_viral.fa")
            generateViralSubset(args.samtools, args.assembly, args.viruses, outfile = "temp_viral.fa")
            vctgfile = "temp_viral.fa"
            
        # Long read alignment to viruses
        workhorse.alignECReads(vctgfile, args.long_read, args.minimap, 
                               args.output + '.algn.viruses', oThresh = args.overhang)
        
        # Align overhangs back to assembly
        workhorse.realignECOverhangs(args.assembly, args.long_read, args.minimap, 
                                     args.output + '.ovlps.fa', args.output + '.lread.vir.graph', 
                                     args.samtools)
        
        print("Finished long read overlap!")
    if hic:
        workhorse.generateHiCLinkTable(args.samtools, args.hic_links, 
                                       args.output + '.hiclinks.tab', 
                                       args.link_thresh)
        print("Finished Hi-C link processing!")
        
    # Merge the results
    workhorse.combineTables(reads, hic)
    
    # Apply taxonomic data to the merged table
    workhorse.loadTaxonomy(args.blob_tools)
    
    # Print out the final table
    workhorse.printOutFinalTable(args.output + '.final.tab')
    
    # By default, try to plot in a network
    if not args.noplot:
        workhorse.generatePlot(args.output + '.nxplot.png')
    
class viralComparison:

    def __init__(self, viralCtgFile : str):
        self.viruses = dict()
        
        # Load contigs and lengths
        with open(viralCtgFile, 'r') as fh:
            for l in fh:
                l = l.rstrip()
                s = l.split()
                self.viruses[s[0]] = int(s[1])
        
        self.ovlpSizes = list()
        self.ovlpEC = defaultdict(list) # {readname} -> [start, end, vctg]
        self.ovlpASM = defaultdict(dict)
        self.hicASM = defaultdict(dict)        
        
        self.totalCount = 0
        self.finalTable = defaultdict(dict) #{virus} -> {host} = VAssoc
        
        self.readErrors = 0
        
        print(f'Loaded {len(self.viruses)} viral contigs for analysis\n')
        
    def isVirus(self, ctg : str):
        return ctg in self.viruses
    
    def generatePlot(self, outfile : str):
        vcolors = ['b', 'g', 'r', 'c', 'y']
        hcolors = ['m', 'w', 'k']
        ecolors = {"Read" : 'b', "HiC" : 'y', "Both" : 'g'}
        graph = nx.Graph()
        vgenus = defaultdict(list)
        hking = defaultdict(list)
        ecat = defaultdict(list)
        for v in self.finalTable:
                for h in self.finalTable[v]:
                    working = self.finalTable[v][h].getListOutput()
                    graph.add_node(working[0])
                    graph.add_node(working[1])
                    vgenus[working[3]].append(working[0])
                    hking[working[4]].append(working[1])
                    graph.add_edge(working[0], working[1])
                    ecat[working[2]].append((working[0], working[1]))
        
        pos = nx.spring_layout(graph)
        # Set viral node styles
        for i in range(len(vgenus.keys())):
            k = list(vgenus.keys())[i]
            if i >= 5:
                i = 4
            print(f'Viral Genus {k} color: {vcolors[i]}')
            nx.draw_networkx_nodes(graph, pos, nodelist=vgenus[k], node_color=vcolors[i], alpha=0.8, node_size=500, node_shape='8')
        
        # Set host node styles
        for i in range(len(hking.keys())):
            k = list(hking.keys())[i]
            if i >= 3:
                i = 2
            print(f'Host kingdom {k} color: {hcolors[i]}')
            nx.draw_networkx_nodes(graph, pos, nodelist=hking[k], node_color=hcolors[i], alpha=0.8, node_size=400, node_shape='o')
            
        # Set edge styles
        for k, v in ecolors.items():
            if not k in ecat:
                print(f'Houston, we have a problem with cat: {k}')
                continue
            nx.draw_networkx_edges(graph, pos, edgelist=ecat[k], width=2, alpha=0.6, edge_color=ecolors[k])
            
        plt.axis('off')
        plt.savefig(outfile)

    
    def printOutFinalTable(self, outtab : str):
        print("Final associations in table {}".format(outtab))
        with open(outtab, 'w') as out:
            out.write("VirusCtg\tHostCtg\tCategory\tVirusGenus\tHostKingdom\tHostGenus\tEvidence\n")
            for v in sorted(self.finalTable):
                for h in sorted(self.finalTable[v]):
                    out.write('\t'.join(self.finalTable[v][h].getListOutput()) + '\n')
    
    def loadTaxonomy(self, blobtable : str):
        print(f'Loading taxonomic information from file {blobtable}')
        # first, get the set of contig names that we need to process
        contigs = set()
        for v, j in self.finalTable.items():
            contigs.add(v)
            for h, c in j.items():
                contigs.add(h)
              
        taxonomy = defaultdict(list) #{contig} -> [kingdom, genus]
        with open(blobtable, 'r') as input:
            # read until first line and get index information
            kingidx = 0
            genusidx = 0
            while(True):
                line = input.readline()
                if line.startswith("##"):
                    continue
                elif line.startswith("#"):
                    line = line.rstrip()
                    segs = line.split(sep="\t")
                    for i in range(len(segs)):
                        if segs[i].startswith("genus.t"):
                            genusidx = i
                        if segs[i].startswith("superkingdom.t"):
                            kingidx = i
                    break
            if kingidx == 0 or genusidx == 0:
                print("Could not ID taxonomic table entries! Did you remember to add genus and superkingdom table headers?")
                return
            
            for l in input.readlines():
                l = l.rstrip()
                segs = l.split(sep="\t")
                if segs[0] in contigs:
                    segs[kingidx] = rundef.sub('', segs[kingidx])
                    segs[genusidx] = rundef.sub('', segs[genusidx])
                    taxonomy[segs[0]] = [segs[kingidx], segs[genusidx]]
        
        # Now assign tax affiliations to the contigs
        taxassign = 0
        for v, j in self.finalTable.items():
            for h, c in j.items():
                if v in taxonomy and h in taxonomy:
                    taxassign += 1
                    c.setTaxonomy(taxonomy[v][1], taxonomy[h][0], taxonomy[h][1])
        print(f'Filled in taxonomic entries for {taxassign} associations out of {self.totalCount} possible entries')
                
    def combineTables(self, reads : bool, hic : bool):
        print("Generating final network table")
        if reads:
            for v, j in self.ovlpASM.items():
                for h, c in j.items():
                    self.totalCount += 1
                    self.finalTable[v][h] = c
           
        overlaps = 0
        if hic:
            for v, j in self.hicASM.items():
                for h, c in j.items():
                    if v in self.finalTable:
                        if h in self.finalTable[v]:
                            overlaps += 1
                            self.finalTable[v][h].combine(c)
                            continue
                    self.totalCount += 1
                    self.finalTable[v][h] = c
        
        if reads and hic:
            print("Identified {} overlapping viral-host associations.".format(overlaps))
    
    def generateHiCLinkTable(self, samtools : str, insam : str, outtab : str, hiThresh : float):
        print("Generating a Hi-C link table from alignment file: {}".format(insam))
        hicLinks = defaultdict(lambda : defaultdict(int))
        selfLinks = defaultdict(float)
        cmd = [samtools, 'view', insam]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, encoding='utf8')
        
        while True:
            l = proc.stdout.readline()
            if l == '' and proc.poll() is not None:
                break
            if not l.startswith('@'):
                segs = l.split()
                if len(segs) < 7:
                    print(f'Error parsing malformed sam at record: {" ".join(segs)}')
                    continue
                if segs[6] == segs[2] or segs[6] == "=":
                    selfLinks[segs[2]] += 0.5   # half a count to avoid double-counting pairs
                if segs[6] != segs[2] and segs[6] != "=" and self.isVirus(segs[2]):
                    hicLinks[segs[2]][segs[6]] += 1
        
        # Let's try k-means clustering to divide the samples
        # assuming that viral intercontig link noise will be in cluster 1 and signal in cluster 2
        vLCounts = list()
        for v, j in hicLinks.items():
            for h, n in j.items():
                vLCounts.append(n)
             
        # What follows is an abandoned attempt to use k-means clustering to determine optimal hi-C threshold signals
        #flattened = np.array(vLCounts)
        #centroids,_ = sp.kmeans(flattened, 2)
        #idx,_ = sp.vq(flattened, centroids)
        
        #minV1 = min([int(x) for x in flattened[idx==0]])
        #minV2 = min([int(x) for x in flattened[idx==1]])
        
        #kThresh = minV1 if minV1 > minV2 else minV2
        #SNR = np.where(flattened.std(axis=0, ddof=0) == 0, 0, 
        #               flattened.mean(0) / flattened.std(axis=0, ddof=0))
        
        #print(f'K-means clustering identified a minimum intercontig link count of {kThresh} from {len(vLCounts)} observations with a signal-to-noise ratio of {SNR[0]}')
        # More simplistic, but hopefully representative: use number of reads above stdevs from the mean
        vLAvg = np.mean(vLCounts)
        vLStd = np.std(vLCounts)
        kThresh = vLAvg + (vLStd * hiThresh)
        
        print(f'Using a threshold of {kThresh} reads to select Hi-C links to host contigs')
        # Now generate intermediary output tab file and fill data structure
        vassocNum = list() # list of number of host contigs associated with viruses
        with open(outtab, 'w') as out:
            for v, j in hicLinks.items():
                vassocNum.append(len(j))
                for h, n in j.items():
                    if n >= kThresh:
                        out.write(f'{v}\t{h}\t{n}\n')
                        self.hicASM[v][h] = VAssoc(h, v, "HiC", n)
                    
        meanHcontig = np.mean(vassocNum)
        print(f'Found valid Hi-C link associations for {len(hicLinks)} viral contigs out of {len(vLCounts)} original candidates.')
        print(f'There were an average of {meanHcontig} host contig associations in this dataset')
    
    def samfaidx(self, samtools : str, reffasta : str, ctglst, outfasta):
        cmd = [samtools, 'faidx', reffasta]
        
        proc = subprocess.Popen(cmd + ctglst, stdout=subprocess.PIPE, encoding='utf8')
        while True:
            l = proc.stdout.readline()            
            if l == '' and proc.poll is not None:
                break
            outfasta.write(l)
    
    def realignECOverhangs(self, asmCtgFasta : str, ecReads : str, minimap : str,
                           outfasta: str, outfile : str, samtools : str, 
                           rcountThresh = 3, minimapOpts =['-x', 'map-pb']):
        print("Generating overhangs in {} fasta file".format(outfasta));
        # Create the overhanging read fasta for alignment
        with open(outfasta, 'w') as fasta:
            container = list()
            for k, f in self.ovlpEC.items():
                if f[1] - f[0] < 2:
                    continue
                container.append(f'{k}:{f[0]}-{f[1]}')
                if len(container) >= 100:
                    self.samfaidx(samtools, ecReads, container, fasta)
                    container = list()
            
            if len(container) > 0:
                self.samfaidx(samtools, ecReads, container, fasta)
                container = list()
        
        # Now, map the reads and filter the results
        print("Aligning overhangs to full genome fasta file")
        redundancies = 0
        overlaps = dict() # temp container for EC read associations
        
        proc = subprocess.Popen([minimap] +  minimapOpts + [asmCtgFasta, outfasta], stdout=subprocess.PIPE, encoding='utf8')
        with open(outfile, 'w') as out:
            while True:
                l = proc.stdout.readline()
                if l == '' and proc.poll() is not None:
                    break
                
                l = l.rstrip()
                segs = l.split()
                
                # Get original alignments
                read = segs[0].split(':')
                if not read[0] in self.ovlpEC:
                    print(f'Could not find read: {read[0]} in previous overlaps! Skipping...')
                    continue

                vir = self.ovlpEC[read[0]]
                if read[0] not in overlaps and not self.isVirus(segs[5]):
                    # Note: this preferentially prints out the first alignment encountered
                    if len(vir) < 3:
                        print('\t'.join(vir + [' ; '] + segs))
                    overlaps[read[0]] = VAssoc(segs[5], vir[2])
                    out.write(f'{read[0]}\t{segs[5]}\t{vir[2]}\n')
                elif self.isVirus(segs[5]):
                    continue
                else:
                    overlaps[read[0]].setRedundant()
                    # Ignore subsequent alignments
                    redundancies += 1
        
        successes = len({s.vctg for k, s in overlaps.items()})                    
        print(f'Found successful associations for {successes} out of {len(self.viruses)} viral contigs and {redundancies} ambiguously aligned reads')
        
        # Loading into the final container
        readcounts = defaultdict(lambda : defaultdict(int))
        for k, s in overlaps.items():
            # Skip over any reads that may have been ambiguous alignments 
            if not s.redund:
                readcounts[s.vctg][s.hostctg] += 1
                
        for v, j in readcounts.items():
            for h, c in j.items():
                if c >= rcountThresh:
                    # The above selects only associations that have at least X read alignments
                    self.ovlpASM[v][h] = VAssoc(h, v, "Read", c)
    
    def alignECReads(self, viralCtgFasta : str, ecReads : str, minimap : str, outfile : str,
                     minimapOpts = ['-x', 'map-pb'], algLen = 500, oThresh = 150):
        cmd = [minimap]
        cmd.extend(minimapOpts)
        cmd.extend([viralCtgFasta, ecReads])
        
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, encoding='utf8')
        with open(outfile, 'w') as out:
            while True:
                l = proc.stdout.readline()
                if l == '' and proc.poll() is not None:
                    break # Terminate at the end of the process
                
                l = l.rstrip()
                segs = l.split()
                
                if int(segs[9]) > algLen and (int(segs[7]) < oThresh or int(segs[6]) - int(segs[8]) < oThresh):
                    self.getOverhang(segs[0], segs[5], oThresh, int(segs[7]), 
                                     int(segs[8]), int(segs[2]), int(segs[3]), 
                                     int(segs[1]), segs[4], out)
        count = len(self.ovlpSizes)
        mean = np.mean(self.ovlpSizes)
        viruscounts = len({s[2] for s in self.ovlpEC})
        print(f'Identified {count} overlapping reads with a mean length of {mean} for {viruscounts} unique viral contigs')
    
    def getOverhang(self, rname : str, vctg : str, othresh : int, vstart : int, 
                    vend : int, rstart : int, rend : int, rlen: int, rorient : str,
                    ofh):
        right_end = rlen - rend < othresh
        read_unmapright = self.viruses[vctg] - vend > self.viruses[vctg] - vstart
        
        if rorient == "+":
            # read: ------
            # ctg:  --
            if right_end and read_unmapright:
                ofh.write(f'{rname}\t{rend}\t{rlen}\t{vctg}\n')
                self.ovlpSizes.append(rlen - rend)
                self.ovlpEC[rname] = [rend, rlen, vctg]
            elif not right_end and not read_unmapright:
                # read: ------
                # ctg:      --
                ofh.write(f'{rname}\t0\t{rstart}\t{vctg}\n')
                self.ovlpSizes.append(rstart)
                self.ovlpEC[rname] = [1, rstart, vctg]
            else:
                self.readErrors += 1
        else:
            # read: -----
            # ctg:  --
            if right_end and not read_unmapright:
                ofh.write(f'{rname}\t0\t{rstart}\t{vctg}\n')
                self.ovlpSizes.append(rstart)
                self.ovlpEC[rname] = [1, rstart, vctg]
            elif not right_end and not read_unmapright:
                # read: ------
                # ctg:     ---
                ofh.write(f'{rname}\t{rend}\t{rlen}\t{vctg}\n')
                self.ovlpSizes.append(rlen - rend)
                self.ovlpEC[rname] = [rend, rlen, vctg]
            else:
                self.readErrors += 1

class VAssoc:
    def __init__(self, hostctg : str, vctg : str, cat = "Read", count = 0):
        self.hostctg = hostctg
        self.vctg = vctg
        self.redund = False
        self.category = cat
        self.count = count
        # Taxonomic placeholders
        self.vTax = "N/A"
        self.hostKing = "N/A"
        self.hostGenus = "N/A"
        
        # used only if combined
        self.complex = defaultdict(int)
    
    def setRedundant(self):
        self.redund = True      
        
    def setTaxonomy(self, vTax : str, hKing : str, hGen : str):
        self.vTax = vTax
        self.hostKing = hKing
        self.hostGenus = hGen
    
    def combine(self, other : "VAssoc"):
        self.complex[self.category] += self.count
        self.complex[other.category] += other.count
        self.category = "Both"
        
    def getEvidence(self) -> str:
        if self.category == "Both":
            return ';'.join([f'{k}:{v}' for k, v in self.complex.items()])
        else:
            return f'{self.category}:{self.count}'
        
    def getListOutput(self) -> list:
        return [self.vctg, self.hostctg, self.category, self.vTax, self.hostKing, 
                self.hostGenus, self.getEvidence()]
  
if __name__ == "__main__":
    args = parse_user_input()
    main(args)
