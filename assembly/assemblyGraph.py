import pandas as pd
import numpy as np 
import subprocess
import os
import bisect

class contig():

    def __init__(self, L):
        self.L = L
        self.breakpoints = [0,L]
        self.p_terminal = None
        self.m_terminal = None 
        
    def __str__(self):
        return str(self.breakpoints)

    def add_segment(self,interval):
        if (interval[0] not in self.breakpoints):
            bisect.insort(self.breakpoints,interval[0])
        if (interval[1] not in self.breakpoints):
            bisect.insort(self.breakpoints,interval[1])

class merger():
    
    def __init__(self, gfa_ref, gfa_query, fasta_ref, fasta_query,aln_loc):

        # Store file names
        self.gfa_r = gfa_ref
        self.gfa_q = gfa_query
        self.fa_r = fasta_ref
        self.fa_q = fasta_query
        self.aln = aln_loc

        # Read in GFA files to build initial assembly graph structures.
        self.rcontig,self.rID = self._readGFA(self.gfa_r)
        self.qcontig,self.qID = self._readGFA(self.gfa_q)

        # Align the assemblies together.
        self._minimap2()
        print self.rcontig[0]
        print self.qcontig[0]
        self._parseMappingSegmentation()
        print self.rcontig[0]
        print self.qcontig[0]
        self._registerMapping()

    def _readGFA(self,gfa):
        contig_id = {}
        n = 0
        contigs = []
        with open(gfa,'r') as GFA:
            for line in GFA:
                words = line.split()
                if (words[0] == "S"):
                    contig_id[words[1]] = n
                    contigs.append(contig( int( words[3].split(":")[2]) ) )
                    n += 1
                elif (words[0] == "L"):
                    contig_1 = contig_id[words[1]]
                    contig_2 = contig_id[words[3]]
                    
                    if (words[2]=="+" and words[4] == "+"):
                        contigs[contig_1].p_terminal = contig_2
                        contigs[contig_2].p_terminal = contig_1
                    elif (words[2]=="+" and words[4] == "-"):
                        contigs[contig_1].p_terminal = contig_2
                        contigs[contig_2].m_terminal = contig_1
                    elif (words[2]=="-" and words[4] == "+"):
                        contigs[contig_1].m_terminal = contig_2
                        contigs[contig_2].p_terminal = contig_1
                    else:
                        contigs[contig_1].m_terminal = contig_2
                        contigs[contig_2].m_terminal = contig_1

        return contigs,contig_id

    def _minimap2(self):
        if (not os.path.exists(self.aln)):
            map_cmd = "minimap2 -x asm5 " + self.fa_r + " " + self.fa_q + " > " + self.aln
            subprocess.call(map_cmd,shell=True)

    def _parseMappingSegmentation(self):
        aln = pd.read_csv(self.aln, sep="\t", header=None)
        for row in aln.iterrows():
            if (row[1][10] > 1000):
                qContig = self.qID[row[1][0].rstrip("_pilon")]
                qInterval = (row[1][2],row[1][2])
                rContig = self.rID[str(row[1][5]).rstrip("_pilon")]
                rInterval = (row[1][7],row[1][8])
                
                self.qcontig[qContig].add_segment(qInterval)
                self.rcontig[rContig].add_segment(rInterval)



if __name__ == "__main__":
    carb_folder = "/scicore/home/neher/GROUP/data/Carbapenemases/carb045/"
    gfa1 = carb_folder + "unicycler_hybrid_nano/assembly.gfa"
    gfa2 = carb_folder + "canu_nano/assembly.contigs.gfa"
    fa_1 = carb_folder + "unicycler_hybrid_nano/assembly.fasta"
    fa_2 = carb_folder + "canu_nano/assembly.contigs.fasta"
    aln_loc = carb_folder + "unicycler_hybrid_nano/aln.paf"

    x = merger(gfa1, gfa2, fa_1, fa_2, aln_loc)