#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 10:40:56 2017

@author: nolln
"""

import numpy as np    
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI
from Bio import SeqIO
import gzip
from itertools import izip
import glob

def returnPlasmids(gbk,name,cardDB):

    for n,contig in enumerate(SeqIO.parse(gbk,'genbank')):
        L = len(contig.seq)
        for feature in contig.features:
            if ('note' in feature.qualifiers and 'ARO' in feature.qualifiers['note']):
                print feature.qualifiers['note']

def resLocations(gbk):
    
    resLocs = []
    contigL = []
    for n,contig in enumerate(SeqIO.parse(gbk,'genbank')):
        contigL.append(len(contig.seq))
        resLocs_contig = []
        for feature in contig.features:
            if ('note' in feature.qualifiers and 'ARO' in feature.qualifiers['note'][0]):
                # Get gene location.
                start = []
                end = []
                for interval in feature.location.parts: 
                    start.append(interval.nofuzzy_start)
                    end.append(interval.nofuzzy_end)
                if (len(start)==1):
                    resLocs_contig.append( (start[0],end[0],feature.qualifiers['note']) )
                else:
                    resLocs_contig.append( (start[1],end[0],feature.qualifiers['note']) )
        resLocs.append(resLocs_contig)
        
    return resLocs,contigL

def groupGenesAcrossPG(pgLib,pattern):
    """ 
    Scans pangenomes located within pgLib and returns all genes with names that
    map to `pattern'. It then concatenates individuals across PGs that share
    the same gene.
    """
    from collections import defaultdict
    import json
    import os,shutil
    import subprocess
    
    panGenomes = glob.glob(pgLib + '/*') # Obtain list of pangenomes to scan.
    genes = defaultdict(list)
    for pG in panGenomes:
        pG_name = pG[len(pgLib):]
        if (os.path.exists(pG + '/vis/geneCluster.json')):
            with open(pG + '/vis/geneCluster.json') as clusterFile:
                gC = json.load(clusterFile)
                for cluster in gC:
                    if (pattern in cluster['ann']):
                        genes[str(cluster['ann'])].append( [pG_name, str(cluster['msa'])] )
    
    # Build the concatenated alignment files for only genes with cross-species data 
    abFolder = pgLib + 'AB_genes/'
    if (not os.path.exists(abFolder) ):
        os.mkdir(abFolder)  
          
    for name in genes.keys():
        if (len(genes[name]) > 1):
            outName = abFolder + name + '.fna'
            with open(outName, 'w') as outFaa:
                records = []
                for gC in genes[name]:
                    gC_folder = pgLib + gC[0] + '/geneCluster/'
                    for seq in SeqIO.parse(gC_folder + gC[1] + '.fna','fasta'):
                        seq.id = seq.id.replace('isolate',gC[0])
                        seq.name = seq.name.replace('isolate',gC[0])
                        seq.description = seq.description.replace('isolate',gC[0])
                        records.append(seq)
                SeqIO.write(records,outFaa,'fasta')
            # Sequences are unaligned. Use Mafft to align.
            maft = 'mafft --reorder --nuc '+ outName +' 1> ' + outName[:-4] + 'aln.fna'
            process = subprocess.Popen(maft,shell=True) ;
            process.communicate()
            shutil.move(outName[:-4] + 'aln.fna',outName)
            
def qualAlongRead(library,readLen):
    """
    Returns the distributions of read quality as a function of read length.
    For quality control purposes. 
    """
    
    SANGER_OFFSET = ord('!')
    q_mapping = dict()
    for letter in range(0,255):
        q_mapping[chr(letter)] = letter - SANGER_OFFSET
        
    if (library[0].endswith('.gz')):
        openf = gzip.open
        readMode = 'rb'
    else:
        openf = open
        readMode = 'r'
      
    quality = [ [ [] for n in xrange(readLen)] for PE in xrange(2) ]
    lenOfRead = [ [] for PE in xrange(2) ]
    with openf(library[0],readMode) as lib1, openf(library[1],readMode) as lib2:
        for (n,reads) in enumerate(izip(FGI(lib1),FGI(lib2))):
            for (pE,read) in enumerate(reads):
                lenOfRead[pE].append(len(read[2]))
                for pos,qletter in enumerate(read[2]):
                    if (pos < readLen):
                        quality[pE][pos].append(q_mapping[qletter])
                    else:
                        print pos
                
    quality = [ [np.array(q) for q in quality[i]] for i in xrange(2)]
    lenOfRead = np.array(lenOfRead)
    return quality,lenOfRead

def readLengths(library):
    """
    Returns the distributions of read quality as a function of read length.
    For quality control purposes. 
    """
        
    if (library.endswith('.gz')):
        openf = gzip.open
        readMode = 'rb'
    else:
        openf = open
        readMode = 'r'
      
    readLen = []
    with openf(library,readMode) as lib1:
        for (n,reads) in enumerate(FGI(lib1)):
            readLen.append(len(reads[1]))
                
    readLen = np.array(readLen)
    return readLen

def plotDistributions(data,title,lineC,f=None,ax=None):
    """
    Plots the distribution of reads lengths (or any other data defined over 
    reads) over both paired end sets
    """ 
    
    from statsmodels.distributions import ECDF
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pylab as plt
    
    if (f is None):
        f, ax = plt.subplots(2,sharey=True,sharex=True)
    
    ecdfPE1 = ECDF(data[0,:])
    ecdfPE2 = ECDF(data[1,:])
    ax[0].step(ecdfPE1.x,ecdfPE1.y,color=lineC,alpha=.5)
    ax[1].step(ecdfPE2.x,ecdfPE2.y,color=lineC,alpha=.5)
    ax[0].set_title('Read 0')
    ax[1].set_title('Read 1')

    f.suptitle('Cumulative Distributions of Read Length')
    return f,ax

def plotQualAlongRead(qual,color='b',f=None,ax=None):
    """
    Plots the distributions of read quality as a function of read length.
    For quality control purposes. 
    """
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pylab as plt
    
    xLocs = np.array(range(len(qual[0])))
    
    if (f is None):
        f, ax = plt.subplots(2,sharey=True,sharex=True)

    for pe in xrange(2):
        quartile1 = np.array( [ np.percentile(1.*q,25) for q in qual[pe] ] )
        medianQ = np.array([ np.percentile(1.*q,50) for q in qual[pe] ] )
        quartile2 = np.array([ np.percentile(1.*q,75) for q in qual[pe] ] )
    
        ax[pe].vlines(xLocs,quartile1-.1,quartile2+.1,color=color,linestyle='-',lw=2,alpha=.25)
        ax[pe].scatter(xLocs,medianQ,marker='o',color=color,edgecolor='k',s=50)
        ax[pe].set_title('Read ' + str(pe))
    plt.suptitle('Phred Quality score vs read length')
    
    return f,ax
    
def getBUSCOStats(folder):
    buscoFile = open(folder + 'short_summary_busco.txt','r')
    bScore = '00.0%'
    for line in buscoFile.readlines():
        if ('C:' in line):
            bScore = line.split('\t')[1].split(',')[0][2:7]
    return bScore

def getAssemblyStats(assemblyFile,alpha):
    
    if (assemblyFile.endswith('.gz')):
        openf = gzip.open
        readMode = 'rb'
    else:
        openf = open
        readMode = 'r'
        
    contigLens = []
    with openf(assemblyFile,readMode) as contigFile:
        for contig in SeqIO.parse(contigFile,'fasta'):
            contigLens.append(len(contig.seq))
    contigLens = np.array(contigLens)
    contigLens = np.sort(contigLens)[::-1]
    prcAssemblyCovered = 1. * np.cumsum(contigLens) / np.sum(contigLens)
    Lalpha = np.where(prcAssemblyCovered > alpha)[0][0]
    Nalpha = contigLens[Lalpha]
    
    return Nalpha,Lalpha,contigLens[0],np.sum(contigLens)

def obtainReadNames(directory,fq=False):
    if (not directory.endswith('/')):
        directory += '/'
    if (not fq):
        PE_1 = glob.glob(directory + '*R1*' + '.fastq.gz' )
        PE_2 = glob.glob(directory + '*R2*' + '.fastq.gz' )
    else:
        PE_1 = glob.glob(directory + '*R1*' + '.fq.gz' )
        PE_2 = glob.glob(directory + '*R2*' + '.fq.gz' )
        
    match2 = [fileName.replace('R1','R2').replace('_val_1','_val_2') for fileName in PE_1]
    matchIdx = [PE_2.index(match) for match in match2]
   
    return [ (PE_1[i],PE_2[m]) for i,m in enumerate(matchIdx)]

def loadPileUp( dirName ):
    import os
    
    refNames = []
    pileUp = []
    for fileName in os.listdir(dirName):
        if ( fileName.endswith('.npz') or fileName.endswith('.npy') ):
            refNames.append(fileName.split('_pileUp.npz')[0])
            pileUp.append(np.load(dirName + fileName)['arr_0'])
            
    return refNames,pileUp
   
def computeSiteEntropy( pileUp ):
    
    S = []
    for n in xrange(len(pileUp)):
        coverage = 1.0 * np.sum(pileUp[n],axis = 0,keepdims=True).T
        coverage[np.where(coverage==0)] = 1
        nucFreq = pileUp[n].T / coverage
        logArg = nucFreq
        logArg[np.where(logArg==0)] = 1
        
        entropy = nucFreq * np.log(logArg)
        S.append( -np.sum( entropy, axis=1) )
    return S

def computeCoverage( pileUp ):
    
    cov = []
    for n in xrange(len(pileUp)):
        cov.append(np.sum(pileUp[n],axis=0))
        
    return cov

def computeMedianCoverage( pileUp ):
    
    cov = []
    for n in xrange(len(pileUp)):
        cov.extend(np.sum(pileUp[n],axis=0).tolist())
        
    return np.median(np.array(cov))

def plotCDF( data, logScale=False ):
    
    from statsmodels.distributions.empirical_distribution import ECDF
    import matplotlib.pylab as plt
    
    ecdf = ECDF(data)
    plt.figure()
    plt.plot(ecdf.x,ecdf.y,linewidth=2)
    plt.ylabel('CDF')
    if logScale:
        plt.xscale('log')
    plt.show()
    
def plotSpatial( data, refNames ):
    import matplotlib.pylab as plt

    plt.figure(figsize = (18,8))
    offset = 0
    ticks = []
    plt.title('spatial')
    for n in xrange(len(refNames)):
        plt.plot(offset+np.arange(refNames.shape[-1]),data[n])
        ticks.append([])
    