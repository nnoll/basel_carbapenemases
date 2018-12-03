from subprocess import call
import sys
import glob
import os
import pandas as pd
# import phylopandas as ph
import subprocess
import numpy as np
from collections import defaultdict
from Bio import SeqIO
# import pysam
import toyplot
# from numba import jit
import toytree

class plotter(object):

    def __init__(self,folder='./mnt/neherGroup/Carbapenemases/',longReadTech='nano',assembler='final',shortMapper='bwa'):

        # Intialization of all relevant data parameters!
        self.topFolder = folder.rstrip('/') + "/"
        self.isolateFolders = glob.glob(self.topFolder+'carb???/')
        self.moleculeFolders = glob.glob(self.topFolder+"blaTrees/molecules/*/")
        self.clusterJSON = folder + "blaTrees/geneClusters.json"
        self.assembler = assembler

        if (assembler != 'unicycler' and assembler != 'canu' and assembler != 'spades' and assembler != 'hybrid' and assembler != 'final'):
            raise ValueError("Only recognized assemblers: 'unicycler' or 'canu' or 'spades' or 'hybrid' ")

        self.longReadTech = longReadTech
        if (longReadTech != 'pacbio' and longReadTech != 'nano'):
            raise ValueError("Only recognized long read techs: 'nano' or 'pacbio'")

        self.shortMapper = shortMapper
        if (shortMapper != 'bwa' and shortMapper != 'bowtie'):
            raise ValueError("Only recognized short read mapper: 'bwa' or 'bowtie'")

        self.longMapper = 'minimap2'

        self.illuminaR1 = 'illumina_r1.fq.gz'
        self.illuminaR2 = 'illumina_r2.fq.gz'
        self.nanopore = 'nanopore_reads.fastq.gz'
        self.pacbio = 'pacbio_reads.fastq.gz'

        self.species = {} # Map of isolate names to their corresponding species
        harmonizeName = {'Kleb': "Kleb.", 'klepne': "Kleb.", 'klpene': "Kleb.", 'pseaer': "Pseud.","entclo":"Entclo.",
                         'Pseu': "Pseud.", 'Acin': "Acin.",'acibau': "Acin.", "Esch" : "Ecoli", "esccol": "Ecoli"}

        self.mlst = {}

        # Load in first file
        tbl = pd.read_csv(folder + "carbapenamase_seq_runs/prioritized_samples_carbapenemansen.csv",sep = ",")
        for index, row in tbl.iterrows():
            if (str(row['Internal #']) != "nan" and str(row['Internal #']) != " "):
                isolate = "carb%03d"%(int(row['Internal #']))
                self.species[isolate] = harmonizeName[row['Species']]

        tbl = pd.read_csv(folder + "carbapenamase_seq_runs/second_sample_list.csv",sep = ",")
        for index, row in tbl.iterrows():
            if (str(row['Internal #']) != "nan" and str(row['Internal #']) != " "):
                isolate = "carb%03d"%(int(row['Internal #']) + 60)
                self.species[isolate] = harmonizeName[row['Species']]

        tbl = pd.read_csv(folder + "carbapenamase_seq_runs/ST_types.tsv",sep = "\t")
        for index, row in tbl.iterrows():
                isolate = row['isolate']
                self.mlst[isolate] = row['ST']

    def assemblyName(self,strainName,cmpr=None):
        if (cmpr is None):
            assembler = self.assembler
        else:
            assembler = cmpr

        if (assembler == 'unicycler'):
            folderName = strainName + "/unicycler_hybrid_" + self.longReadTech + "/"
            assembleName = 'assembly.fasta'
        elif (assembler == 'canu'):
            if (os.path.exists(strainName+"/canu_" + self.longReadTech + "/pilon/")):
                folderName =  strainName + "/canu_" + self.longReadTech + "/pilon/"
                assembleName = "assembly.fasta"
            else:
                folderName = strainName + "/canu_" + self.longReadTech + "/"
                assembleName = "assembly.contigs.fasta"
        elif (assembler == 'hybrid'):
            if (os.path.exists(strainName + '/merged/pilon/')):
                folderName = strainName + '/merged/pilon/'
                assembleName = "assembly_filtered.fasta"
            else:
                folderName = strainName + '/merged/'
                assembleName = "assembly.fasta"
        elif (assembler == 'final'):
            if (os.path.exists(strainName + '/final/pilon/')):
                folderName = strainName + '/final/pilon/'
                assembleName = "assembly.fasta"
            else:
                folderName = strainName + '/final/'
                assembleName = "assembly.fasta"
        else:
            raise ValueError("Not a recognized assembly name")

        fileName = folderName + assembleName
        return fileName,assembleName,folderName

    def plot_assemblyEntropy(self):
        canuN = []
        uniN = []
        fileSize = []
        labels = []

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            isolateName = strainName.split('/')[-1]

            canuName,_,_ = self.assemblyName(strainName,'canu')
            uniName,_,_ = self.assemblyName(strainName,'unicycler')

            if (os.path.exists(canuName) and os.path.exists(uniName)):
                fileSize.append(os.path.getsize(strainName+"/"+self.nanopore)/1e7)

                with open(canuName, 'r') as input_handle:
                    contigLen = []
                    for record in SeqIO.parse(input_handle,"fasta"):
                        contigLen.append(1.*len(record.seq))
                    contigLen = np.array(contigLen)
                    contigLen = contigLen / (1. * np.sum(contigLen))
                    canuN.append(-np.sum( contigLen * np.log(contigLen)))

                with open(uniName, 'r') as input_handle:
                    contigLen = []
                    for record in SeqIO.parse(input_handle,"fasta"):
                        contigLen.append(1.*len(record.seq))
                    contigLen = np.array(contigLen)
                    contigLen = contigLen / (1. * np.sum(contigLen))
                    uniN.append(-np.sum( contigLen * np.log(contigLen )))

                labels.append(isolateName)
        return np.array(canuN),np.array(uniN)

    def plot_compare_assemblyPileUp(self):
        import os
        import pickle as pickle

        if (not os.path.exists('assembly_plot_data.pkl')):
            import json

            c_mean_cov = []
            c_rel_std_cov = []

            c_frac_uncov = []
            c_frac_unmapped = []

            c_avgSiteEntropy = []

            u_mean_cov = []
            u_rel_std_cov = []

            u_frac_uncov = []
            u_frac_unmapped = []

            u_avgSiteEntropy = []

            for strain in self.isolateFolders:
                # print strain
                strainName = strain.rstrip('/')
                _,_,canuFolder = self.assemblyName(strainName,cmpr='canu')
                _,_,uniFolder = self.assemblyName(strainName,cmpr='unicycler')

                outCDir = canuFolder + self.shortMapper + "/pileup/"
                c_pileUp = outCDir + "pileUp.json"
                c_unmapped = outCDir + "unmapped.txt"

                outUDir = uniFolder + self.shortMapper + "/pileup/"
                u_pileUp = outUDir + "pileUp.json"
                u_unmapped = outUDir + "unmapped.txt"

                if (os.path.exists(c_pileUp) and os.path.exists(c_unmapped) and os.path.exists(u_pileUp) and os.path.exists(u_unmapped) ):
                    with open(c_pileUp,'r') as pcUp, open(c_unmapped,'r') as cUm, \
                        open(u_pileUp,'r') as puUp, open(u_unmapped,'r') as uUm:
                        # print um.read()
                        c_frac_unmapped.append( float(cUm.read()) )
                        u_frac_unmapped.append( float(uUm.read()) )

                        ### CANU ###
                        cov = json.load(pcUp)

                        contig_len = np.zeros(len(cov))
                        mean_cov = np.zeros(len(cov))
                        site_entropy = np.zeros(len(cov))

                        for n in xrange(len(cov)):
                            contig_len[n] = len(cov[n]['cov'])
                            mean_cov[n] = cov[n]['avg_cov']
                            site_entropy[n] = cov[n]['site_entropy']

                        c_mean_cov.append( np.average(mean_cov,weights = contig_len) )
                        c_avgSiteEntropy.append( np.average(site_entropy,weights = contig_len) )

                        ### UNICYCLER ###
                        cov = json.load(puUp)

                        contig_len = np.zeros(len(cov))
                        mean_cov = np.zeros(len(cov))
                        site_entropy = np.zeros(len(cov))

                        for n in xrange(len(cov)):
                            contig_len[n] = len(cov[n]['cov'])
                            mean_cov[n] = cov[n]['avg_cov']
                            site_entropy[n] = cov[n]['site_entropy']

                        u_mean_cov.append( np.average(mean_cov,weights = contig_len) )
                        u_avgSiteEntropy.append( np.average(site_entropy,weights = contig_len) )

            c_mean_cov = np.array(c_mean_cov)
            c_rel_std_cov = np.array(c_rel_std_cov)
            c_frac_uncov = np.array(c_frac_uncov)
            c_frac_unmapped = np.array(c_frac_unmapped)
            c_avgSiteEntropy = np.array(c_avgSiteEntropy)

            u_mean_cov = np.array(u_mean_cov)
            u_rel_std_cov = np.array(u_rel_std_cov)
            u_frac_uncov = np.array(u_frac_uncov)
            u_frac_unmapped = np.array(u_frac_unmapped)
            u_avgSiteEntropy = np.array(u_avgSiteEntropy)

            return (c_mean_cov,c_frac_unmapped,c_avgSiteEntropy),(u_mean_cov,u_frac_unmapped,u_avgSiteEntropy)

        else:
            with open('assembly_plot_data.pkl','r') as pkl:
                x = pickle.load(pkl)
            return x

    def plot_final_assemblyPileUp(self):
        import os
        import json

        ill_sEnt = []
        nan_sEnt = []

        ill_mean_cov = []
        nan_mean_cov = []

        ill_std_cov = []
        nan_std_cov = []

        ill_unmapped = []
        nan_unmapped = []

        ill_uncov = []
        nan_uncov = []

        labels = []
        nan_size = []
        for strain in self.isolateFolders:
            # print strain
            strainName = strain.rstrip('/')
            _,_,folder = self.assemblyName(strainName,cmpr='final')

            outIllDir = folder + self.shortMapper + "/pileup/"
            ill_pileUp = outIllDir + "pileUp.json"
            ill_summ = outIllDir + "summary.json"
            ill_txt = outIllDir + "unmapped.txt"

            outNanoDir = folder + self.longMapper + "/pileup/"
            nan_pileUp = outNanoDir + "pileUp.json"
            nan_summ = outNanoDir + "summary.json"
            nan_txt = outNanoDir + "unmapped.txt"

            if (os.path.exists(nan_pileUp) and os.path.exists(ill_pileUp)):
                labels.append(strainName.split("/")[-1])
                with open(ill_pileUp,'r') as iUp, open(ill_summ,'r') as iSum, open(ill_txt,'r') as iTxt, \
                     open(nan_pileUp,'r') as nUp, open(nan_summ,'r') as nSum, open(nan_txt,'r') as nTxt:

                    ill_unmapped.append( float(iTxt.read()) )
                    nan_unmapped.append( float(nTxt.read()) )

                    ### ILLUMINA ###
                    ill = json.load(iSum)

                    contig_len = np.zeros(len(ill))
                    mean_cov = np.zeros(len(ill))
                    std_cov = np.zeros(len(ill))
                    site_entropy = np.zeros(len(ill))
                    uncov = np.zeros(len(ill))

                    for n in xrange(len(ill)):
                        contig_len[n] = ill[n]['contig_len']
                        mean_cov[n] = ill[n]['avg_cov']
                        std_cov[n] = ill[n]['std_cov']
                        uncov[n] = ill[n]['num_uncov']
                        site_entropy[n] = ill[n]['site_entropy']

                    ill_mean_cov.append( np.average(mean_cov,weights = contig_len) )
                    ill_std_cov.append( np.average(std_cov,weights = contig_len) )
                    ill_sEnt.append( np.average(site_entropy,weights = contig_len) )
                    ill_uncov.append( np.sum(uncov) / (1.0 * np.sum(contig_len)) )

                    ### NANOPORE ###
                    nan_size.append(os.path.getsize(folder + self.longMapper + '/mappedNanopore.bam'))
                    nan = json.load(nSum)

                    contig_len = np.zeros(len(nan))
                    mean_cov = np.zeros(len(nan))
                    std_cov = np.zeros(len(nan))
                    site_entropy = np.zeros(len(nan))
                    uncov = np.zeros(len(nan))

                    for n in xrange(len(nan)):
                        contig_len[n] = nan[n]['contig_len']
                        mean_cov[n] = nan[n]['avg_cov']
                        std_cov[n] = nan[n]['std_cov']
                        uncov[n] = nan[n]['num_uncov']
                        site_entropy[n] = nan[n]['site_entropy']

                    # if
                    nan_mean_cov.append( np.average(mean_cov, weights = contig_len) )
                    nan_std_cov.append( np.average(std_cov, weights = contig_len) )
                    nan_sEnt.append( np.average(site_entropy, weights = contig_len) )
                    nan_uncov.append( np.sum(uncov) / (1.0 * np.sum(contig_len)) )

        ill_mean_cov = np.array(ill_mean_cov)
        ill_std_cov = np.array(ill_std_cov)
        ill_sEnt = np.array(ill_sEnt)
        ill_unmapped = np.array(ill_unmapped)
        ill_uncov = np.array(ill_uncov)

        nan_mean_cov = np.array(nan_mean_cov)
        nan_std_cov = np.array(nan_std_cov)
        nan_sEnt = np.array(nan_sEnt)
        nan_unmapped = np.array(nan_unmapped)
        nan_uncov = np.array(nan_uncov)

        return (ill_mean_cov,ill_std_cov,ill_sEnt,ill_unmapped,ill_uncov),(nan_mean_cov,nan_std_cov,nan_sEnt,nan_unmapped,nan_uncov),np.array(labels),np.array(nan_size)

    # @jit
    # def computeMinorVariantFraction(self,cov):


    #     uncoveredLength = 0
    #     totalLength = 0

    #     for contig in cov:
    #         minorVariantFraction = np.array([ 1.0*np.max(x)/np.sum(x) for x in contig['cov'] if np.sum(x) > 0 ])
    #     #     uncoveredLength += sum(1 for x in contig['cov'] if np.sum(x)==0)
    #     #     totalLength += len(contig['cov'])
    #     # print(totalLength)
    #     return "hi"

    def storeFinalAssemblyCovStats(self, outFile="cov.pkl"):
        import os
        import json
        import pickle
        # import itertools
        # from collections import Counter

        assembly_cov = {}
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            _,_,folder = self.assemblyName(strainName,cmpr='final')
            covFile = folder + "covStats.json"

            if ('pilon' in folder and os.path.exists(covFile)):
                isoName = strainName.split("/")[-1]
                with open(covFile,'r') as inFile:
                    assembly_cov[isoName] = json.load(inFile)

        with open(outFile,'wb+') as outH:
            pickle.dump(assembly_cov,outH,protocol=pickle.HIGHEST_PROTOCOL)

        # if (not os.path.exists(outFile)):
        #     assembly_cov = {}
        #     for strain in self.isolateFolders:
        #         strainName = strain.rstrip('/')
        #         _,_,folder = self.assemblyName(strainName,cmpr='final')

        #         outIllDir = folder + self.shortMapper + "/pileup/"
        #         ill_pileUp = outIllDir + "pileUp.json"
        #         ill_txt = outIllDir + "unmapped.txt"

        #         outNanoDir = folder + self.longMapper + "/pileup/"
        #         nan_pileUp = outNanoDir + "pileUp.json"
        #         nan_txt = outNanoDir + "unmapped.txt"

        #         if (os.path.exists(nan_pileUp) and os.path.exists(ill_pileUp) and "pilon" in folder):  #If we have mapped data for this isolate, continue
        #             isoName = strainName.split("/")[-1]
        #             assembly_cov[isoName] = {}
        #             assembly_cov[isoName]['species'] = self.species[isoName]
        #             assembly_cov[isoName]['ill'] = {}
        #             assembly_cov[isoName]['nan'] = {}

        #             with open(ill_pileUp,'r') as iUp, open(ill_txt,'r') as iTxt, \
        #                  open(nan_pileUp,'r') as nUp, open(nan_txt,'r') as nTxt:

        #                 assembly_cov[isoName]['ill']['frac_unmapped'] = float(iTxt.read())
        #                 # assembly_cov[isoName]['nan']['frac_unmapped'] = float(nTxt.read())

        #                 ### ILLUMINA PILEUP ###
        #                 print "Illumina"
        #                 cov = json.load(iUp)
        #                 print "Illumina Loaded"
        #                 minorVariantFraction = [ [ 1.0*np.max(x)/np.sum(x) for x in contig['cov'] if np.sum(x) > 0 ] for contig in cov ]
        #                 uncoveredLength = sum([ sum(1 for x in contig['cov'] if np.sum(x)==0) for contig in cov ])
        #                 totalLength = sum([ len(contig['cov']) for contig in cov ])
        #                 assembly_cov[isoName]['ill']['minor'] = Counter(itertools.chain.from_iterable(minorVariantFraction))
        #                 assembly_cov[isoName]['ill']['uncovered'] = uncovered

        #                 ### NANOPORE PILEUP ###
        #                 print "Nanopore"
        #                 cov = json.load(nUp)
        #                 print "Nanopore Loaded"
        #                 nanoCov = [ [ np.sum(x) for x in contig['cov'] ] for contig in cov ]
        #                 assembly_cov[isoName]['nan']['cov'] = Counter(itertools.chain.from_iterable(nanoCov))

    def number_of_closed(self):
        from collections import defaultdict

        ntot = 0
        nclosed = 0
        closed_strains = defaultdict(lambda: False)
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')

            name,_,_ = self.assemblyName(strainName)
            isolate = strainName.split('/')[-1]
            if (os.path.exists(name)):
                with open(name,'r') as fasta:
                    genomeHeader = fasta.readline()
                    if ("circular=true" in genomeHeader.lower()):
                        nclosed += 1
                        closed_strains[isolate] = True
                    else:
                        closed_strains[isolate] = False
                ntot += 1

        print( nclosed / (1. * ntot) )
        return closed_strains

    def geneLength_compare_to_uniprot(self,canvas):
        from statsmodels.distributions.empirical_distribution import ECDF
        # import matplotlib.pylab as plt

        binEdges = np.array( [0] + np.linspace(.5,.83,5).tolist() + np.linspace(.86,1.05,20).tolist() + np.linspace(1.1,1.2,4).tolist() + [2.0] )
        # x90 = .5*(binEdges[0:-1] + binEdges[1:])
        # x90 = np.argmin( np.abs(.9-x90) )

        speciesColors = toyplot.color.brewer.map("Set1",count=4)
        spIndex = {"Ecoli":0,"Acin.":1,"Pseud.":2,"Kleb.":3}

        binWidth = (1.05-.86)/20
        axes = canvas.cartesian(margin=60,padding=2,yscale='log',xmax=1.1,xmin=.5,ymin=.02,ymax=1,xlabel='Ratio of gene length to Uniprot match',ylabel='Fraction of genes in bin') #,label="Quality of Annotation")

        ecoliL = None
        acinL = None
        pseuL = None
        klebL = None

        ecoli90 = []
        acin90 = []
        pseu90 = []
        kleb90 = []
        for nS,strain in enumerate(self.isolateFolders):
            strainName = strain.rstrip('/')
            _,_,folderName = self.assemblyName(strainName)
            uniProtName = glob.glob(folderName+'prokka/swissProt.m8')
            queryFile = glob.glob(folderName+'prokka/*.faa')

            if (len(uniProtName) == 1 and len(queryFile) == 1):
                blastTable = pd.read_csv(uniProtName[0],sep='\t',header=None,
                             names=['query','db','pmatch','length','nMismatch','nGapOpen','qStart','qEnd','sStart','sEnd','eVal','bitScore'])
                queryTable = ph.read_fasta(queryFile[0])
                blastTable = blastTable.set_index("query")
                # queryTable.set_index("name",drop=True,inplace=True)
                # queryTable = queryTable.to_dict(orient="index")

                fracCov = []
                numDiff = 0

                for index,row in queryTable.iterrows():
                    if (row['name'] in blastTable.index):
                        mapLengths = blastTable.loc[row['name'],'length']
                        prcIds = blastTable.loc[row['name'],'pmatch']
                        # print row['name']
                        # print mapLengths
                        # print prcIds

                        if (mapLengths.size > 1):
                            mapLengths = mapLengths.values
                            prcIds = prcIds.values
                            ind = np.argmax( np.array( [ float(p) for p in prcIds] ) * np.array([ float(l) for l in mapLengths ]) )
                            mapLengths = mapLengths[ind]
                            prcIds = prcIds[ind]

                    # if (blastTable[blastTable['query'] == row['name']].size > 0):
                    #     topMatch = blastTable[blastTable['query'] == row['name']].iloc[0]
                        if (prcIds > 80):
                            fracCov.append(1.*mapLengths/len(row['sequence']))
                            if (abs(fracCov[-1] -1) > .1):
                                numDiff += 1

                if (len(fracCov) > 0):
                    # print numDiff
                    # print numDiff/(1.*index)
                    # print strain

                    fracCov = np.array(fracCov)
                    Y,X = np.histogram(fracCov,binEdges,density=True)
                    Y *= binWidth

                    # ind = np.argmin(np.abs(np.array(X)-.9))
                    # print index*Y[int(ind)]
                    isoName = strainName.split('/')[-1]
                    h = axes.plot(X[:-1],Y+.01,opacity=.5,color=speciesColors.colors(spIndex[self.species[isoName]]))

                    ecdf = ECDF(fracCov)
                    X,Y = self.returnECDFLine(ecdf)
                    X90 = np.array(X)
                    X90 = np.argmin( np.abs(X90-.9) )
                    fG = Y[X90]
                    if (spIndex[self.species[isoName]] == 0):
                        ecoli90.append(fG)
                    elif (spIndex[self.species[isoName]] == 1):
                        acin90.append(fG)
                    elif (spIndex[self.species[isoName]] == 2):
                        pseu90.append(fG)
                    elif(spIndex[self.species[isoName]] == 3):
                        kleb90.append(fG)

                    # print isoName
                    # h = axes.plot(X[::20],Y[::20],opacity=.5,color=speciesColors.colors(spIndex[self.species[isoName]]))
                    if (spIndex[self.species[isoName]] == 0 and ecoliL is None):
                        ecoliL = h
                    elif (spIndex[self.species[isoName]] == 1 and acinL is None):
                        acinL = h
                    elif (spIndex[self.species[isoName]] == 2 and pseuL is None):
                        pseuL = h
                    elif(spIndex[self.species[isoName]] == 3 and klebL is None):
                        klebL = h

                # if (nS > 60):
                #     break

        species_markers = [("Ecoli",ecoliL),("Acin.",acinL),("Pseud.",pseuL),("Kleb.",klebL)]
        _,bin_edges = np.histogram(ecoli90+acin90+pseu90+kleb90,bins=12)
        hE,_ = np.histogram(ecoli90,bins=bin_edges)
        hA,_ = np.histogram(acin90,bins=bin_edges)
        hP,_ = np.histogram(pseu90,bins=bin_edges)
        hK,_ = np.histogram(kleb90,bins=bin_edges)

        barH = np.column_stack((hE,hA,hP,hK))
        embed = canvas.cartesian(corner=("top-left", 120, 175, 175),xlabel="Prc genes lt .9",ylabel="# of strains")
        embed.x.ticks.locator = toyplot.locator.Explicit(locations=[1,2,3,4,5],labels=["1","2","3","4","5"]) #,format="{:.2f}")
        embed.bars(100*bin_edges[:-1],barH,color=[speciesColors.colors(nn) for nn in range(4)])
        # legend = canvas.legend(species_markers,corner=("top-right",70,10,50))
        legend = None
        return axes,embed,legend

    def number_of_completed(self,X=.9):
        # import bisect

        contiguity = []
        num_contigs = []
        gen_len = []
        species = []

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            name,_,folder = self.assemblyName(strainName)

            if (os.path.exists(name) and "pilon" in folder):
                isoName = strainName.split('/')[-1]
                ncontig = 0
                lcontig = 0
                contig_lens = []
                with open(name,'r') as fasta:
                    genomeHeader = fasta.readline()
                    while genomeHeader:
                        if (">" == genomeHeader[0]):
                            ncontig += 1
                            if (lcontig > 0):
                                contig_lens.append(lcontig)
                            lcontig = 0
                        else:
                            lcontig += len(genomeHeader.strip())
                        genomeHeader = fasta.readline()

                contig_lens.append(lcontig)
                contig_lens = np.array(contig_lens)
                Lgenome = np.sum(contig_lens)
                lcumsum = np.cumsum(contig_lens) / (1. * Lgenome)

                contiguity.append(lcumsum[0])
                num_contigs.append(ncontig)
                gen_len.append(Lgenome)
                species.append(self.species[isoName])

        return np.array(num_contigs), np.array(contiguity), np.array(gen_len), np.array(species)

    def plotStructureTrees(self,forceComp=False,save=False):
        import json
        from Bio import SeqIO
        import toyplot.pdf
        import toyplot.png

        key_prefix = "aln_geneFamily/"
        gbkLoc = "/final/pilon/prokka/entero.gbk"
        spIndex = {"Ecoli":0,"Acin.":1,"Pseud.":2,"Kleb.":3}

        saveDir = "/home/nicholas/clusterData/papers/carbapenem_manuscript/figs/synteny/core_dump/"

        for molecule in self.moleculeFolders:
            # if ("KPC" in molecule):
            if ("TEM" in molecule or "CTX-M-15" in molecule or "SHV" in molecule):
                print(molecule)
                synTree = molecule+"syntenyCluster.newick"
                phyTree = molecule+"/nc/RAxML_bestTree.nc_rax.newick"
                if (not os.path.exists(phyTree)):
                    phyTree = molecule+"/nc/nc_ft.newick"

                divClusters = self.cutLongBranches(phyTree)
                tree = toytree.tree(synTree, format=0)

                tips = list(tree.get_tip_labels())
                label = tree.get_node_values(feature='name', show_root=True, show_tips=True)
                rows = np.flip(np.array( [x.split('||')[0] for x in label if x in tips] ),0)

                species = np.array( [spIndex[self.species[row.split('_')[0]]] for row in rows] )
                cluster = np.array( [divClusters[row] for row in rows] )

                mlsts = [ self.mlst[row.split('_')[0]] for row in rows ]
                mlsts_id = { st:n for n,st in enumerate(set(mlsts))}
                mlsts = np.array( [mlsts_id[st] for st in mlsts] )

                family = molecule.split('/')[-2]
                genus = family.split('_')[0].split('-')[0]
                family = family.replace("_", " ")

                contigJSON = molecule + 'contigs.json'
                if (not os.path.exists(contigJSON) or forceComp):
                    with open(self.clusterJSON,'r') as GC:
                        gClusters = json.load(GC)
                        members = gClusters[key_prefix + genus][family]['members']
                        contigs = {}
                        numIso = defaultdict(lambda: 0)
                        for member in members:
                            iso = member[0]
                            locus = member[1].replace(" ","_")
                            with open(self.topFolder + iso + gbkLoc,'r') as gbk:
                                for n,contig in enumerate(SeqIO.parse(gbk,'genbank')):
                                    for segment in contig.features:
                                        if (segment.type == 'CDS' and segment.qualifiers['locus_tag'][0] == locus):
                                            index = numIso[iso]
                                            # index = np.where(loci==segment.qualifiers['locus_tag'][0])[0][0]
                                            isoName = iso + "_%02d"%index
                                            contigs[isoName] = (n,len(contig.seq))
                                            numIso[iso] += 1

                    with open(contigJSON,'w+') as contigFile:
                        json.dump(contigs,contigFile)
                else:
                    with open(contigJSON,'r') as rJson:
                        contigs = json.load(rJson)

                # print contigs.keys()
                # with open(self.clusterJSON,'r') as GC:
                #     gClusters = json.load(GC)
                #     members = gClusters[key_prefix + genus][family]['members']
                #     print members
                contig = np.array( [1 * ( (contigs[row][1] > 1000) * (contigs[row][1] < 5e5) ) for row in rows] )
                dataMatrix = np.vstack((cluster,contig,mlsts,species)).transpose()
                labels = ['Clu.','Pla.','ST','Sp.','']

                numCluster = np.max(cluster) + 1
                # if (numCluster < 8):
                #     clusterMap = toyplot.color.brewer.map("Set2",count=8)
                # else:
                #     if (numCluster < 24):
                #         C = 3
                #     elif (numCluster < 32):
                #         C = 4
                #     else:
                #         C = 5
                # basicMap = toyplot.color.brewer.map("Set2",count=8)
                # clusterMap = toyplot.color.spread(basicMap.color(0),lightness=0.5,count=C) + \
                #             toyplot.color.spread(basicMap.color(1),lightness=0.5,count=C) + \
                #             toyplot.color.spread(basicMap.color(2),lightness=0.5,count=C) + \
                #             toyplot.color.spread(basicMap.color(3),lightness=0.5,count=C) + \
                #             toyplot.color.spread(basicMap.color(4),lightness=0.5,count=C) + \
                #             toyplot.color.spread(basicMap.color(5),lightness=0.5,count=C) + \
                #             toyplot.color.spread(basicMap.color(6),lightness=0.5,count=C) + \
                #             toyplot.color.spread(basicMap.color(7),lightness=0.5,count=C)
                binNum = int( np.ceil(numCluster / 7.) )

                if (binNum > 1):
                    clusterMap = toyplot.color.spread("Indigo",lightness=0.5,count=binNum) + \
                                toyplot.color.spread("LightSteelBlue",lightness=0.5,count=binNum) + \
                                toyplot.color.spread("MediumAquamarine",lightness=0.5,count=binNum) + \
                                toyplot.color.spread("LimeGreen",lightness=0.5,count=binNum) + \
                                toyplot.color.spread("Gold",lightness=0.5,count=binNum) + \
                                toyplot.color.spread("OrangeRed",lightness=0.5,count=binNum) + \
                                toyplot.color.spread("Red",lightness=0.5,count=binNum)
                else:
                    clusterMap = toyplot.color.Palette( [toyplot.color.css("DarkOrchid"), \
                                                            toyplot.color.css("SteelBlue"), \
                                                            toyplot.color.css("Teal"), \
                                                            toyplot.color.css("LightGreen"), \
                                                            toyplot.color.css("Gold"), \
                                                            toyplot.color.css("OrangeRed"), \
                                                            toyplot.color.css("Red")]  )

                plasmidMap = toyplot.color.brewer.map("Set2",count=3)
                speciesMap = toyplot.color.brewer.map("Set1",count=4)
                stMap_tmp = toyplot.color.brewer.palette("Spectral",count=8)
                binNum = len(set(mlsts))/8.
                if (binNum <= 1):
                    stMap = stMap_tmp
                else:
                    C = int(np.ceil(binNum))
                    stMap = toyplot.color.spread(stMap_tmp.color(0),lightness=0.5,count=C) + \
                            toyplot.color.spread(stMap_tmp.color(1),lightness=0.5,count=C) + \
                            toyplot.color.spread(stMap_tmp.color(2),lightness=0.5,count=C) + \
                            toyplot.color.spread(stMap_tmp.color(3),lightness=0.5,count=C) + \
                            toyplot.color.spread(stMap_tmp.color(4),lightness=0.5,count=C) + \
                            toyplot.color.spread(stMap_tmp.color(5),lightness=0.5,count=C) + \
                            toyplot.color.spread(stMap_tmp.color(6),lightness=0.5,count=C) + \
                            toyplot.color.spread(stMap_tmp.color(7),lightness=0.5,count=C)


                # print mlsts
                # print mlsts_id
                # print len(set(mlsts))
                # print stMap
                colors = [clusterMap,plasmidMap,stMap,speciesMap]
                if ('NDM' in molecule or "OXA-48" in molecule):
                    canvas = self.plotTree_w_attributes(tree,dataMatrix,labels,colors,H=500,colGap=4)
                else:
                    canvas = self.plotTree_w_attributes(tree,dataMatrix,labels,colors,H=20*dataMatrix.shape[0],colGap=4)

                # with open(molecule+"synClusters.pdf","wb+") as img:
                name = molecule.split("/")[-2]
                # print molecule
                # print name
                if (save):
                    toyplot.pdf.render(canvas,saveDir+name+".pdf")
                    toyplot.png.render(canvas,saveDir+name+".png")
                    # break

    def cutLongBranches(self,treeFile,branch_cutoff=.001):
        tree = toytree.tree(treeFile, format=0)

        assignedLeaves = set([])
        clusters = []
        for node in tree.tree.traverse( strategy= "postorder" ):
            if (node.dist > branch_cutoff):
                leaves = set([leaf.name for leaf in node.get_leaves()]) - assignedLeaves
                assignedLeaves.update(leaves)
                clusters.append(leaves)
        allLeaves = set( [leaf.name for leaf in node.get_leaves()] )
        remainingLeaves = allLeaves
        remainingLeaves.difference_update(assignedLeaves)
        clusters.append(remainingLeaves)

        cAssign = {}
        for n,clusterMembers in enumerate(clusters):
            for member in clusterMembers:
                cAssign[member.split('||')[0]] = n

        return cAssign

    def plotTree_w_attributes(self,tree,matrix,tlabel,colors,W=400,H=400,M=50,frac=.68,rowH=.9,colGap=20):

        canvas = toyplot.Canvas(width=W,height=H)
        tree_axes = canvas.cartesian(bounds=(M,frac*(W-M),M,H-M),xlabel="Normalized edit distance")
        for node in tree.tree.traverse():
            if node.is_leaf():
                node.color= 'Grey'
            else:
                node.color= 'WhiteSmoke'
        nodeColors = tree.get_node_values('color', show_root=1, show_tips=1)

        _,tree_axes = tree.draw(axes=tree_axes,
                                node_labels=None,
                                node_size=None,
                                node_color=nodeColors,
                                node_style={"stroke":"black"},
                                use_edge_lengths=True,
                                tip_labels_align=True,
                                tip_labels=False,
                                padding=None,
                                show_tips=True,
                                show_root=False,
                                )

        tree_axes.y.show = False
        verts = tree.verts

        # Subsample verts for just leaves.
        tips = set(tree.get_tip_labels())
        label = tree.get_node_values(feature='name', show_root=True, show_tips=True)
        indx = np.array( [x in tips for x in label] )
        verts = verts[indx,:]

        tree_ymin = tree_axes.__dict__['_ymin_range'] #+ tree_axes.__dict__['_padding']/2
        tree_ymax = tree_axes.__dict__['_ymax_range'] #- tree_axes.__dict__['_padding']/2

        # Linearly interpolate tree vertex position into tree domain
        vert_max = np.max(verts[:,1])
        vert_min = np.min(verts[:,1])

        vertY = ((tree_ymax - tree_ymin) / (vert_max - vert_min) ) * ( verts[:,1] - vert_min ) + tree_ymin
        rowH = rowH * np.mean( np.abs(np.diff(vertY)))

        topRow = vertY - rowH/2
        botRow = vertY + rowH/2

        # print(topRow)
        # print(botRow)

        tableMin = np.min(topRow) - rowH
        tableMax = np.max(botRow)

        # N = verts.shape[0]
        # matrix = np.random.rand( N,4 )

        # Enlarge matrix.
        matrix = toyplot.require.scalar_matrix(matrix)
        # en_matrix = np.zeros( (2*matrix.shape[0]-1,matrix.shape[1]+1) )
        # for i in np.arange(matrix.shape[0]):
        #     en_matrix[2*i,:-1] = matrix[i,:]

        # colormap = toyplot.color.brewer.map("BlueRed", domain_min=matrix.min(), domain_max=matrix.max())
        # colors = colormap.colors(matrix)

        # Build up the attribute table incrementally
        xmin_range, xmax_range, _, _ = toyplot.layout.region(
            0, canvas._width, 0, canvas._height, bounds=(frac*(W-M),W-M,M,H-M), rect=None, corner=None, grid=None, margin=None)
        xmin_range += 10

        table = toyplot.coordinates.Table(
                             rows=2*matrix.shape[0]-1,
                             columns=matrix.shape[1]+1,
                             annotation=False,
                             brows=0,
                             trows=1,
                             rcolumns=0,
                             lcolumns=0,
                             xmax_range=xmax_range,
                             xmin_range=xmin_range,
                             ymax_range=tableMax,
                             ymin_range=tableMin,
                             parent=canvas,
                             label=None,
                             filename=None
                            )
        # Insert labels on top
        # tlabel = ["Col 1","Col 2","Col 3","Col 4"]

        for n,tRow in enumerate(table.top.column):
            # if (n%2 ==0):
            tRow.data = tlabel[n]
            tRow.height = rowH

        # table.body.cells.data = en_matrix
        table.body.cells.format = toyplot.format.NullFormatter()

        # Color each cell.
        for i in np.arange(matrix.shape[0]):
            for j in np.arange(matrix.shape[1]):
                cell = table.body.cell[2*i, j]
                # cell.style = {"stroke": "black", "stroke-width":.5, "fill": toyplot.color.to_css(colors[i, j])}
                if (j != 1):
                    cell.style = {"stroke": "black", "stroke-width":.5, "fill": toyplot.color.to_css(colors[j].color(matrix[i,j]))}
                else:
                    if (matrix[i,j] == 0):
                        cell.style = {"stroke": "black", "stroke-width":.5, "fill": "Wheat"}
                    else:
                        cell.style = {"stroke": "black", "stroke-width":.5, "fill": "Maroon"}

                cell.title = matrix[i, j]

        J = matrix.shape[1]
        tip_names = tree.get_tip_labels()
        for i in np.arange(matrix.shape[0]):
            cell = table.body.cell[2*i,J]
            cell.data = tip_names[i].replace("carb","ci").split("_")[0]
            # cell.style = {"fill":"white"}
            cell.lstyle = {"fill":"black","font-size":8}
            cell.format = toyplot.format.FloatFormatter("{:.1f}")

        # Change each cell's geometry.
        for i in np.arange(2*matrix.shape[0]-1):
            for j in np.arange(matrix.shape[1]+1):
                cell = table.body.cell[i, j]
                if (i % 2 == 0):
                    cell.height = rowH
                else:
                    cell.height = topRow[i/2+1] - botRow[i/2]

        table.body.gaps.columns[...] = colGap
        canvas._children.append(table)

        return canvas

    def gather_and_plot_basic_strain_data(self,W=600,H=300):
        import csv
        from collections import defaultdict, Counter
        import json

        import numpy as np
        import toyplot

        STinfo = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/ST_types.tsv"
        plasmidInfo = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/plasmid_types_detailed_v2.json"
        seqFile1 = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/prioritized_samples_carbapenemansen.csv"
        seqFile2 = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/second_sample_list.csv"
        geneFile = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/blaGenes.txt"

        strain_data = {"kpneumoniae":{"count":Counter(),"plasmid":Counter(),"plasmidNum":Counter(),"date":Counter(),"gene":Counter(),"closed":0,"total":0}, \
                       "ecoli":{"count":Counter(),"plasmid":Counter(),"plasmidNum":Counter(),"date":Counter(),"gene":Counter(),"closed":0,"total":0}, \
                       "paeruginosa":{"count":Counter(),"plasmid":Counter(),"plasmidNum":Counter(),"date":Counter(),"gene":Counter(),"closed":0,"total":0}, \
                       "abaumannii_2":{"count":Counter(),"plasmid":Counter(),"plasmidNum":Counter(),"date":Counter(),"gene":Counter(),"closed":0,"total":0} }
        species = {}
        displayName = {"kpneumoniae": "Kleb.", "ecoli":"Ecoli","paeruginosa":"Pseud.","abaumannii_2":"Acin."}

        spIndex = {"Ecoli":0,"Acin.":1,"Pseud.":2,"Kleb.":3}
        speciesColors = toyplot.color.brewer.map("Set1",count=4)

        with open(STinfo,'r') as metaData:
            rd = csv.reader(metaData,delimiter="\t")
            for n,row in enumerate(rd):
                if (n == 0):
                    species_index = row.index("species")
                    st_index = row.index("ST")
                    isolate_index = row.index("isolate")
                else:
                    if (row[st_index] != "-"):
                        strain_data[ row[species_index] ]['count'].update([int(row[st_index])])
                    else:
                        strain_data[ row[species_index] ]['count'].update([-1])
                    species[ row[isolate_index] ] = row[species_index]

        with open(plasmidInfo,'r') as plasmidData:
            plasmids = json.load(plasmidData)
            for isolate in plasmids.keys():
                strain_species = species[isolate]
                for iso_plasmids,dat in plasmids[isolate].items():
                    # for gene in dat["bla"]:
                        # strain_data[strain_species][gene] += 1
                    strain_data[strain_species]["plasmid"].update([x.split("(")[0] for x in dat["type"]])
                strain_data[strain_species]["plasmidNum"].update( [len(plasmids[isolate])] )

        with open(seqFile1,'r') as seqFile:
            rd = csv.reader(seqFile,delimiter=",")
            for n,row in enumerate(rd):
                if (n > 0):
                    if (row[0] != ' ' and row[0] != ""):
                        isolate = "carb%03d"%int(row[0])
                        if (isolate in species):
                            sp = species[isolate]
                            date = row[4]
                            if (date == "20xx"):
                                date = "2012"
                            strain_data[sp]["date"].update([date])

        with open(seqFile2,'r') as seqFile:
            rd = csv.reader(seqFile,delimiter=",")
            for n,row in enumerate(rd):
                if (n > 0):
                    if (row[0] != ' ' and row[0] != ""):
                        isolate = "carb%03d"%(60+int(row[0]))
                        if (isolate in species):
                            sp = species[isolate]
                            date = "20" + row[7].split(' ')[0].split('.')[-1]
                            strain_data[sp]["date"].update([date])

        with open(geneFile,'r') as gFile:
            for lineNum, line in enumerate(gFile):
                if ("db.fasta:" in line):
                    isolate = line[0:7]
                    gene = line.split("db.fasta:")[-1].split('_')[0].split('-')[0]
                    sp = species[isolate]
                    strain_data[sp]["gene"].update([gene])

        closedStrains = self.number_of_closed()
        for strain,closed in closedStrains.items():
            sp = species[strain]
            if (closed):
                strain_data[sp]["closed"] += 1
            strain_data[sp]["total"] += 1

        canvas = toyplot.Canvas(width=W,height=H)
        table = canvas.table(rows=4,columns=5,trows=1)

        table.top.cells.style = {"fill":"grey", "opacity":0.1}
        table.top.row[0].height = H/10
        table.top.grid.hlines[1,...] = "single"

        table.top.cell[0,0].data = "Species"
        table.top.cell[0,1].data = "# collected"
        table.top.cell[0,2].data = "Date collected"
        table.top.cell[0,3].data = "3 main MLSTs"
        table.top.cell[0,4].data = "Common bla genes"

        table.cells.column[0].width=W/10
        table.cells.column[1].width=W/10
        table.cells.column[2].width=W/4
        table.cells.column[3].width=W/5
        # table.cells.column[4].width=W/5

        table.body.grid.vlines[...,2:4] = "double"
        sample_dates = [ str(x) for x in range(2010,2018) ]
        for n,species in enumerate(strain_data.keys()):
            table.body.cell[n,0].data = displayName[species]
            table.body.cell[n,1].data = sum(strain_data[species]['count'].values())

            num_samples = np.array( [strain_data[species]['date'][x] for x in sample_dates ] )
            axes = table.body.cell[n,2].cartesian()
            axes.bars(num_samples,color=speciesColors.color(spIndex[displayName[species]]))

            types = [ x for x in strain_data[species]['count'].most_common(4) if x[0]>-1 ]
            table.body.cell[n,3].data = ','.join([str(x[0])+"("+str(x[1])+")" for x in types[0:3]])
            table.body.cell[n,4].data = ','.join([str(x[0]).lstrip('bla') for x in strain_data[species]['gene'].most_common(6)])

            # table.body.cell[n,3].data =
        # table.body.grid.hlines[0:4,...] = "single"

        # dat = np.random.normal(loc=1,size=25)
        # axes = table.body.cell[0,2].cartesian()
        # h = axes.bars(dat)

        # dat = np.random.normal(loc=1,size=25)
        # axes = table.body.cell[1,2].cartesian()
        # h = axes.bars(dat)

        # dat = np.random.normal(loc=1,size=25)
        # axes = table.body.cell[2,2].cartesian()
        # h = axes.bars(dat)

        # dat = np.random.normal(loc=1,size=25)
        # axes = table.body.cell[3,2].cartesian()
        # h = axes.bars(dat)

        return canvas,table,strain_data
        # df = pd.DataFrame(strain_data)
        # df.to_csv("/home/nolln/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/strain_overview.tsv",sep="\t")


    def buildPlasmidTable(self,strain_data,W=600,H=300):
        import numpy as np
        import toyplot

        displayName = {"kpneumoniae": "Kleb.", "ecoli":"Ecoli","paeruginosa":"Pseud.","abaumannii_2":"Acin."}

        spIndex = {"Ecoli":0,"Acin.":1,"Pseud.":2,"Kleb.":3}
        speciesColors = toyplot.color.brewer.map("Set1",count=4)

        canvas = toyplot.Canvas(width=W,height=H)
        table = canvas.table(rows=4,columns=4,trows=1)

        table.top.cells.style = {"fill":"grey", "opacity":0.1}
        table.top.row[0].height = H/10
        table.top.grid.hlines[1,...] = "single"

        table.top.cell[0,0].data = "Species"
        table.top.cell[0,1].data = "# closed"

        table.top.cell[0,2].data = "# plasmids/isolate"
        table.top.cell[0,3].data = "Dominant Inc groups"

        table.cells.column[0].width=W/10
        table.cells.column[1].width=W/10

        table.cells.column[2].width=W/5
        # table.cells.column[2].width=2*W/5

        table.body.grid.vlines[...,2:4] = "double"

        sample_num = range(7)
        for n,species in enumerate(strain_data.keys()):
            table.body.cell[n,0].data = displayName[species]
            table.body.cell[n,1].data =  str(strain_data[species]["closed"]) +  "/" + str(strain_data[species]["total"])
            num_samples = np.array( [strain_data[species]['plasmidNum'][x] for x in sample_num ] )
            axes = table.body.cell[n,2].cartesian()
            axes.bars(num_samples,color=speciesColors.color(spIndex[displayName[species]]))

            table.body.cell[n,3].data = ','.join([str(x[0]).lstrip('bla') for x in strain_data[species]['plasmid'].most_common(6)])

        return canvas,table

    def generate_inc_plasmid_figures(self,W=500,H=500):
        import csv
        import json
        from collections import Counter
        import itertools

        import numpy as np
        import seaborn as sns

        STinfo = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/ST_types.tsv"
        plasmidInfo = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/plasmid_types_detailed_v2.json"
        geneFile = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/blaGenes.txt"

        species = {}
        displayName = {"kpneumoniae": "Kleb.", "ecoli":"Ecoli","paeruginosa":"Pseud.","abaumannii_2":"Acin."}

        spIndex = {"Ecoli":0,"Acin.":1,"Pseud.":2,"Kleb.":3}
        speciesColors = toyplot.color.brewer.map("Set1",count=4)

        with open(STinfo,'r') as metaData:
            rd = csv.reader(metaData,delimiter="\t")
            for n,row in enumerate(rd):
                if (n == 0):
                    species_index = row.index("species")
                    isolate_index = row.index("isolate")
                else:
                    species[ row[isolate_index] ] = row[species_index]

        inc_matrix = Counter()
        inc_types = set()
        with open(plasmidInfo,'r') as plasmidData:
            plasmids = json.load(plasmidData)
            for isolate in plasmids.keys():
                if (len(plasmids[isolate].items()) > 1):
                    for typeA,typeB in itertools.combinations(plasmids[isolate].values(),r=2):
                        A = [x.split("(")[0] for x in typeA["type"]]
                        B = [x.split("(")[0] for x in typeB["type"]]
                        inc_types.update(A)
                        inc_types.update(B)
                        inc_matrix.update(list(itertools.product(A,B)))

                      # inc_matrix.update([x.split("(")[0] for x in dat["type"]])
        inc_types = {iType:n for n,iType in enumerate(inc_types)}
        plasmid_coccur = np.zeros( (len(inc_types),len(inc_types)) )

        for pair,amount in inc_matrix.items():
            plasmid_coccur[inc_types[pair[0]],inc_types[pair[1]]] = amount
            plasmid_coccur[inc_types[pair[1]],inc_types[pair[0]]] = amount

        return plasmid_coccur

    def gene_cooccurence(self,W=500,H=500,molecule=True):
        import json
        import numpy as np
        import toyplot
        import scipy.cluster.hierarchy as hcl
        import scipy.spatial.distance as ssd

        from collections import Counter
        import itertools

        # STinfo = "/home/nolln/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/ST_types.tsv"
        contigFile = self.topFolder + "carbapenamase_seq_runs/contig_meta.json"
        plasmidInfo = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/plasmid_types_detailed_v2.json"

        with open(contigFile,'r') as contigF:
            contigs = json.load(contigF)

        # species = {}
        # displayName = {"kpneumoniae": "Kleb.", "ecoli":"Ecoli","paeruginosa":"Pseud.","abaumannii_2":"Acin."}

        # spIndex = {"Ecoli":0,"Acin.":1,"Pseud.":2,"Kleb.":3}
        # speciesColors = toyplot.color.brewer.map("Set1",count=4)

        geneO = Counter()
        geneNames = Counter()
        with open(plasmidInfo,'r') as plasmidData:
            plasmids = json.load(plasmidData)
            for isolate,iso_plasmids in plasmids.items():
                if (molecule):
                    for p in iso_plasmids.keys():
                        genes = contigs[isolate][p.split(' ')[0]]["blaGenes"]
                        geneNames.update(genes)
                        for pair in itertools.combinations(genes,2):
                            geneO.update([tuple(sorted(pair))])
                else:
                    genes = set([])
                    for p in iso_plasmids.keys():
                        genes.update(contigs[isolate][p.split(' ')[0]]["blaGenes"])
                    geneNames.update(genes)
                    for pair in itertools.combinations(genes,2):
                        geneO.update([tuple(sorted(pair))])

        D = np.zeros( (len(geneNames), len(geneNames)) )
        gene_index = {name:n for n,name in enumerate(geneNames.keys())}
        for pair,value in geneO.items():
            D[ gene_index[pair[0]], gene_index[pair[1]] ] = value / np.sqrt( geneNames[pair[0]] * geneNames[pair[1]] )
            D[ gene_index[pair[1]], gene_index[pair[0]] ] = value / np.sqrt( geneNames[pair[0]] * geneNames[pair[1]] )

        # for key,value in geneNames.items():
        #     D[ gene_index[key], gene_index[key] ] = value

        Ddist = np.max(D)-D
        Ddist -= np.diag(np.diag(Ddist))
        Dsq = ssd.squareform(Ddist)
        Z = hcl.linkage(Dsq,method='average')
        dend = hcl.dendrogram(Z)

        index = dend['leaves']
        D = D[index,:]
        D = D[:,index]
        D += np.eye(D.shape[0])

        simplifiedNames = {'KPC-3':'KPC','NDM-4':'NDM', "CARB-5":"CARB","VIM-1":"VIM","PAM-1":"PAM","SHV-11":"SHV","DHA-7":"DHA"}
        simpleNames = np.array([simplifiedNames[x.split('_')[0].replace("bla","")] if x.split('_')[0].replace("bla","") 
            in simplifiedNames else x.split('_')[0].replace("bla","")  for x,val in geneNames.items()])
        tlocator = toyplot.locator.Explicit(np.arange(0,len(geneNames)),simpleNames[index])
        canvas,table = toyplot.matrix(D,tlocator=tlocator,llocator=tlocator,label="Molecular Cooccurence",colorshow=True)
        table.top.cells.angle = 90
        table.top.row[1].height = 75

        return canvas,table

    def plasmid_length_distribution(self,W=500,H=500,numBins=12):
        import csv
        import json
        import numpy as np
        import toyplot

        STinfo = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/ST_types.tsv"
        contigFile = self.topFolder + "carbapenamase_seq_runs/contig_meta.json"
        plasmidInfo = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/plasmid_types_detailed_v2.json"

        with open(contigFile,'r') as contigF:
            contigs = json.load(contigF)

        species = {}
        displayName = {"kpneumoniae": "Kleb.", "ecoli":"Ecoli","paeruginosa":"Pseud.","abaumannii_2":"Acin."}

        # spIndex = {"Ecoli":0,"Acin.":1,"Pseud.":2,"Kleb.":3}
        speciesColors = toyplot.color.brewer.map("Set1",count=4)

        with open(STinfo,'r') as metaData:
            rd = csv.reader(metaData,delimiter="\t")
            for n,row in enumerate(rd):
                if (n == 0):
                    species_index = row.index("species")
                    isolate_index = row.index("isolate")
                else:
                    species[ row[isolate_index] ] = row[species_index]

        ecoli_lens = []
        kleb_lens = []
        acin_lens = []
        pseu_lens = []
        with open(plasmidInfo,'r') as plasmidData:
            plasmids = json.load(plasmidData)
            for isolate,iso_plasmids in plasmids.items():
                for p in iso_plasmids.keys():
                    L = contigs[isolate][p.split(' ')[0]]["length"]
                    if (L < 1e6):
                        if ( displayName[species[isolate]] == "Ecoli" ):
                            ecoli_lens.append(L)
                        elif ( displayName[species[isolate]] == "Acin." ):
                            acin_lens.append(L)
                        elif ( displayName[species[isolate]] == "Pseud." ):
                            pseu_lens.append(L)
                        else:
                            kleb_lens.append(L)

        _,bin_edges = np.histogram(np.log10(np.array(ecoli_lens+kleb_lens+acin_lens+pseu_lens)),bins=numBins)
        bin_edges = np.array([0] + list(np.power(10,bin_edges)))
        hist_ecoli,_ = np.histogram(ecoli_lens,bins=bin_edges)
        hist_acin,_ = np.histogram(acin_lens,bins=bin_edges)
        hist_pseu,_ = np.histogram(pseu_lens,bins=bin_edges)
        hist_kleb,_ = np.histogram(kleb_lens,bins=bin_edges)

        heights = np.column_stack((hist_ecoli,hist_kleb,hist_acin,2*hist_pseu))
        colors = [speciesColors.color(0),speciesColors.color(3),speciesColors.color(1),speciesColors.color(2)]
        canvas = toyplot.Canvas(width=W,height=H)
        axes = canvas.cartesian(xscale='log',xlabel='Plasmid Size (kB)',ylabel='Number of plasmids found')
        bars = axes.bars(bin_edges[1:]/1e3,heights,color=colors)

        canvas.legend([("Ecoli/Kleb.",bars)],corner=("top-left",20,100,100))
        return canvas


    def plot_table(self, csv_file=None, df=None, sep="\t", label=None, W=200, H=200, hist_data=None, hist_row=None):
        import pandas as pd
        import toyplot

        if (csv_file is not None):
            data = pd.read_csv(csv_file,sep=sep)
        elif (df is not None):
            data = df
        else:
            raise ValueError("Need to enter table as csv file (x)OR dataframe")

        canvas = toyplot.Canvas(width=W,height=H)
        if (label is not None):
            table = canvas.table(data,label=label)
        else:
            table = canvas.table(data)
        return canvas,table

    def store_contig_info(self):
        from Bio import SeqIO
        import json
        from collections import defaultdict

        contigFile = self.topFolder + "carbapenamase_seq_runs/contig_meta.json"
        clusterFile = self.topFolder + "blaTrees/geneClusters.json"

        contigs = defaultdict(lambda: defaultdict(lambda: {}))
        with open(clusterFile,'r') as families:
            geneFamilies = json.load(families)

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            fileName,assembleName,folderName = self.assemblyName(strainName,cmpr='final')
            isolate = strainName.split('/')[-1]

            genbank = folderName + "prokka/entero.gbk"
            if (os.path.exists(genbank)):
                for contig in SeqIO.parse(genbank,"genbank"):
                    contigs[isolate][contig.name]["length"] = len(contig.seq)
                    contigs[isolate][contig.name]["blaGenes"] = set([])
                    for segment in contig.features:
                        if (segment.type == 'CDS' and 'bla' in segment.qualifiers['product'][0]):
                            dict_key = "aln_geneFamily/" + segment.qualifiers['product'][0].split("-")[0].split("_")[0]
                            if (dict_key in geneFamilies):
                                familyName = geneFamilies[dict_key][isolate][segment.qualifiers["locus_tag"][0].replace('_', " ")]
                                # print(familyName)
                                # print(geneFamilies[dict_key][isolate].keys())
                                if (len(geneFamilies[dict_key][familyName]['members']) > 1 ):
                                    contigs[isolate][contig.name]["blaGenes"].add(familyName.replace(' ', "_"))
                    contigs[isolate][contig.name]["blaGenes"] = list(contigs[isolate][contig.name]["blaGenes"])

        with open(contigFile,'w+') as contigInfo:
            json.dump(contigs,contigInfo)


    def returnECDFLine(self,ecdf):
        X = ecdf.x[1:].tolist() + ecdf.x[2:].tolist()
        X[::2] = ecdf.x[1:]
        X[1::2] = ecdf.x[2:]
        Y = ecdf.y[1:].tolist() + ecdf.y[1:-1].tolist()
        Y[::2] = ecdf.y[1:]
        Y[1::2] = ecdf.y[1:-1]

        return X,Y

    def partitionBySpecies(self,data,species,hist=False,nBins=12,binRange=None):
        from itertools import izip

        partition = defaultdict(list)
        for (datum,index) in izip(data,species):
            partition[index].append(datum)
        if (hist):
            # Get the bin edges.
            if (binRange is None):
                _,bin_edges = np.histogram(data)
            else:
                _,bin_edges = np.histogram(data,range=binRange)
            heights = []
            for index,data in partition.iteritems():
                H,_ = np.histogram(data,bins=bin_edges)
                heights.append(H)
            return np.column_stack(heights),bin_edges,partition
        else:
            return partition

    def find_repeat_regions(self):
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            _,_,folder = self.assemblyName(strainName,cmpr='final')

            polMap = {'+1':1,'-1':-1}
            nucmerFile = folder + "chrom_out.aln"
            if (os.path.exists(nucmerFile)):
                matches = defaultdict(list)
                with open(nucmerFile,'r') as aln:
                    for line in aln:
                        if ('BEGIN alignment' in line):
                            words = line.split(" ")
                            rev1 = polMap[words[4]]
                            if (rev1 == 1):
                                intervalOne = (int(words[5]),int(words[7]),rev1)
                            else:
                                intervalOne = (int(words[7]),int(words[5]),rev1)

                            rev2 = polMap[words[9]]
                            if (rev2 == 1):
                                intervalTwo = (int(words[10]),int(words[12]),polMap[words[9]])
                            else:
                                intervalTwo = (int(words[12]),int(words[10]),polMap[words[9]])

                            if (intervalOne[0] < intervalTwo[0] and intervalOne[1] < intervalTwo[1]):
                                if (intervalTwo not in matches[intervalOne]):
                                    matches[intervalOne].append(intervalTwo)
                            elif (intervalTwo[0] < intervalOne[0] and intervalTwo[1] < intervalOne[1]):
                                if (intervalOne not in matches[intervalTwo]):
                                    matches[intervalTwo].append(intervalOne)
                break
        return matches
