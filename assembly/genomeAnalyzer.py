from __future__ import print_function
from subprocess import call
import sys
import glob
import os
import pandas as pd
import phylopandas as ph
import subprocess
import numpy as np
from collections import defaultdict
from Bio import SeqIO
# import matplotlib.pylab as plt 
import pysam

def siteEntropy(x):
    x = np.array(x)
    C = np.sum(x)
    if (C > 0):
        x /= np.sum(x)
        return -np.sum( x * np.log(x) )
    else:
        return 0

class genomeAnalyzer():

    def __init__(self,rawNanoLoc='/scicore/home/neher/GROUP/data/nanopore/',rawIlluminaLoc='/scicore/home/neher/GROUP/data/Carbapenemases_illumina/',
                 folder='/scicore/home/neher/GROUP/data/Carbapenemases/',longReadTech='nano',assembler='final',shortMapper='bwa'):

        # Intialization of all relevant data parameters!
        self.rawNanoLoc = rawNanoLoc.rstrip('/') + "/"
        self.rawIlluminaLoc = rawIlluminaLoc.rstrip('/') + "/"

        self.topFolder = folder.rstrip('/') + "/"
        self.csvFiles = [self.topFolder + "carbapenamase_seq_runs/prioritized_samples_carbapenemansen.csv", self.topFolder + "carbapenamase_seq_runs/second_sample_list.csv"]

        self.isolateFolders = glob.glob(self.topFolder+'carb???/')  
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

        self.plasmidMin = 2000
        self.plasmidMax = 2000000

        self.species = {} # Map of isolate names to their corresponding species
        self.dates = {} # Map of isolate names to their corresponding species
        harmonizeName = {'Kleb': "Kleb.", 'klepne': "Kleb.", 'klpene': "Kleb.", 'pseaer': "Pseud.","entclo":"Entclo.",
                         'Pseu': "Pseud.", 'Acin': "Acin.",'acibau': "Acin.", "Esch" : "Ecoli", "esccol": "Ecoli"}

        # Load in first file
        tbl = pd.read_csv(folder + "carbapenamase_seq_runs/prioritized_samples_carbapenemansen.csv",sep = ",")
        for index, row in tbl.iterrows():
            if (str(row['Internal #']) != "nan" and str(row['Internal #']) != " "):
                isolate = "carb%03d"%(int(row['Internal #']))
                self.species[isolate] = harmonizeName[row['Species']]

                if (row['Date'] == '20xx'):
                    self.dates[isolate] = "2010"
                else:
                    self.dates[isolate] = row["Date"]

        tbl = pd.read_csv(folder + "carbapenamase_seq_runs/second_sample_list.csv",sep = ",")
        for index, row in tbl.iterrows():
            if (str(row['Internal #']) != "nan" and str(row['Internal #']) != " "):
                isolate = "carb%03d"%(int(row['Internal #']) + 60)
                self.species[isolate] = harmonizeName[row['Species']]
                self.dates[isolate] = "20" + row["Date"].split(".")[1]

        self.genbankSpecies = {'Klebsiella pneumoniae': 'Kleb.','E.coli and Shigella':'Ecoli','Acinetobacter baumannii':'Acin.','Klebsiella oxytoca':"Kleb."}

    def updateTopFolder(self,verbose=False):
        # Link nanopore.
        offSet = [0,60]
        for n,csv in enumerate(self.csvFiles):
            table = pd.read_csv( csv, sep=",")
            for ri, row in table.iterrows():
                if (str(row.sequenced) != 'nan' and str(row["Internal #"]) != 'nan' and row.barcode != ' '):

                    nanoName = str(int(row.sequenced))
                    if (nanoName[3]=='8'):
                        nanoName = nanoName + '_Carbapenemases'
                    else:
                        nanoName = nanoName + '_USB_ESBL'
                    
                    illuminaName = row["Labor Number"].lstrip("B")

                    if (len(glob.glob( self.rawIlluminaLoc + illuminaName + '*_R1_001_val_1.fq.gz'))==0):
                        print(illuminaName)
                    else:
                        read1 = os.path.abspath(str(glob.glob( self.rawIlluminaLoc + illuminaName + '*_R1_001_val_1.fq.gz')[0]))
                        read2 = os.path.abspath(str(glob.glob( self.rawIlluminaLoc + illuminaName + '*_R2_001_val_2.fq.gz')[0]))

                        isolateNum = offSet[n] + int(row["Internal #"])   
                        isolateFolder = self.topFolder+'carb%03d/'%(isolateNum)

                        if (str(row["resequenced"]) == 'nan'): # If we did not resequence  

                            fname = self.rawNanoLoc+nanoName+'/porechop/BC%02d.fastq.gz'%int(row.barcode)
                            if (not os.path.exists(isolateFolder) and os.path.isfile(fname)):
                                os.makedirs(isolateFolder)

                            nanoTarget = isolateFolder +self.nanopore
                                                    
                            if (os.path.isfile(fname) and not os.path.exists(nanoTarget)):
                                cmd = "ln -s %s %s"%(fname, nanoTarget)
                                os.system(cmd)
                            elif (os.path.exists(nanoTarget) and verbose):
                                print(fname)
                        else: 
                            nanoName2 = str(row.resequenced)[:-2] + '_Carbapenemases'
                            fname1 =  self.rawNanoLoc+nanoName+'/porechop/BC%02d.fastq.gz'%int(row.barcode)
                            fname2 =  self.rawNanoLoc+nanoName2+'/porechop/BC%02d.fastq.gz'%int(row.rebarcode)

                            if (not os.path.exists(isolateFolder) and os.path.isfile(fname1) and os.path.isfile(fname2)):
                                os.makedirs(isolateFolder)

                            nanoTarget = isolateFolder+self.nanopore

                            # if (os.path.exists(nanoTarget)):
                            #     os.system('rm ' + nanoTarget)

                            # if (os.path.isfile(fname1) and os.path.isfile(fname2) and not os.path.exists(nanoTarget)):
                            if (os.path.isfile(fname2) and not os.path.exists(nanoTarget)):
                                # cmd = "cat %s %s > %s"%(fname1,fname2,nanoTarget)
                                cmd = "ln -s %s %s"%(fname2, nanoTarget)
                                os.system(cmd)
                            elif (not os.path.exists(nanoTarget)):
                                print(fname1)
                                print(fname2)

                        link1 = isolateFolder+self.illuminaR1
                        link2 = isolateFolder+self.illuminaR2
                        if (os.path.exists(isolateFolder) and os.path.isfile(read1) and os.path.isfile(read2) and not os.path.exists(link1) and not os.path.exists(link2)):
                            cmd = "ln -s %s %s"%(read1, link1)
                            os.system(cmd) 
                            cmd = "ln -s %s %s"%(read2, link2)
                            os.system(cmd)
                        elif (os.path.exists(link1) and os.path.exists(link1) and verbose):
                            print(link1)
                            print(link2)

    def assemble(self):
        for folder in self.isolateFolders:
            if (self.assembler == 'unicycler'):
                subFile = "/scicore/home/neher/nolln/genomeAssembly/submitScripts/unicycler_hybrid.sh"

                shortRead1 = folder + self.illuminaR1
                shortRead2 = folder + self.illuminaR2
                if (self.longReadTech == 'nano'):
                    longRead = folder + self.nanopore
                else:
                    longRead = folder +  self.pacbio
                outDir = folder + 'unicycler_' + 'hybrid_' + self.longReadTech + '/'
                sbatchCall = "sbatch --export=pair1='%s',pair2='%s',long='%s',outDir='%s' %s"%(shortRead1,shortRead2,longRead,outDir,subFile)
                contin = os.path.exists(shortRead1) and os.path.exists(shortRead2) and os.path.exists(longRead)
            elif (self.assembler == 'spades'):
                subFile = "/scicore/home/neher/nolln/genomeAssembly/submitScripts/unicycler_short.sh"

                shortRead1 = folder + self.illuminaR1
                shortRead2 = folder + self.illuminaR2
                
                outDir = folder + 'unicycler_' + 'short' + '/'
                sbatchCall = "sbatch --export=pair1='%s',pair2='%s',outDir='%s' %s"%(shortRead1,shortRead2,outDir,subFile)
                contin = os.path.exists(shortRead1) and os.path.exists(shortRead2)
            elif (self.assembler == 'canu'):
                subFile = "/scicore/home/neher/nolln/genomeAssembly/submitScripts/canu.sh"
                outDir = folder + 'canu_' + self.longReadTech + '/'

                if (self.longReadTech == 'nano'):
                    longRead = folder + self.nanopore
                else:
                    longRead = folder +  self.pacbio

                if (os.path.exists(folder + 'unicycler_hybrid_nano/assembly.fasta')):
                    genomeSize = os.stat(folder + 'unicycler_hybrid_nano/assembly.fasta').st_size
                else:
                    genomeSize = 5200000
                genomeSize = str(genomeSize/100000)
                genomeSize = genomeSize[0] + '.' + genomeSize[1] + 'M'
                sbatchCall = "sbatch --export=long='%s',genomeSize='%s',outDir='%s' %s"%(longRead,genomeSize,outDir,subFile)
                contin = os.path.exists(longRead) 
            else:
                raise ValueError("Hybrid can't be assembled directly as of now.")
            if (contin and not os.path.exists(outDir)): # We run assembly if all files exist and the assembly hasn't already been run!
                os.mkdir(outDir)
                call(sbatchCall, shell=True)
    
    def merge_assemblies(self):
        import shutil
        from interval import interval
        
        for strain in self.isolateFolders:
            print(strain)
            strainName = strain.rstrip('/')

            canuName,_,_ = self.assemblyName(strainName,cmpr="canu")
            uniName,_,_ = self.assemblyName(strainName,cmpr="unicycler")

            canuGFAName = canuName.rstrip(".fasta").replace("/pilon/","/").rstrip(".contigs.") + ".contigs.gfa"
            uniGFAName = uniName.rstrip(".fasta") + ".gfa"

            finalFolder = strain + "final/"
            finalAssembly = finalFolder + "assembly.fasta"

            if (not os.path.exists(finalFolder)):
                os.mkdir(finalFolder)

            if (os.path.exists(canuName) and not os.path.exists(uniName)):
                # shutil.copyfile(canuName, finalAssembly)
                canu_circles = set([])
                with open(canuGFAName, 'r') as input_handle:
                    contigs = {}
                    linkage = {}
                    for line in input_handle:
                        line = line.split("\t")
                        if (line[0]=="S"):
                            contigs[line[1]] = int(line[3].split(":")[2])
                            linkage[line[1]] = set([line[1]])
                        elif (line[0] == "L"):
                            linkage[line[1]].add(line[3])
                            linkage[line[3]].add(line[1])
                            if (line[1] == line[3]):
                                canu_circles.add(line[1]+"_pilon")

                with open(canuName, 'r') as input_handle, open(finalAssembly, 'w+') as output_handle:
                    for record in SeqIO.parse(input_handle,"fasta"):
                        if (record.name in canu_circles):
                            record.description += " circular=True"
                        SeqIO.write(record,output_handle,'fasta')


            elif (not os.path.exists(canuName) and os.path.exists(uniName)):
                shutil.copyfile(uniName, finalAssembly)
            elif (os.path.exists(canuName) and os.path.exists(uniName)):
                # Compute assembly entropy to see which method assembled better to use as reference.
                canuN = 0
                uniN = 0
                canu_circles = set([])
                uni_circles = set([])

                with open(canuGFAName, 'r') as input_handle:
                    contigs = {}
                    linkage = {}
                    for line in input_handle:
                        line = line.split("\t")
                        if (line[0]=="S"):
                            contigs[line[1]] = int(line[3].split(":")[2])
                            linkage[line[1]] = set([line[1]])
                        elif (line[0] == "L"):
                            linkage[line[1]].add(line[3])
                            linkage[line[3]].add(line[1])
                            if (line[1] == line[3]):
                                canu_circles.add(line[1]+"_pilon")

                    c_components = set([])
                    for val in linkage.itervalues():
                        c_components.add(frozenset(val))

                    S = 0
                    for component in c_components:
                        Lcomp = 0
                        Lcont = []
                        for contig in component:
                            Lcont.append(contigs[contig])
                            Lcomp += 1.0 * Lcont[-1]
                        Lcont = np.array(Lcont) / Lcomp
                        S += np.exp( -np.sum(Lcont * np.log(Lcont)) )
                    canuN = np.log(S)

                with open(uniGFAName, 'r') as input_handle:
                    contigs = {}
                    linkage = {}
                    for line in input_handle:
                        line = line.split("\t")
                        if (line[0]=="S"):
                            contigs[line[1]] = int(line[3].split(":")[2])
                            linkage[line[1]] = set([line[1]])
                        elif (line[0] == "L"):
                            linkage[line[1]].add(line[3])
                            linkage[line[3]].add(line[1])
                            if (line[1] == line[3]):
                                uni_circles.add(line[1])

                    c_components = set([])
                    for val in linkage.itervalues():
                        c_components.add(frozenset(val))

                    S = 0
                    for component in c_components:
                        Lcomp = 0
                        Lcont = []
                        for contig in component:
                            Lcont.append(contigs[contig])
                            Lcomp += 1.0 * Lcont[-1]
                        Lcont = np.array(Lcont) / Lcomp
                        S += np.exp(-np.sum(Lcont * np.log(Lcont)))

                    uniN = np.log(S)

                canu_contigs = []
                uni_contigs = []
                canu_len = []
                uni_len = []
                with open(canuName, 'r') as input_handle:
                    # contigLen = []
                    for record in SeqIO.parse(input_handle,"fasta"):
                        # contigLen.append(1.*len(record.seq))
                        canu_len.append(1.*len(record.seq))
                        canu_contigs.append(record.name)

                    # contigLen = np.array(contigLen)
                    # contigLen = contigLen / (1. * np.sum(contigLen))
                    # canuN.append(-np.sum( contigLen * np.log(contigLen)))

                with open(uniName, 'r') as input_handle:
                    # contigLen = []
                    for record in SeqIO.parse(input_handle,"fasta"):
                        # contigLen.append(1.*len(record.seq))
                        uni_len.append(1.*len(record.seq))
                        uni_contigs.append(record.name)

                    # contigLen = np.array(contigLen)
                    # contigLen = contigLen / (1. * np.sum(contigLen))
                    # uniN.append(-np.sum( contigLen * np.log(contigLen )))

                # canuN = np.array(canuN)
                # uniN = np.array(uniN)
                # print canuN
                # print uniN
                if (canuN <= uniN-.1):
                    ref = canuName
                    query = uniName
                    map_interval = { name:[] for name in uni_contigs }
                    q_len = { name:uni_len[n] for n,name in enumerate(uni_contigs) }
                    q_circles = uni_circles
                    r_circles = canu_circles
                else:
                    ref = uniName
                    query = canuName
                    map_interval = { name:[] for name in canu_contigs }
                    q_len = { name:canu_len[n] for n,name in enumerate(canu_contigs) }
                    q_circles = canu_circles
                    r_circles = uni_circles

                pafName = finalFolder + "aln.paf"
                minimap_call = "minimap2 -x asm5 " + ref + " " + query + " > " +  pafName
                subprocess.call(minimap_call, shell=True)
                paf_file = pd.read_csv(pafName,sep="\t",header=None)
                
                for row in paf_file.itertuples():
                    q_contig = row[1]
                    map_interval[str(q_contig)].append( interval([row[3],row[4]]) )
                
                frac_covered = {}
                for key in map_interval.iterkeys():
                    frac_covered[key] = 0
                    if (len(map_interval[key]) > 1):
                        I = map_interval[key][0]
                        for n in xrange(1,len(map_interval[key])):
                            I = I | map_interval[key][n]
                        map_interval[key] = I
                    elif (len(map_interval[key]) == 1):
                        map_interval[key] = map_interval[key][0]

                    for intval in map_interval[key]:
                        frac_covered[key] += (intval[1] - intval[0])

                    frac_covered[key] /= q_len[key]

                merged_records = []
                with open(ref, 'r') as input_handle:
                    for record in SeqIO.parse(input_handle,"fasta"):
                        if (record.name in r_circles):
                            # record.id += " circular=True"
                            record.description.replace("circular=True","")
                            record.description += " circular=True"
                        merged_records.append(record)

                with open(query, 'r') as input_handle:
                    for record in SeqIO.parse(input_handle,"fasta"):
                        if (frac_covered[record.name] < .5 and len(record.seq) > 3000 and len(record.seq) < 1e6):
                            if (record.name in q_circles):
                                # record.id += " circular=True"
                                record.description.replace("circular=True","")
                                record.description += " circular=True"
                            merged_records.append(record)

                merged_records.sort(key=lambda x: len(x.seq),reverse=True)
                with open(finalAssembly,'w+') as outF:
                    for record in merged_records:
                        SeqIO.write(record,outF,'fasta')

    def number_of_closed(self):
        ntot = 0
        nclosed = 0
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')

            name,_,_ = self.assemblyName(strainName)
            if (os.path.exists(name)):
                with open(name,'r') as fasta:
                    genomeHeader = fasta.readline()
                    if ("circular=true" not in genomeHeader.lower()):
                        nclosed += 1
                ntot += 1
        print(1 - nclosed / (1. * ntot))

    def number_of_completed(self,X=.9):
        # import bisect

        contiguity = []
        num_contigs = []
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')

            name,_,_ = self.assemblyName(strainName)
            if (os.path.exists(name)):
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

        return np.array(num_contigs), np.array(contiguity)

    def prune_merger(self):
        import json

        if (self.assembler is not "final"):
            raise ValueError("Incorrect assembly for this function")

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')

            name,assembly,folderName = self.assemblyName(strainName)
            outDir = folderName + self.shortMapper + "/pileup/"
            pileUp = outDir + "pileUp.json"

            if (os.path.exists(pileUp)): 
                mean_cov = {}
                with open(pileUp,'r') as pUP:
                    cov = json.load(pUP)
                    for n in xrange(len(cov)):
                        mean_cov[cov[n]["name"]] = cov[n]['avg_cov']

                prune_contig = set([])
                for key,val in mean_cov.iteritems():
                    # print val
                    if (val < .1):
                        prune_contig.add(key)

                # seqRecords = []
                # with open(name,'r') as fasta:
                #     for record in SeqIO.parse(fasta,'fasta'):
                #         if (record.name not in prune_contig):
                #             seqRecords.append(record)
                # newName = name.rstrip(".fasta") + "_pruned.fasta"

                # with open(newName,'w+') as newfasta:
                #     for record in seqRecords:
                #         SeqIO.write(record,newfasta,'fasta')
                if (len(prune_contig) > 0):
                    print(prune_contig)
                    print(strainName)
                    print(mean_cov)
    

    def annotate(self):
        subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/annotate.sh'
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            fileName,assembleName,folderName = self.assemblyName(strainName)
            if (os.path.exists(fileName)):
                sbatchCall = "sbatch --export=strain='%s',assembly='%s' %s"%(folderName,assembleName,subFile)
                # print sbatchCall
                call(sbatchCall, shell=True)

    def reload_bak(self):
        import shutil

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            fileName,assemblyName,folderName = self.assemblyName(strainName)
            bakName = fileName + ".bak"
            if (os.path.exists(bakName)):
                shutil.copyfile(bakName,fileName)

    def harmonize_contig_names(self):
        import fileinput
        import re, sys

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            fileName,assemblyName,folderName = self.assemblyName(strainName)
            prePilonName = fileName.replace('/pilon','')

            if (os.path.exists(fileName)):
                f = fileinput.FileInput(fileName,inplace=True,backup=".bak")
                circle = []
                with open(prePilonName,"r") as oldF:
                    for line in oldF:
                        if (line[0] == ">"):
                            circle.append("circular=true" in line.lower() or "circle=true" in line.lower())

                nContig = 0
                for line in f:
                    if (line[0] == ">"):
                        newLine = ">contig%03d"%nContig
                        if ("circular=true" in line.lower() or "circle=true" in line.lower() or circle[nContig] ):
                            newLine += " circular=True"
                        newLine += "\n"
                        nContig += 1
                    else:
                        newLine = line
                    print(newLine,end='')
                f.close()

    def harmonize_annotation_contig_names(self):
        import fileinput
        from StringIO import StringIO
        import re

        for strain in self.isolateFolders:
            stderr = StringIO()
            strainName = strain.rstrip('/')
            _,_,folderName = self.assemblyName(strainName)
            annotateName = glob.glob(folderName + '/prokka/*.gbf')
            if (len(annotateName)==1):
                annotateName = annotateName[0]
                f = fileinput.FileInput(annotateName,inplace=True,backup=".bak")
                for line in f:
                    if (line[0:5]=="LOCUS"):
                        newLine = line.replace("_pilon_pilon","_pilon\t").replace("_pilon","\t")
                        words = newLine.split()
                        # stderr.write(len(words))
                        stderr.write(newLine)
                        if (len(words) == 6):
                            stderr.write(words)
                            stderr.write(re.search('(\d+)$', words[1]).group(0) )
                    else:
                        newLine = line
                    print(newLine,end='')

                f.close()
                os.remove(annotateName+".bak")
            print(stderr.getvalue())
            stderr.close()

    def curate_annotation(self):
        # import json
        # # Read in card database
        # with open('/scicore/home/neher/nolln/libraries/annotationDB/card_db/aro.json','r') as db:
        #     cardDB = json.load(db)
        # cardDB = {c['accession']:c['description'] for c in cardDB}

        # harmonize_annotation_contig_names()
        
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            _,_,folderName = self.assemblyName(strainName)
            annotateName = glob.glob(folderName + '/prokka/*.gbf')
            if (len(annotateName)==1):
                annotateName = annotateName[0]
                curatedName = folderName + '/prokka/entero.gbk'
                with open(annotateName,'r') as gbk, open(curatedName,'w+') as new_gbk:
                    for contig in SeqIO.parse(gbk,'genbank'):
                        for feature in contig.features:
                            if (feature.type == 'CDS' and any(['ISEscan' in x for x in feature.qualifiers['inference']])):
                                feature.qualifiers['product'] = "Putative transposase (identified by ISEscan HMM)"
                                # if ('ARO' in feature.qualifiers['note'][0]):
                                #     feature.qualifiers['product'] = feature.qualifiers['note']
                                #     feature.qualifiers['note'] = "See CARD database" #(cardDB[feature.qualifiers['note'][0]]).encode('utf-8')
                        SeqIO.write(contig,new_gbk,'genbank')
                                
    def clean_contig_names(self):
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            isolateName = strainName.split('/')[-1]
            # print strainName
            # fileName,assembleName,folderName = self.assemblyName(strainName,cmpr='hybrid')
            fileName = strain + "merged/pilon/assembly.fasta"
            folderName = strain + "merged/pilon/"
            assembleName = 'assembly.fasta'
            if (os.path.exists(fileName)):
                clean_assembleName = assembleName.rstrip('.fasta') + '_clean.fasta'
                # print folderName + assembleName
                # print folderName + clean_assembleName
                with open(folderName + clean_assembleName,'w+') as outFile, open(folderName + assembleName,'r') as inFile:
                    records = SeqIO.parse(inFile,'fasta')
                    for n,record in enumerate(records):
                        record.id = 'contig' + str(n)
                        record.description = 'contig' + str(n)
                        SeqIO.write(record,outFile,'fasta')

    def assemblyName(self,strainName,cmpr=None):
        if (cmpr is None):
            assembler = self.assembler
        else: 
            assembler = cmpr

        if (assembler == 'unicycler'):
            folderName = strainName + "/unicycler_hybrid_" + self.longReadTech + "/"
            assembleName = 'assembly.fasta'
        elif (assembler == 'spades'):
            folderName = strainName + "/unicycler_short/"
            assembleName = "assembly.fasta"
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
                folderName = strainName + '/final/'
                assembleName = "assembly.fasta"
        elif (assembler == 'final'):
            if (os.path.exists(strainName + '/final/pilon/')):
                folderName = strainName + '/final/pilon/'
                assembleName = "assembly.fasta"
            elif (os.path.exists(strainName + '/final/assembly_pruned.fasta')):
                folderName = strainName + '/final/'
                assembleName = "assembly_pruned.fasta"
            else:
                folderName = strainName + '/final/'
                assembleName = "assembly.fasta"
        else:
            raise ValueError("Not a recognized assembly name")

        fileName = folderName + assembleName
        return fileName,assembleName,folderName

    def gatherPlasmids(self):
        # import shutil

        plasmidFolder = self.topFolder + 'plasmid/'
        if (not os.path.exists(plasmidFolder)):
            os.mkdir(plasmidFolder)

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            isolateName = strainName.split('/')[-1]

            _,_,folderName = self.assemblyName(strainName,cmpr='final')
            annotateName = glob.glob(folderName + "prokka/*.gbk")
            if (len(annotateName) == 1):
                fname = annotateName[0]
                with open(fname, 'r') as input_handle:
                    for nn,record in enumerate(SeqIO.parse(input_handle,"genbank")):
                        if (len(record.seq) > self.plasmidMin and len(record.seq) < self.plasmidMax):
                            record.name = isolateName + '_tig%02d'%nn
                            outName = plasmidFolder + isolateName + '_tig%02d.gbk'%nn
                            with open(outName,'w') as output_handle:
                                SeqIO.write(record,output_handle,"genbank")

        # for strain in self.isolateFolders:
        #     strainName = strain.rstrip('/')
        #     isolateName = strainName.split('/')[-1]
        #     # First try canu.
        #     fileName,assembleName,folderName = self.assemblyName(strainName,cmpr='hybrid')

            # annotateName = glob.glob(folderName + "prokka/*.gbf")
            # if (len(annotateName) == 1):
            #     fname = annotateName[0]
            #     n = 0
            #     with open(fname, 'r') as input_handle:
            #         for record in SeqIO.parse(input_handle,"genbank"):
            #             n += 1
            #     if (n < 30):
            #         m = 0
            #         with open(fname, 'r') as input_handle:
            #             for record in SeqIO.parse(input_handle,"genbank"):
            #                 # print len(record.seq)
            #                 if (len(record.seq) > 200 and m > 0):
            #                     outName = plasmidFolder + isolateName + '_pp' + record.name + '.gbk'
            #                     with open(outName,'w') as output_handle:
            #                         SeqIO.write(record,output_handle,"genbank")
            #                 m += 1
            #     else:  # If canu file is too big, then use unicycler.
            #         fileName,assembleName,folderName = self.assemblyName(strainName,cmpr='canu')
            #         n = 0
            #         with open(fname, 'r') as input_handle:
            #             for record in SeqIO.parse(input_handle,"genbank"):
            #                 n += 1
            #         if (n < 20):
            #             m = 0
            #             with open(fname, 'r') as input_handle:
            #                 for record in SeqIO.parse(input_handle,"genbank"):
            #                     # print len(record.seq)
            #                     if (len(record.seq) > 200 and m > 0):
            #                         outName = plasmidFolder + isolateName + '_pp' + record.name + '.gbk'
            #                         with open(outName,'w') as output_handle:
            #                             SeqIO.write(record,output_handle,"genbank")

    def concatenatePlasmids(self):
        from Bio.SeqRecord import SeqRecord
        # from Bio.Seq import Seq 
        # from Bio.Alphabet import HasStopCodon
        # from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
        plasmidFolder = self.topFolder + 'plasmid/translated/'

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            isolateName = strainName.split('/')[-1]

            fileName,assembleName,folderName = self.assemblyName(strainName,cmpr='hybrid')
            annotateName = glob.glob(folderName + "prokka/*.gbf")
            if (len(annotateName) == 1):
                fname = annotateName[0]
                writeName = plasmidFolder + 'carb' + strainName.lstrip(self.topFolder) + '.faa'
                # print writeName
                # concatenatedCDS = Seq("",HasStopCodon(ExtendedIUPACProtein(), '*'))
                with open(fname, 'r') as input_handle, open(writeName,'w+') as outFile:
                    for record in SeqIO.parse(input_handle,"genbank"):
                        if (len(record.seq) > self.plasmidMin and len(record.seq) < self.plasmidMax):
                            translatedCDS = record.seq.translate()
                            proteinRec = SeqRecord(translatedCDS,id=record.id,description=record.description)
                            SeqIO.write(proteinRec,outFile,'fasta')

    def makeNanoporeDB(self):
        with open(self.topFolder + 'nanoporeDB.fasta','w') as db:
            seqs = []
            for strain in self.isolateFolders:
                strainName = strain.rstrip('/')
                isolateName = strainName.split('/')[-1]
                fileName,assembleName,folderName = self.assemblyName(strainName)

                if (os.path.exists(fileName)):
                    with open(fileName, 'r') as input_handle:
                        for record in SeqIO.parse(input_handle,"fasta"):
                            record.name = isolateName + "_" + record.name
                            record.id = isolateName + "_" +  record.id
                            seqs.append(record)
            SeqIO.write(seqs,db,'fasta')

    def makeSpadesDB(self):
        with open(self.topFolder + 'spadesDB.fasta','w') as db:
            seqs = []
            for strain in self.isolateFolders:
                strainName = strain.rstrip('/')
                isolateName = strainName.split('/')[-1]
                # fileName,assembleName,folderName = self.assemblyName(strainName)
                fileName = strain + 'unicycler_short/assembly.fasta'
                if (os.path.exists(fileName)):
                    with open(fileName, 'r') as input_handle:
                        for record in SeqIO.parse(input_handle,"fasta"):
                            record.name = isolateName + "_" + record.name
                            record.id = isolateName + "_" +  record.id
                            seqs.append(record)
            SeqIO.write(seqs,db,'fasta')

    def typeCall(self):
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            _,_,folderName = self.assemblyName(strainName)
            if (os.path.exists(folderName + "prokka/entero.gbk")):
                typeCall = 'mlst ' + folderName + "prokka/entero.gbk" + " > "+ folderName + "prokka/sp_type.tsv"
                # os.system(typeCall)
                os.system(typeCall)

    def plasmid_typeCall(self):
        import shutil
        submitScript = "/scicore/home/neher/nolln/genomeAssembly/submitScripts/plasmidFinder.sh"
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            fileName,_,folderName = self.assemblyName(strainName)
            if (os.path.exists(fileName)):
                plasmidDir = folderName + "plasmidType/"
                if (os.path.exists(plasmidDir)):
                    shutil.rmtree(plasmidDir)
                os.mkdir(plasmidDir)
                print(strainName)
                subprocess.call(['bash',submitScript,fileName,plasmidDir])
                # typeCall = 'mlst ' + folderName + "prokka/entero.gbk" + " > "+ folderName + "prokka/sp_type.tsv"
                # # os.system(typeCall)
                # os.system(typeCall)
                # break

    def store_typeCall(self):
        typeCall = []
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            _,_,folderName = self.assemblyName(strainName)
            if (os.path.exists(folderName + "prokka/sp_type.tsv")):
                with open(folderName + "prokka/sp_type.tsv",'r') as STfile:
                    line = STfile.readline()
                    words = line.split('\t')
                    print(strainName)
                    print(words)
                    STinfo = {'isolate':strainName.split('/')[-1],'species':words[1],'ST':words[2],'Description':(''.join(words[3:])).replace('\n','')}
                    typeCall.append(STinfo)
        df = pd.DataFrame(typeCall)
        df.to_csv(self.topFolder + "carbapenamase_seq_runs/ST_types.tsv", sep="\t")

    def store_plasmid_typeCall(self):
        import json

        plasmid_calls = defaultdict(lambda: defaultdict(list))
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            _,_,folderName = self.assemblyName(strainName)
            plasmidFile = folderName + "plasmidType/results_tab.txt"
            if (os.path.exists(plasmidFile)):
                df = pd.read_csv(plasmidFile,sep="\t")
                name = strainName.split('/')[-1]
                for index,row in df.iterrows():
                    plasmid_calls[name][row['Contig']].append(row['Plasmid'])

        with open(self.topFolder + "carbapenamase_seq_runs/plasmid_types.json",'w+') as outFile:
            json.dump(plasmid_calls,outFile)

    def store_plasmid_typeCall_detailed(self):
        import json
        import pandas as pd
        import re

        iso_metaData = pd.read_csv(self.topFolder + "carbapenamase_seq_runs/ST_types.tsv",sep="\t")
        plasmid_calls = defaultdict(lambda: defaultdict(lambda: defaultdict(list) ) )
        for strain in self.isolateFolders:
            # Get list of plasmids for each strain
            strainName = strain.rstrip('/')
            assmble,_,folderName = self.assemblyName(strainName)

            name = strainName.split('/')[-1]
            species = iso_metaData[iso_metaData['isolate']==name]['species']
            if len(species) == 0:
                species = None
            else:
                species = list(species)[0]

            gbkFile = folderName + "prokka/entero.gbk"
            plasmidFile = folderName + "plasmidType/results_tab.txt"

            if (os.path.exists(plasmidFile) and species is not None and (species == "kpneumoniae" or species =='ecoli')):
                df = pd.read_csv(plasmidFile,sep="\t")
                for index,row in df.iterrows():
                    plasmid_calls[name][row['Contig']]['type'].append(row['Plasmid'])

                # For each identified plasmid, get the list of blaGenes, their positions, the length of contig, and any IS.
                plasmid_contigs = set(plasmid_calls[name].keys())
                for contig in SeqIO.parse(gbkFile,'genbank'):
                    contigName = contig.name

                    contin = False
                    if (contig.name in plasmid_contigs):
                        contin = True
                    elif (contig.name + ' circular=True' in plasmid_contigs):
                        contin = True
                        contigName += ' circular=True'

                    if (contin):
                        plasmid_calls[name][contigName]['bla'] = []
                        plasmid_calls[name][contigName]['blaLoc'] = []
                        plasmid_calls[name][contigName]['isLoc'] = []
                        for segment in contig.features:
                            if (segment.type == "CDS" and "bla" in segment.qualifiers['product'][0]):
                                    plasmid_calls[name][contigName]['bla'].append(segment.qualifiers['product'][0])
                                    plasmid_calls[name][contigName]['blaLoc'].append(
                                            sorted( ( int(segment.location.start),int(segment.location.end)) ) )
                            elif (segment.type == "CDS" and "ISE" in segment.qualifiers['product'][0]):
                                    plasmid_calls[name][contigName]['isLoc'].append(
                                            sorted( ( int(segment.location.start),int(segment.location.end)) ) )

            elif (species is not None and (species != "kpneumoniae" and species != "ecoli")):

                contigNames = np.array([contig.description for contig in SeqIO.parse(assmble,'fasta')])
                for n,contig in enumerate(SeqIO.parse(gbkFile,'genbank')):
                    contigName = contig.name

                    contin = False
                    if (('circular' in contigNames[n].lower()) and len(contig) > 5e3 and len(contig) < 3e5):
                        contin = True

                    if (contin):
                        plasmid_calls[name][contigName]['type']= "N/A"
                        plasmid_calls[name][contigName]['length']= len(contig)
                        plasmid_calls[name][contigName]['bla'] = []
                        plasmid_calls[name][contigName]['blaLoc'] = []
                        plasmid_calls[name][contigName]['isLoc'] = []
                        for segment in contig.features:
                            if (segment.type == "CDS" and "bla" in segment.qualifiers['product'][0]):
                                    plasmid_calls[name][contigName]['bla'].append(segment.qualifiers['product'][0])
                                    plasmid_calls[name][contigName]['blaLoc'].append(
                                            sorted( ( int(segment.location.start),int(segment.location.end)) ) )
                            elif (segment.type == "CDS" and "ISE" in segment.qualifiers['product'][0]):
                                    plasmid_calls[name][contigName]['isLoc'].append(
                                            sorted( ( int(segment.location.start),int(segment.location.end)) ) )

        with open(self.topFolder + "carbapenamase_seq_runs/plasmid_types_detailed_v2.json",'w+') as outFile:
            json.dump(plasmid_calls,outFile)

    def makeCarbGeneTrees(self,delta=0):
        import gzip
        import glob
        from Bio.SeqRecord import SeqRecord

        genBankFolder = "/scicore/home/neher/nolln/libraries/genbank_pathogens/"
        isolateGeneFolder = "/scicore/home/neher/GROUP/data/Carbapenemases/blaTrees/molecules/"

        # geneFolders = glob.glob(genBankFolder + "*/ncbi-genomes-2018-06-29/")
        geneFolders = glob.glob(genBankFolder + "kpc/ncbi-genomes-2018-06-29/")

        searchName = {"oxa23":"OXA-23","oxa48":"OXA-48","ctxM15":"CTX-M-15","ctxM24":"CTX-M-24","kpc":"KPC","ndm":"NDM","shv":"SHV","tem1d":"blaTEM","vim":"VIM" }
        geneFolderName = {"oxa23":"blaOXA-23_1_HQ700358","oxa48":"blaOXA-48_2_AY236073","ctxM15":"blaCTX-M-15_70_FJ815277",
                          "ctxM24":"blaCTX-M-24_8_EF374096","kpc":"blaKPC-3_1_HM769262","ndm":"blaNDM-4_1_JQ348841",
                          "shv":"blaSHV-11_15_DQ219473","tem1d":"blaTEM-1D_83_AF188200","vim":"blaVIM-1_1_Y18050" }

        for gene in geneFolders:
            geneName = gene.split("/")[-3]
            gbk_files = glob.glob(gene+"*.gbff.gz")
            
            gene_sequences = []
            nameCollide = defaultdict(lambda: 0)

            for gbk in gbk_files:
                isolateName = gbk.split('/')[-1]
                isolateName = isolateName.split('.')
                isolateName = isolateName[0] + "." + isolateName[1][0]

                geneLoc = None
                with gzip.open(gbk,'rb') as inGBK:
                    for contig in SeqIO.parse(inGBK,'genbank'):
                        for nS,segment in enumerate(contig.features):
                            if (segment.type == 'CDS' and searchName[geneName] in segment.qualifiers['product'][0]):
                                if (delta == 0):
                                    gene_sequences.append( SeqRecord(segment.extract(contig.seq),id=recName,name=recName ) )
                                else:
                                    geneLoc = nS
                        if (delta > 0 and geneLoc is not None):
                            recName = isolateName + "_%02d"%nameCollide[isolateName]
                            nameCollide[isolateName] += 1

                            minPos = max(geneLoc-delta,0)
                            maxPos = min(geneLoc+delta,len(contig.features))
                            seq = contig.features[minPos].extract(contig.seq)
                            for nF in xrange(minPos+1,maxPos):
                                seq += contig.features[nF].extract(contig.seq)
                            gene_sequences.append(SeqRecord(seq,id=recName,name=recName))
                            geneLoc = None

            # print(gene_sequences)
            # isoFastaFile = isolateGeneFolder + geneFolderName[geneName] + "/nc/nucleotide.fna"
            # with open(isoFastaFile,'r') as fna:
            #     for record in SeqIO.parse(fna,'fasta'):
            #         gene_sequences.append(record)
            
            if (delta > 0):
                catFastaFile = gene[:-24] + "alignments/delta%02d.fna"%delta
            else:
                catFastaFile = gene[:-24] + "gene.fna"%delta

            with open(catFastaFile,'w+') as allGenes:
                SeqIO.write(gene_sequences,allGenes,'fasta')

    def makeCarbPanGenome(self,runGene="vim",species='Kleb.'):
        import glob
        
        subFile = "/scicore/home/neher/nolln/genomeAssembly/submitScripts/panX.sh"
        genBankFolder = "/scicore/home/neher/nolln/libraries/genbank_pathogens/"
        isolateGeneFolder = "/scicore/home/neher/GROUP/data/Carbapenemases/blaTrees/molecules/"
        isolateFolder = "/scicore/home/neher/GROUP/data/Carbapenemases/"

        # geneFolders = glob.glob(genBankFolder + "*/ncbi-genomes-2018-06-29/")
        geneFolders = glob.glob(genBankFolder + runGene + "/ncbi-genomes-2018-06-29/")

        searchName = {"oxa23":"OXA-23","oxa48":"OXA-48","ctxM15":"CTX-M-15","ctxM24":"CTX-M-24","kpc":"KPC","ndm":"NDM","shv":"SHV","tem1d":"blaTEM","vim":"VIM" }
        geneFolderName = {"oxa23":"blaOXA-23_1_HQ700358","oxa48":"blaOXA-48_2_AY236073","ctxM15":"blaCTX-M-15_70_FJ815277",
                          "ctxM24":"blaCTX-M-24_8_EF374096","kpc":"blaKPC-3_1_HM769262","ndm":"blaNDM-4_1_JQ348841",
                          "shv":"blaSHV-11_15_DQ219473","tem1d":"blaTEM-1D_83_AF188200","vim":"blaVIM-1_1_Y18050" }
        
        for gene in geneFolders:
            geneName = gene.split("/")[-3]
            metaData = pd.read_csv(genBankFolder+geneName+'_genbank.csv',sep=',')
            gbk_files = glob.glob(gene+"*.gbff.gz")
            
            dirName = genBankFolder + geneName + '/panGenome/'
            # print(dirName)
            if (not os.path.exists(dirName)):
                os.mkdir(dirName)

            panX_metaData = [] #pd.DataFrame(columns=["accession","collection_date","country","region"])

            for gbk in gbk_files:
                isolateName = gbk.split('/')[-1]
                isolateName = isolateName.split('.')
                isolateName = isolateName[0] + "." + isolateName[1][0]
                isolateName = isolateName.replace("GCF","GCA")

                cmd = "ln -s " + gbk + " " + dirName + isolateName.replace('.','_') + '.gbk.gz'
                row = (metaData.index[metaData['Assembly'] == isolateName]).tolist()
                if (len(row) > 0):
                    row = row[0]
                else:
                    row = None
                if (row is not None and self.genbankSpecies[metaData.iloc[row]['#Organism Group']] == species):
                    os.system(cmd)
                    loc = metaData.iloc[row]['Location']
                    loc = loc.split(':')
                    if (len(loc)==2):
                        country = loc[0]
                        region = loc[1]
                    else:
                        country = loc[0]
                        region = None
                    date = metaData.iloc[row]["Collection Date"]
                    date = date.split('-')

                    if (len(date)==1):
                        year = date[0]
                        month = None
                    else:
                        year = date[0]
                        month = date[1]

                    accession = metaData.iloc[row]["Assembly"].replace('.','_')
                    # print(accession)
                    panX_metaData.append( {"accession": accession, "collection_year":year, "collection_month":month, "country":country, "region":region} )

            ourIsolates = set( [x.split('/')[-1].split('_')[0] for x in glob.glob( isolateGeneFolder + geneFolderName[geneName] + "/input_GenBank/*gbk" )] )
            for isolate in ourIsolates:
                gbkFile = isolateFolder + isolate + "/final/pilon/prokka/entero.gbk"
                cmd = "ln -s " + gbkFile + " " + dirName + isolate + '.gbk'
                if (self.species[isolate] == species):
                    os.system(cmd)
                    panX_metaData.append( {"accession": isolate, "collection_year":self.dates[isolate], "collection_month":None, "country":"Switzerland", "region":"Basel"} )

            panX_metaData = pd.DataFrame(panX_metaData)
            metaData_category = [ {"meta_category":"accession", "data_type":"discrete", "display":"no", "log_scale":"no"}, 
                                {"meta_category":"collection_year", "data_type":"discrete", "display":"yes", "log_scale":"no"},
                                {"meta_category":"collection_month", "data_type":"discrete", "display":"no", "log_scale":"no"},
                                {"meta_category":"country", "data_type":"discrete", "display":"yes", "log_scale":"no"},
                                {"meta_category":"region", "data_type":"discrete", "display":"no", "log_scale":"no"}
                                ]
            catInfo = pd.DataFrame(metaData_category)

            with open(dirName + "metainfo.tsv","w+") as metaInfo, open(dirName + "meta_config.tsv","w+") as metaConfig:
                panX_metaData.to_csv(metaInfo,"\t")
                catInfo.to_csv(metaConfig,"\t")


            pGF = genBankFolder + runGene + "/panGenome"
            mI = dirName + "metainfo.tsv"
            mC = dirName + "meta_config.tsv"
            name = 'entero_' + runGene
            endLoc = genBankFolder + runGene + "/panGenome"
            subCMD = "sbatch --export=panGenomeFolder='%s',metaInfo='%s',metaConfig='%s',name='%s',endLoc='%s' %s"%(pGF,mI,mC,name,endLoc,subFile)
            # print(subCMD)
            subprocess.call(subCMD, shell=True)
            # with open(isoFastaFile,'r') as fna:
            #     for record in SeqIO.parse(fna,'fasta'):
            #         gene_sequences.append(record)
            
            # if (delta > 0):
            #     catFastaFile = gene[:-24] + "alignments/delta%02d.fna"%delta
            # else:
            #     catFastaFile = gene[:-24] + "gene.fna"%delta

            # with open(catFastaFile,'w+') as allGenes:
            #     SeqIO.write(gene_sequences,allGenes,'fasta')


    def assemblyQC(self):
        import matplotlib.pylab as plt
        import mpld3
        import numpy as np
        mpld3.enable_notebook()

        canuN = []
        uniN = []
        fileSize = []
        labels = []
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            isolateName = strainName.split('/')[-1]

            canuName,_,_ = self.assemblyName(strainName,'canu')
            uniName,_,_ = self.assemblyName(strainName,'unicycler')

            # Convert to gfa file
            canuName = canuName.rstrip(".fasta").replace("/pilon/","/").rstrip(".contigs.") + ".contigs.gfa"
            uniName = uniName.rstrip(".fasta") + ".gfa"

            if (os.path.exists(canuName) and os.path.exists(uniName)):
                fileSize.append(os.path.getsize(strainName+"/"+self.nanopore)/1e7)

                with open(canuName, 'r') as input_handle:
                    contigs = {}
                    linkage = {}
                    for line in input_handle:
                        line = line.split("\t")
                        if (line[0]=="S"):
                            contigs[line[1]] = int(line[3].split(":")[2])
                            linkage[line[1]] = set([line[1]])
                        elif (line[0] == "L"):
                            linkage[line[1]].add(line[3])
                            linkage[line[3]].add(line[1])

                    c_components = set([])
                    for val in linkage.itervalues():
                        c_components.add(frozenset(val))

                    S = 0
                    for component in c_components:
                        Lcomp = 0
                        Lcont = []
                        for contig in component:
                            Lcont.append(contigs[contig])
                            Lcomp += 1.0 * Lcont[-1]
                        Lcont = np.array(Lcont) / Lcomp
                        S += np.exp( -np.sum(Lcont * np.log(Lcont)) )
                    canuN.append(np.log(S))

                with open(uniName, 'r') as input_handle:
                    contigs = {}
                    linkage = {}
                    for line in input_handle:
                        line = line.split("\t")
                        if (line[0]=="S"):
                            contigs[line[1]] = int(line[3].split(":")[2])
                            linkage[line[1]] = set([line[1]])
                        elif (line[0] == "L"):
                            linkage[line[1]].add(line[3])
                            linkage[line[3]].add(line[1])

                    c_components = set([])
                    for val in linkage.itervalues():
                        c_components.add(frozenset(val))

                    S = 0
                    for component in c_components:
                        Lcomp = 0
                        Lcont = []
                        for contig in component:
                            Lcont.append(contigs[contig])
                            Lcomp += 1.0 * Lcont[-1]
                        Lcont = np.array(Lcont) / Lcomp
                        S += np.exp(-np.sum(Lcont * np.log(Lcont)))
                    # print uniName
                    # print np.log(S)
                    uniN.append(np.log(S))
                labels.append(isolateName)

        fig,ax = plt.subplots()
        scatter = ax.scatter(np.array(canuN),np.array(uniN),fileSize)
        ax.plot(np.linspace(0,3.5),np.linspace(0,3.5),'r--')
        tooltip = mpld3.plugins.PointLabelTooltip(scatter,labels=labels)
        mpld3.plugins.connect(fig,tooltip)
        mpld3.save_html(fig,'qc.html')

        return np.array(canuN),np.array(uniN)
        
    def diamondAnnotation(self):
        subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/diamond.sh'
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            _,_,folderName = self.assemblyName(strainName)
            annotationName = glob.glob(folderName+'prokka/*.faa')
            if (len(annotationName) == 1):
                # print annotationName[0]
                sbatchCall = "sbatch --export=input='%s',output='%s' %s"%(annotationName[0],folderName+'prokka/swissProt.m8',subFile)
                call(sbatchCall, shell=True)

    def parseUniProt(self):
        from statsmodels.distributions.empirical_distribution import ECDF
        import matplotlib.pylab as plt

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
                for index,row in queryTable.iterrows():
                    if (row['name'] in blastTable.index):
                        mapLengths = blastTable.loc[row['name'],'length']
                        if (mapLengths.size > 1):
                            mapLengths = mapLengths[0]
                    # if (blastTable[blastTable['query'] == row['name']].size > 0):
                    #     topMatch = blastTable[blastTable['query'] == row['name']].iloc[0]
                        fracCov.append(1.*mapLengths/len(row['sequence']))

                fracCov = np.array(fracCov)
                ecdf = ECDF(fracCov)
                # plt.clf()
                plt.step(ecdf.x,ecdf.y)
                if (nS > 5):
                    break

            plt.xlabel('Fractional Overlap with Top UniProt Match')
            plt.ylabel('CDF')
            plt.savefig('/scicore/home/neher/nolln/uniprot_overlap.png', bbox_inches='tight')

    def polish(self):
        if (self.assembler == "canu"):
            subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/pilon.sh'
            for strain in self.isolateFolders:
                strainName = strain.rstrip('/')
                baseFolder = strain + 'canu_' + self.longReadTech + '/'
                ref = baseFolder + "assembly.contigs.fasta"

                bamFile = baseFolder + self.shortMapper + "/mappedIllumina.bam"
                sortedBamFile = baseFolder + self.shortMapper + "/mappedIllumina.sorted.bam"
                outDir = baseFolder + "pilon/"
   
                if ( not os.path.exists(outDir+'assembly.fasta') and os.path.exists(ref) and os.path.exists(bamFile)):
                    if (not os.path.exists(outDir)):
                        os.mkdir(outDir)
                    pilonCall = "sbatch --export=assembly='%s',bam='%s',sortedBam='%s',oDir='%s' %s"%(ref,bamFile,sortedBamFile,outDir,subFile)
                    call(pilonCall, shell=True)
                # else:
                    # print strainName
        elif (self.assembler == 'hybrid'):
            subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/pilon.sh'
            for strain in self.isolateFolders:
                strainName = strain.rstrip('/')
                baseFolder = strain + 'merged/'
                ref = baseFolder + "assembly.fasta"

                bamFile = baseFolder + self.shortMapper + "/mappedIllumina.bam"
                sortedBamFile = baseFolder + self.shortMapper + "/mappedIllumina.sorted.bam"
                outDir = baseFolder + "pilon/"

                if ( not os.path.exists(outDir+'assembly.fasta') and os.path.exists(ref) and os.path.exists(bamFile)):
                    if (not os.path.exists(outDir)):
                        os.mkdir(outDir)
                    pilonCall = "sbatch --export=assembly='%s',bam='%s',sortedBam='%s',oDir='%s' %s"%(ref,bamFile,sortedBamFile,outDir,subFile)
                    call(pilonCall, shell=True)
                else:
                    print(strainName)
        elif (self.assembler == 'final'):
            subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/pilon.sh'
            for strain in self.isolateFolders:
                strainName = strain.rstrip('/')
                baseFolder = strain + 'final/'
                ref = baseFolder + "assembly.fasta"

                bamFile = baseFolder + self.shortMapper + "/mappedIllumina.bam"
                sortedBamFile = baseFolder + self.shortMapper + "/mappedIllumina.sorted.bam"
                outDir = baseFolder + "pilon/"
   
                if ( not os.path.exists(outDir+'assembly.fasta') and os.path.exists(ref) and os.path.exists(bamFile)):
                    if (not os.path.exists(outDir)):
                        os.mkdir(outDir)
                    pilonCall = "sbatch --export=assembly='%s',bam='%s',sortedBam='%s',oDir='%s' %s"%(ref,bamFile,sortedBamFile,outDir,subFile)
                    call(pilonCall, shell=True)
                else:
                    print(strainName)
        else:
            raise ValueError("Canu assembly not selected!")

    def groupPanGenome(self,speciesPartition=False):
        import shutil
        if (speciesPartition):
            print("Not implemented yet!")
        else: 
            # Make the folder containing all genomes
            panGenomeFolder = self.topFolder + 'panGenome/'
            # if (os.path.exists(panGenomeFolder)):
            #     shutil.rmtree(panGenomeFolder)
            os.mkdir(panGenomeFolder)

            for strain in self.isolateFolders:
                strainName = strain.rstrip('/')
                isolateName = strainName.split('/')[-1]
                fileName,assembleName,folderName = self.assemblyName(strainName)

                annotateName = glob.glob(folderName + "prokka/*.gbf")
                if (len(annotateName) == 1):
                    fname = annotateName[0]
                    target = panGenomeFolder + isolateName + '.gbk'
                    cmd = "ln -s %s %s"%(fname, target)
                    os.system(cmd)
        
    def generateMetaMIC(self,speciesPartition=False):
        MICdata = pd.read_csv('/scicore/home/neher/GROUP/data/Carbapenemases/carbapenamase_seq_runs/MIC_data.csv') 
        if (speciesPartition):
            print("Not implemented yet!")
        else: 
            # Remove uninformative columns
            del MICdata['Labor Number']
            del MICdata['Isolate #']
            del MICdata['Species']
            del MICdata['Box']
            del MICdata['Platz']
            del MICdata['MIC (ug/mL)']
            del MICdata['Localisation']
            del MICdata['DNA concentration']
            del MICdata['Date']
            # del MICdata['Unnamed']
            del MICdata['comment']
            # print MICdata.axes
            # MICdata.rename(columns={'Unnamed: 0':'strain'},inplace=True)
            strings = MICdata.select_dtypes(['object'])
            MICdata[strings.columns] = strings.apply(lambda x: x.str.lstrip('>= ').str.lstrip('<= ') )

            MICdata.to_csv('/scicore/home/neher/GROUP/data/Carbapenemases/carbapenamase_seq_runs/metaInfo.tsv',sep='\t')

    def runPanX(self,speciesPartition=False,plasmid=False):
        if (speciesPartition):
            print("Not implemented yet!")
        else: 
            if (plasmid):
                loc = self.topFolder + 'panGenome/'
                subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/panX_plasmid.sh'
            else:
                loc = self.topFolder + 'plasmid/'
                subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/panX.sh'

            mapCall = "sbatch --export=panGenomeFolder='%s', %s"%(loc,subFile)
            subprocess.call(mapCall, shell=True)

    def storeFinalAssemblyCovStats(self):
        import os

        subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/covStats.sh'
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            isoName = strainName.split("/")[-1]
            species = self.species[isoName]

            _,_,folder = self.assemblyName(strainName,cmpr='final')

            outIllDir = folder + self.shortMapper + "/pileup/"
            outNanoDir = folder + self.longMapper + "/pileup/"

            ill_pileUp = outIllDir + "pileUp.json"
            # ill_txt = outIllDir + "unmapped.txt"

            nan_pileUp = outNanoDir + "pileUp.json"
            # nan_txt = outNanoDir + "unmapped.txt"

            if (os.path.exists(nan_pileUp) and os.path.exists(ill_pileUp) and "pilon" in folder): 
                mapCall = "sbatch --export=outIllDir='%s',outNanoDir='%s',outDir='%s',species='%s' %s"%(outIllDir,outNanoDir,folder,species,subFile)
                subprocess.call(mapCall, shell=True)
                # break

    def mapIlluminaToAssembly(self):
        if (self.shortMapper == 'bwa'):
            subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/bwa_deNovo.sh'
        else:
            subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/bowtie2_deNovo.sh'

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')

            ref,assembleName,baseFolder = self.assemblyName(strainName)

            r1File = strain + self.illuminaR1
            r2File = strain + self.illuminaR2

            if (os.path.exists(ref) and os.path.exists(r1File) and os.path.exists(r2File)):
                mapCall = "sbatch --export=ref='%s',read1='%s',read2='%s',folder='%s' %s"%(ref,r1File,r2File,baseFolder,subFile)
                subprocess.call(mapCall, shell=True)

    def mapCanuToUni(self):

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')

            refC,assembleNameC,baseFolderC = self.assemblyName(strainName,cmpr='canu')
            refU,assembleNameU,baseFolderU = self.assemblyName(strainName,cmpr='unicycler')

            r1File = strain + self.illuminaR1
            r2File = strain + self.illuminaR2

            if (os.path.exists(refC) and os.path.exists(refU)):
                # mapCall = "sbatch --export=ref='%s',read1='%s',read2='%s',folder='%s' %s"%(ref,r1File,r2File,baseFolder,subFile)
                finalFolder = strain + "final/"
                if (not os.path.exists(finalFolder)):
                    os.mkdir(finalFolder)
                mapCall = "minimap2 -x asm5 " + refU + " " + refC + " > " + finalFolder + "aln.paf"
                print(mapCall)
                # subprocess.call(mapCall, shell=True)

    def sortAndIndexBam(self):
        subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/sam_sort_index.sh'
        if (self.shortMapper == 'bwa'):
            mapF = 'bwa/'
        else:
            mapF = 'bowtie2/'

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')

            ref,_,baseFolder = self.assemblyName(strainName)
            if (os.path.exists(ref)):
                mapCall = "sbatch --export=folder='%s' %s"%(baseFolder,subFile)
                subprocess.call(mapCall, shell=True)

    def mapIlluminaToDB(self):
        if (self.shortMapper == 'bwa'):
            subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/bwa.sh'
        else:
            raise ValueError("Not implmented yet!")
            # subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/bowtie2.sh'

        folder = self.topFolder + 'nanoporeAssemblies/'
        ref = folder + 'nanoporeDB.fasta'
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            _,_,baseFolder = self.assemblyName(strainName)

            r1File = strain + self.illuminaR1
            r2File = strain + self.illuminaR2
            # print 
            if (os.path.exists(ref) and os.path.exists(r1File) and os.path.exists(r2File)):
                mapCall = "sbatch --export=ref='%s',read1='%s',read2='%s',folder='%s',name='%s' %s"%(ref,r1File,r2File,folder,strainName[-7:],subFile)
                subprocess.call(mapCall, shell=True)

    def submitKmerCompare(self,canu=True):
        if (canu):
            pre = '/scicore/home/neher/nolln/kmerCompare/ill_to_canu'
            subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/kmer_canu.sh'
        else:
            pre = '/scicore/home/neher/nolln/kmerCompare/ill_to_spades'
            subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/kmer_spades.sh'

        for strain in self.isolateFolders:
            if (os.path.exists(strain+self.illuminaR1)):
                outPrefix = pre + '_carb' + strain.lstrip(self.topFolder).rstrip('/')
                # print outPrefix
                mapCall = "sbatch --export=outPrefix='%s',isoDir='%s' %s"%(outPrefix,strain,subFile)
                subprocess.call(mapCall, shell=True)

    def allToAllAssemblyAlign(self,mode='hybrid'):
        subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/minimap2_align.sh'
        outFolder = '/scicore/home/neher/GROUP/data/Carbapenemases/assemblyMaps/'
        if (mode == 'hybrid'):
            outFolder += 'shortToLong/'
        elif (mode == 'long'):
            outFolder += 'longToLong/'
        else:
            outFolder += 'shortToShort/'

        for strain in self.isolateFolders:
            if (mode == 'hybrid' or mode == 'long'):
                refAssembly = strain + 'unicycler_hybrid_nano/assembly.fasta'#'canu_nano/assembly.contigs.fasta'
            else:
                refAssembly = strain + 'unicycler_short/assembly.fasta'
            if (os.path.exists(refAssembly)):
                mapCall = "sbatch --export=outFolder='%s',refIso='%s',canu='%s' %s"%(outFolder,strain.lstrip(self.topFolder).rstrip('/'),refAssembly,subFile)
                subprocess.call(mapCall, shell=True)

    def oneToOneAssemblyAlign(self):
        subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/minimap2_alignOne.sh'
        outFolder = '/scicore/home/neher/GROUP/data/Carbapenemases/assemblyMaps/longToLong'

        for strain in self.isolateFolders:
            uniAssembly = strain + 'unicycler_hybrid_nano/assembly.fasta'#'canu_nano/assembly.contigs.fasta'
            canuAssembly = strain + 'canu_nano/pilon/assembly.fasta'#'canu_nano/assembly.contigs.fasta'
            strainName = "/" + strain.split('/')[-2]
            if (os.path.exists(canuAssembly) and os.path.exists(uniAssembly)):
                mapCall = "sbatch --export=out='%s',uni='%s',canu='%s' %s"%(outFolder+strainName+'.sam',uniAssembly,canuAssembly,subFile)
                subprocess.call(mapCall, shell=True)

    def parseAlltoAllAlign(self,canu=True):
        outFolder = '/scicore/home/neher/GROUP/data/Carbapenemases/assemblyMaps/'
        if (canu):
            outFolder += 'shortToLong/'
        else:
            outFolder += 'shortToShort/'

        samFiles = glob.glob(outFolder+"*.sam")
        isolates = set(['carb'+s.lstrip(outFolder)[0:3] for s in samFiles])
        index = { name : n for n,name in enumerate(sorted(isolates)) }

        alignScore = np.zeros((len(index),len(index)))

        for sf in samFiles:
            mapIso = 'carb' + sf.lstrip(outFolder)[0:3]
            refIso = 'carb' + sf[-7:-4]
            # print mapIso
            # print refIso
            f = pysam.AlignmentFile(sf,"r")

            totalQLength = 0
            mapQ = 0
            for read in f:
                totalQLength += read.qlen
                if (not read.is_unmapped):
                    mapQ += read.get_overlap(read.reference_start,read.reference_end) * read.mapping_quality

            mapQ = 1. * mapQ / totalQLength
            if (mapQ > 60):
                print(mapQ)
                print(mapIso)
                print(refIso)

            alignScore[index[mapIso],index[refIso]] = mapQ
        
        np.save(outFolder+'matrix',alignScore)
        np.save(outFolder+'names',index)

    def findDominantMatch(self,illumina=True):
        folder = self.topFolder + 'nanoporeAssemblies/'
        if (illumina):
            # Get all reads in all BAM files. 
            bamFiles = glob.glob(folder+"???????.bam")
            # index = { 'carb' + bF.lstrip(folder).rstrip('.bam') : n for n,bF in enumerate(bamFiles) }

            refIsolates = set()
            with open(folder+'nanoporeDB.fasta','r') as inHandle:
                for contig in SeqIO.parse(inHandle,'fasta'):
                    refIsolates.add(contig.name[:7])

            index = { name : n for n,name in enumerate(sorted(refIsolates)) }
            # Make a map of name -> index
            numMaps = np.zeros((len(refIsolates),len(refIsolates)))

            for n,bamFile in enumerate(bamFiles):
                print(n)
                bF = pysam.AlignmentFile(bamFile,"rb")
                isolate = 'carb' + bamFile.lstrip(folder).rstrip('.bam')
                if (isolate in refIsolates):
                    r = index[isolate]
                    refNames = bF.references
                    for read in bF:
                        if (not read.is_unmapped and not read.is_secondary and not read.is_supplementary):
                            refName = refNames[read.reference_id]
                            mappedIsolate = refName[:7]
                            c = index[mappedIsolate]
                            numMaps[r,c] += 1
                else:
                    print(isolate)
                    # print refNames

            np.save(folder+'matrix',numMaps)
            np.save(folder+'names',index)
            # plt.matshow(numMaps)
            # plt.show()  
        else:
            import matplotlib.pylab as plt
            
            samFile = pysam.AlignmentFile(folder+'longOvlap.sam','r')
            refNames = samFile.references
            refIsolates = set([n[0:7] for n in refNames])
            index = { name : n for n,name in enumerate(sorted(refIsolates)) }
            # print index
            numOvlaps = np.zeros( (len(index), len(index)) )
            for read in samFile: 
                if (not read.is_unmapped):
                    qName = read.query_name.split('_')[0]
                    if (qName in index):
                        rName = refNames[read.reference_id][:7]
                        numOvlaps[index[qName],index[rName]] += 1. * read.alen / read.qlen

            # numOvlaps = numOvlaps / numOvlaps.sum(axis=1)[:,None] 
            plt.matshow(numOvlaps)
            plt.show()

    def mapNanoporeToAssembly(self):
        subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/minimap2.sh'
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            ref,assembleName,baseFolder = self.assemblyName(strainName)
            reads = strain+self.nanopore
            if (os.path.exists(ref) and os.path.exists(reads)):
                mapCall = "sbatch --export=ref='%s',reads='%s',folder='%s' %s"%(ref,reads,baseFolder,subFile)
                subprocess.call(mapCall, shell=True)

    def bamName(self,strainName,illumina):
        _,_,baseFolder = self.assemblyName(strainName)
        if (illumina):
            bamFile = baseFolder + self.shortMapper +  '/mappedIllumina.bam'
            sortFile = baseFolder + self.shortMapper +  '/mappedIllumina.sorted.bam'
        else:
            bamFile = baseFolder + self.longMapper + '/mappedNanopore.bam'
            sortFile = baseFolder + self.longMapper +  '/mappedNanopore.sorted.bam'
        return bamFile,sortFile

    def sortBAMFiles(self,illumina=True):
        subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/bam_sort_index.sh'
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            bamFile,sortFile = self.bamName(strainName,illumina)

            if (os.path.exists(bamFile) and not os.path.exists(sortFile) ): # Only submit if the pileup hasn't been made already
                sbatchCall = "sbatch --export=bamFile='%s',sortFile='%s' %s"%(bamFile,sortFile,subFile)
                call(sbatchCall, shell=True)

    def gatherErrorStats(self,illumina=True):
        # import pysam
        subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/errorCollect.sh'
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            assembleName,_,baseFolder = self.assemblyName(strainName)
            _,sortFile = self.bamName(strainName,illumina)

            if (os.path.exists(sortFile)):
                if (illumina):
                    sbatchCall = "sbatch --export=sortFile='%s',assemble='%s',folder='%s' %s"%(sortFile,assembleName,baseFolder+self.shortMapper+'/',subFile)
                else:
                    sbatchCall = "sbatch --export=sortFile='%s',assemble='%s',folder='%s' %s"%(sortFile,assembleName,baseFolder+self.longMapper+'/',subFile)
                call(sbatchCall, shell=True)
                # collectErrors(sortFile,assembleName,baseFolder)

    def parseErrors(self,errors):
        import matplotlib.pylab as plt

        errorRate = np.zeros(errors[0].keys()[-1]-1)
        predRate = np.zeros(errors[0].keys()[-1]-1)

        for k in errors[0].keys():
            if (k > 1 and len(errors[0][k]) > 0):
                totErrors = 0
                totPass = 0
                for n in errors[0][k].keys():
                    totErrors += errors[0][k][n]
                for n in errors[1][k].keys():
                    totPass += errors[1][k][n]    
                errorRate[k-2] = 1. * totErrors/(totErrors + totPass)
                predRate[k-2] = np.power(10.,(-k/10.))
                
        plt.scatter(predRate,errorRate)
        plt.xlabel('Predicted Error Rate')
        plt.ylabel('Measured Error Rate')
        plt.show()

    def getPileUp(self,illumina=True):
        subFile = '/scicore/home/neher/nolln/genomeAssembly/submitScripts/pileUp.sh'
        N = 0
        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            ref,assembleName,baseFolder = self.assemblyName(strainName)
            if (illumina):
                bamFile = baseFolder + self.shortMapper +  '/mappedIllumina.bam'
                outDir = baseFolder + self.shortMapper + "/pileup/"
                booleanVar = "1"
            else:
                bamFile = baseFolder + self.longMapper + '/mappedNanopore.sorted.bam'
                outDir = baseFolder + self.longMapper + "/pileup/"
                booleanVar = "0"

            if (os.path.exists(bamFile)): # and not os.path.exists(outDir)): # Only submit if the pileup hasn't been made already
                if (not os.path.exists(outDir)):
                    os.mkdir(outDir)
                if (not os.path.exists(outDir + "unmapped.txt") or os.stat(outDir + "unmapped.txt").st_size <= 1):
                    sbatchCall = "sbatch --export=bamFile='%s',outDir='%s',illumina='%s' %s"%(bamFile,outDir,booleanVar,subFile)
                    N += 1
                    # print(strainName)
                    # print(outDir + "unmapped.txt")
                    call(sbatchCall, shell=True)
                    # print sbatchCall
        print(N)
    def getAMR_stats(self, get_cov=False):
        # import matplotlib.pylab as plt
        import json
        # from statsmodels.distributions.empirical_distribution import ECDF

        with open('/scicore/home/neher/nolln/libraries/annotationDB/card_db/aro.json','r') as db:
            cardDB = json.load(db)
        cardDB = {c['accession']:c['description'] for c in cardDB}

        db = defaultdict(lambda: defaultdict(list))

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            _,_,folderName = self.assemblyName(strainName)

            annotateName = glob.glob(folderName + 'prokka/*.gbf')

            outDir = folderName + self.shortMapper + "/pileup/"
            pileUp = outDir + "pileUp.json"

            if (len(annotateName)==1 and ( not get_cov or os.path.exists(pileUp) )):
                annotateName = annotateName[0]
                print(annotateName)

                with open(annotateName,'r') as gbk: #with open(pileUp,'r') as pUp, 

                    if (get_cov):
                        cov = json.load(pUp)
                        contig_cov = []
                        contig_len = []

                        for n in xrange(len(cov)):
                            tmp_cov = np.array( [ sum(c) for c in cov[n]['cov'] ] )
                            frac_uncov = np.sum(tmp_cov==0) / (1. * len(tmp_cov))

                            if (frac_uncov < .1):
                                contig_cov.append(tmp_cov)
                                contig_len.append(len(tmp_cov))

                    for n,contig in enumerate(SeqIO.parse(gbk,'genbank')):
                        for feature in contig.features:
                            if (feature.type == 'CDS' and 'bla' in feature.qualifiers['product'][0]):
                                geneName = feature.qualifiers['product'][0]
                                if (get_cov):
                                    startPos = []
                                    endPos = []
                                    for interval in feature.location.parts: 
                                        startPos.append(interval.nofuzzy_start)
                                        endPos.append(interval.nofuzzy_end)
                                    minIndex = min(startPos)-500
                                    maxIndex = max(endPos)+500
                                    db[geneName]['cov'].append(contig_cov[n][minIndex:maxIndex])
                                    db[geneName]['L'].append(contig_len[n])
                                else:
                                    db[geneName]['L'].append(len(contig.seq))
        return db

    def filter_contigs(self):
        import json

        for strain in self.isolateFolders:
            strainName = strain.rstrip('/')
            ref,assembleName,baseFolder = self.assemblyName(strainName)
            outDir = baseFolder + self.shortMapper + "/pileup/"
            strainName = strain.rstrip('/')
            isolateName = strainName.split('/')[-1]

            # fileName,assembleName,folderName = self.assemblyName(strainName,cmpr='hybrid')
            fileName = strain + "merged/pilon/assembly_clean.fasta"
            folderName = strain + "merged/pilon/"
            assembleName = 'assembly_clean.fasta'
            pileUp = outDir + "pileUp.json"

            if (os.path.exists(fileName) and os.path.exists(pileUp)):
                clean_assembleName = assembleName.rstrip('_clean.fasta') + '_filtered.fasta'
                with open(folderName + clean_assembleName,'w+') as outFile, open(folderName + assembleName,'r') as inFile, open(pileUp,'r') as pUp:
                    cov = json.load(pUp)
                    frac_uncov = np.zeros(len(cov))
                    contig_len = np.zeros(len(cov))
                    for n in xrange(len(cov)):
                        contig_cov = np.array( [ sum(c) for c in cov[n]['cov'] ] )
                        contig_len[n] = len(contig_cov)
                        frac_uncov[n] = np.sum(contig_cov==0) / (1. * contig_len[n])
                    
                    chrmIndx = np.argmax(contig_len)
                    records = SeqIO.parse(inFile,'fasta')

                    contig_num = 0
                    for n,record in enumerate(records):
                        if (n == chrmIndx or frac_uncov[n] < .1):
                            record.id = 'contig' + str(contig_num)
                            record.description = 'contig' + str(contig_num)
                            SeqIO.write(record,outFile,'fasta')
                            contig_num += 1
                    


    def compile_basic_assembly_stats(self):
        import json
        import matplotlib.pylab as plt 
        from statsmodels.distributions.empirical_distribution import ECDF
        from scipy.stats import iqr

        # Statistics to record for each 
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
            print(strain)
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
                    std_cov = np.zeros(len(cov))
                    frac_uncov = np.zeros(len(cov))
                    site_entropy = np.zeros(len(cov))

                    for n in xrange(len(cov)):
                        covArr = np.array(cov[n]['cov'])
                        contig_cov = np.array( [ sum(c) for c in covArr ] )
                        site_entropy[n] = sum( [ siteEntropy(c) for c in covArr ] ) / len(covArr)

                        contig_len[n] = len(contig_cov)

                        mean_cov[n] = np.median(contig_cov)
                        std_cov[n] = iqr(contig_cov)
                        frac_uncov[n] = np.sum(contig_cov==0) / (1. * contig_len[n])

                    c_mean_cov.append( np.mean(mean_cov) )
                    c_rel_std_cov.append( np.mean(std_cov/mean_cov) )
                    c_frac_uncov.append( np.mean(frac_uncov) )
                    c_avgSiteEntropy.append( np.mean(site_entropy) )

                    ### CANU ###
                    cov = json.load(puUp)

                    contig_len = np.zeros(len(cov))
                    mean_cov = np.zeros(len(cov))
                    std_cov = np.zeros(len(cov))
                    frac_uncov = np.zeros(len(cov))
                    site_entropy = np.zeros(len(cov))

                    for n in xrange(len(cov)):
                        covArr = np.array(cov[n]['cov'])
                        contig_cov = np.array( [ sum(c) for c in covArr ] )
                        site_entropy[n] = sum( [ siteEntropy(c) for c in covArr ] ) / len(covArr)

                        contig_len[n] = len(contig_cov)

                        mean_cov[n] = np.median(contig_cov)
                        std_cov[n] = iqr(contig_cov)
                        frac_uncov[n] = np.sum(contig_cov==0) / (1. * contig_len[n])

                    u_mean_cov.append( np.mean(mean_cov) )
                    u_rel_std_cov.append( np.mean(std_cov/mean_cov) )
                    u_frac_uncov.append( np.mean(frac_uncov) )
                    u_avgSiteEntropy.append( np.mean(site_entropy) )
                break

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

        with open(self.topFolder + "assemblyCompare.npz") as outfile:
            np.savez(outfile,c_mean_cov,c_rel_std_cov,c_frac_uncov,c_frac_unmapped,c_avgSiteEntropy,
                             u_mean_cov,u_rel_std_cov,u_frac_uncov,u_frac_unmapped,u_avgSiteEntropy)

        #             mean_chrom_cov.append(mean_cov[chrm_indx])
        #             rel_std_chrom_cov.append(std_cov[chrm_indx]/mean_cov[chrm_indx])
        #             frac_chrom_uncov.append(frac_uncov[chrm_indx])

        #             mean_sub_cov.extend(mean_cov[sub_indx].tolist())
        #             rel_std_sub_cov.extend( (std_cov[sub_indx] / mean_cov[sub_indx] ).tolist())
        #             frac_sub_uncov.extend(frac_uncov[sub_indx].tolist())

        # mean_chrom_cov = np.array(mean_chrom_cov)
        # rel_std_chrom_cov = np.array(rel_std_chrom_cov )
        
        # mean_sub_cov = np.array(mean_sub_cov)
        # rel_std_sub_cov = np.array(rel_std_sub_cov)

        # frac_chrom_uncov = np.array(frac_chrom_uncov)
        # frac_sub_uncov = np.array(frac_sub_uncov)

        # frac_unmapped = np.array(frac_unmapped)

        # return mean_chrom_cov,rel_std_chrom_cov,mean_sub_cov,rel_std_sub_cov,frac_chrom_uncov,frac_sub_uncov,frac_unmapped


        # plt.figure(1)
        # ecdf1 = ECDF(mean_chrom_cov)
        # ecdf2 = ECDF(mean_sub_cov)
        # plt.step(ecdf1.x,ecdf1.y,label='Average Chromosome Coverage')
        # plt.step(ecdf2.x,ecdf2.y,label='Average SubContig Coverage')
        # plt.legend()
        # plt.savefig(self.topFolder + 'average_cov.png')

        # plt.figure(2)
        # ecdf1 = ECDF(rel_std_chrom_cov)
        # ecdf2 = ECDF(rel_std_sub_cov)
        # plt.step(ecdf1.x,ecdf1.y,label='Relative std of Chromosome Coverage')
        # plt.step(ecdf2.x,ecdf2.y,label='Relative std of SubContig Coverage')
        # plt.legend()
        # plt.savefig(self.topFolder + 'std_cov.png')

        # plt.figure(3)
        # ecdf1 = ECDF(frac_chrom_uncov)
        # ecdf2 = ECDF(frac_sub_uncov)
        # plt.step(ecdf1.x,ecdf1.y,label='Fraction Chromosome Uncovered')
        # plt.step(ecdf2.x,ecdf2.y,label='Fraction SubContig Uncovered')
        # plt.legend()
        # plt.savefig(self.topFolder + 'frac_uncov.png')
                    


