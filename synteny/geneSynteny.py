import glob
import utils
import numpy as np
# import pandas as pd
import csv 
from collections import defaultdict,OrderedDict 
from Bio import SeqIO
import math

import os
import json 

# def alignGeneOrder(A,B,flip=False):
#     if (flip):
#         B = B[::-1]
#     A_a,B_a,eD = align_utils.seqToSeq(A,B)
#     prc = eD/(1.*len(A_a))
#     return np.array(A_a),np.array(B_a),prc

def linkageToNewick(node,newick,parentDist,leaf_names):
    if (node.is_leaf()):
        return "%s:%.4f%s" % (leaf_names[node.id], parentDist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.4f%s" % (parentDist - node.dist, newick)
        else:
            newick = ");"
        newick = linkageToNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = linkageToNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

class wordBlock(object):

    def __init__(self,gene="KPC",globS=True,overwrite=False):
        if (globS):
            self.topFolder = glob.glob("./molecules/bla*"+gene+"*/")[0]
        else:
            self.topFolder = gene

        self.fullGeneName = self.topFolder.lstrip("./").lstrip("molecules/").split("/")[0]
        self.fullGeneName = self.fullGeneName.rstrip("_v2").rstrip("_v3").rstrip("_v4").rstrip("_v5")
        print(self.fullGeneName)

        self.clusterFile = self.topFolder + 'protein_faa/diamond_matches/allclusters_final.tsv'
        if (not os.path.exists(self.clusterFile)):
            self.clusterFile = self.topFolder + 'protein_faa/diamond_matches/allclusters.tsv'

        self.fileName = self.topFolder + "clusterString.json"
        self.matxFile = self.topFolder + "distanceMatrix.npy"
        self.synFile = self.topFolder + "syntenyCluster.newick"
        self.ed = None

        self.align = utils.seqAlign()

        if (not os.path.exists(self.fileName) or overwrite):

            with open(self.clusterFile,'r') as gcFile:
                db = csv.reader(gcFile,delimiter="\t")
                clusters = { n : np.array( [ np.array(member.split('|')) for member in row if not member == '']) for n,row in enumerate(db) } 
            self.cluster_mD = [ {'annotation':set([]),'length':[]} for n in range(len(clusters))]

            # Get list of genbank files
            plasmids = glob.glob(self.topFolder + "/input_GenBank/*.gbk")

            # Build a dictionary for quick lookup.
            plasmid_id = { x.split('/')[-1].rstrip('.gbk'):n for n,x in enumerate(plasmids) }
            self.isolateNames = np.array([ x.split('/')[-1].rstrip('.gbk') for x in plasmids ])

            geneOrder = []
            isolate_metaData = []
            for plasmid in plasmids:
                genePartition = OrderedDict()
                geneMetaData = defaultdict(lambda: defaultdict(list))
                with open(plasmid,'r') as gbkFile:
                    for contig in SeqIO.parse(gbkFile,'genbank'):
                        for feature in contig.features:
                            if (feature.type == 'CDS'):
                                geneName = feature.qualifiers['locus_tag'][0]
                                genePartition[geneName] = -1 # Initialize to no cluster ID
                                geneMetaData[geneName]['annotation'] = feature.qualifiers['product'][0]

                geneOrder.append(genePartition)
                isolate_metaData.append(geneMetaData)

            for cluster,isolates in clusters.items():
                # print(isolates)
                for isolate in isolates:
                    geneOrder[plasmid_id[isolate[0]]][isolate[1]] = cluster
                    self.cluster_mD[cluster]['annotation'].add( isolate_metaData[plasmid_id[isolate[0]]][isolate[1]]['annotation'] )
            
            self.originCluster = [ c for c in clusters.keys() if self.fullGeneName in self.cluster_mD[c]['annotation']]
            if (len(self.originCluster) > 0):
                self.originCluster = self.originCluster[0]
            else:
                self.originCluster = list(clusters.keys())[0]

            self.locusTags = [ np.array( [locus for locus in molecule.keys()] ) for molecule in geneOrder ]
            self.geneOrder = [ np.array( [cluster for cluster in molecule.values()] ) for molecule in geneOrder ]
            # self.geneOrder = [np.array(molecule) for molecule in self.geneOrder]
            self.origins = np.array([ np.where(molecule==self.originCluster)[0][0] if len(np.where(molecule==self.originCluster)[0]) > 0 else 0 for molecule in self.geneOrder])

            self.shift_and_rotate()
            self.__save__()

        else:
            self.__load__()

    def __load__(self):
        if (not os.path.exists(self.fileName)):
            raise ValueError("File never created")
        else:
            with open(self.fileName,'r') as db:
                load_dict = json.load(db)

                self.geneOrder = load_dict['geneOrder']
                self.geneOrder = [np.array(member) for member in self.geneOrder]

                self.locusTags = load_dict['locusTags']
                self.locusTags = [np.array(member) for member in self.locusTags]

                self.cluster_mD = load_dict['clusterMD']
                for cluster in self.cluster_mD:
                    cluster["annotation"] = set(cluster["annotation"])

                self.originCluster = load_dict["originCluster"]
                self.origins = np.array(load_dict["origins"])

                self.isolateNames = np.array(load_dict["isolates"])

    def __save__(self):
        cpy_cluster = self.cluster_mD.copy()
        for cluster in cpy_cluster:
            cluster['annotation'] = list(cluster['annotation'])
        
        cpy_geneOrder = self.geneOrder.copy()
        cpy_geneOrder = [member.tolist() for member in cpy_geneOrder]

        cpy_locusTags = self.locusTags.copy()
        cpy_locusTags = [member.tolist() for member in cpy_locusTags]

        save_dict = {"geneOrder":cpy_geneOrder, "clusterMD":cpy_cluster, \
                     "originCluster":self.originCluster,"origins":self.origins.tolist(), \
                     "locusTags":cpy_locusTags,"isolates":self.isolateNames.tolist()}

        with open(self.fileName,'w+') as db:
            json.dump(save_dict,db)
    
    def shift_and_rotate(self):
        # rev_motif = motif
        # index = self.originCluster
        # reflect = [ False ] ** len(self.geneOrder)
        for n,molecule in enumerate(self.geneOrder):
            xMid= int(math.floor(len(molecule)/2))
            delta = xMid - self.origins[n]
            self.geneOrder[n] = np.roll(molecule,delta,axis=0)
            self.locusTags[n] = np.roll(self.locusTags[n],delta,axis=0)

            self.origins[n] = xMid
            if (self.geneOrder[n][xMid-1:xMid+1] == [0,1,2]):
                self.geneOrder[n] = self.geneOrder[n][::-1]


    def grabGeneOrder(self,index,L=None,delta=0,flip=False):
        A = self.geneOrder[index].copy()
        if (L is not None):
            if (delta > L):
                raise ValueError("Invalid shift")

            xMin = max(0,int(math.floor(self.origins[index]-L/2))-delta)
            xMax = min(len(A),xMin + L)
            if (not flip):
                A = A[xMin:xMax] + 1
            else:
                xOrigin = self.origins[index]
                yMax = min(len(A),2*xOrigin - xMin)
                yMin = max(0, yMax-L)
                A = A[yMax:yMin:-1] + 1
        else:
            if (flip):
                A = A[::-1] + 1
            else:
                A = A + 1
        return A

    def all_to_all_circle_align(self, save=True, overwrite=False):

        matchScore = 10
        mismatchScore = -matchScore
        gapOpen = -2
        gapExtend = -1

        if (not os.path.exists(self.matxFile) or overwrite):
            self.ed = np.zeros( (len(self.geneOrder),len(self.geneOrder)) )
            for n in range(self.ed.shape[0]):
                A = self.grabGeneOrder(n)

                for nn in range(n,self.ed.shape[1]):
                    B_f = self.grabGeneOrder(nn) 
                    L_min = min(len(B_f),len(A))

                    d_f = np.zeros(L_min)
                    for m in range(L_min):
                        if (L_min == len(A)):
                            Aroll = np.roll(A,m,axis=0)
                            self.align = utils.seqAlign(A=Aroll,B=B_f)
                        else:
                            Broll = np.roll(B_f,m,axis=0)
                            self.align = utils.seqAlign(A=A,B=Broll)
                        d_f[m] = self.align.localAlign(matchScore,mismatchScore,gapOpen,gapExtend,return_align=False)

                    B_r = self.grabGeneOrder(nn,flip=True)
                    d_r = np.zeros(L_min)
                    for m in range(L_min):
                        if (L_min == len(A)):
                            Aroll = np.roll(A,m,axis=0)
                            self.align = utils.seqAlign(A=Aroll,B=B_r)
                        else:
                            Broll = np.roll(B_r,m,axis=0)
                            self.align = utils.seqAlign(A=A,B=Broll)
                        d_r[m] = self.align.localAlign(matchScore,mismatchScore,gapOpen,gapExtend,return_align=False)
                    
                    Df = np.max(d_f)
                    Dr = np.max(d_r)
                    Df /=  matchScore * L_min
                    Dr /=  matchScore * L_min

                    self.ed[n,nn] = max(Dr,Df)
                    self.ed[nn,n] = self.ed[n,nn]
                    
            if (save or overwrite):
                with open(self.matxFile,'wb') as outFile:
                    np.save(outFile,self.ed)

        else:
            with open(self.matxFile,'rb') as inFile:
                self.ed = np.load(inFile)

    def synteny_tree(self,overwrite=True):
        import scipy.cluster.hierarchy as hcl
        import scipy.spatial.distance as ssd

        if (self.ed is None):
            self.all_to_all_circle_align()
        else:
            D = ssd.squareform(1-self.ed)
            Z = hcl.linkage(D,method='average')
            tree = hcl.to_tree(Z,False)
            leaf_names = { n:name for n,name in enumerate(self.isolateNames) }
            newickString = linkageToNewick(tree,"",tree.dist,leaf_names)
            with open(self.synFile,'w+') as tree:
                tree.write(newickString)
                

    def all_to_all_align(self,L=None,delta=0,lazy_circularize=True):

        matchScore = 10
        mismatchScore = -matchScore
        gapOpen = -2
        gapExtend = -1

        self.ed = np.zeros( (len(self.geneOrder),len(self.geneOrder)) )
        for n in range(self.ed.shape[0]):
            A = self.grabGeneOrder(n,L=L,delta=delta)
            A_doub = np.concatenate([A,A],axis=0)
            # A = np.concatenate([A,A],axis=0)
            # self.align.setA(A)
            for nn in range(self.ed.shape[1]):
                B_f = self.grabGeneOrder(nn,L=L,delta=delta)
                if (len(A) > len(B_f) and lazy_circularize):
                    self.align = utils.seqAlign(A=A_doub,B=B_f)
                    # self.align.setA(A_doub)
                    # self.align.setB(B_f)
                elif (lazy_circularize):
                    B_fd = np.concatenate([B_f,B_f],axis=0)
                    self.align = utils.seqAlign(A=A,B=B_fd)
                    # self.align.setA(A)
                    # self.align.setB(B_f)
                else:
                    self.align = utils.seqAlign(A=A,B=B_f)
                    # self.align.setA(A)
                    # self.align.setB(B_f)

                d_f = self.align.localAlign(matchScore,mismatchScore,gapOpen,gapExtend,return_align=False)

                B_r = self.grabGeneOrder(nn,L=L,delta=delta,flip=True)
                if (len(A) > len(B_r) and lazy_circularize):
                    self.align = utils.seqAlign(A=A_doub,B=B_r)
                    # self.align.setA(A_doub)
                    # self.align.setB(B_r)
                elif (lazy_circularize):
                    B_rd = np.concatenate([B_r,B_r],axis=0)
                    self.align = utils.seqAlign(A=A,B=B_rd)
                    # self.align.setA(A)
                    # self.align.setB(B_r)
                else:
                    self.align = utils.seqAlign(A=A,B=B_r)
                    # self.align.setA(A)
                    # self.align.setB(B_r)

                d_r = self.align.localAlign(matchScore,mismatchScore,gapOpen,gapExtend,return_align=False)

                d_f /=  matchScore * min( len(A), len(B_f) )
                d_r /=  matchScore * min( len(A), len(B_r) )

                self.ed[n,nn] = max(d_r,d_f)
                # self.ed[nn,n] = self.ed[n,nn]

    def pairwise_align(self, iso1, iso2, L=None, delta=0, circularize=True):
        ind1 = np.where(self.isolateNames == iso1)[0][0]
        ind2 = np.where(self.isolateNames == iso2)[0][0]

        A = self.grabGeneOrder(ind1,L=L,delta=delta)
        # self.align.setA(A)

        B_f = self.grabGeneOrder(ind2,L=L,delta=delta)
        # self.align.setB(B_f)

        if (len(A) >= len(B_f) and circularize):
            Ad = np.concatenate([A,A],axis=0)
            self.align = utils.seqAlign(A=Ad,B=B_f)
            # self.align.setA(Ad)
            # self.align.setB(B_f)
        elif (circularize):
            B_fd = np.concatenate([B_f,B_f],axis=0)
            self.align = utils.seqAlign(A=A,B=B_fd)
            # self.align.setA(A)
            # self.align.setB(B_fd)
        else:
            self.align = utils.seqAlign(A=A,B=B_f)
            # self.align.setA(A)
            # self.align.setB(B_f)

        d_f,Aa_f,Ba_f = self.align.localAlign(return_align=True)
        B_r = self.grabGeneOrder(ind2,L=L,delta=delta,flip=True)
        # self.align.setB(B_r)
        if (len(A) >= len(B_r) and circularize):
            self.align = utils.seqAlign(A=Ad,B=B_r)
            # self.align.setA(Ad)
            # self.align.setB(B_r)
        elif (circularize):
            B_rd = np.concatenate([B_r,B_r],axis=0)
            self.align = utils.seqAlign(A=A,B=B_rd)
            # self.align.setA(A)
            # self.align.setB(B_rd)
        else:
            self.align = utils.seqAlign(A=A,B=B_r)
            # self.align.setA(A)
            # self.align.setB(B_r)
        d_r,Aa_r,Ba_r = self.align.localAlign(return_align=True)

        print(d_f)
        print(d_r)
        # print(",".join([str(b) for b in B_f]))
        # print(",".join([str(b) for b in B_r]))
        if (d_f > d_r):
            print(",".join([str(a) for a in A]))
            print(ind1)
            print("\n")
            print(",".join([str(b) for b in B_f]))
            print(ind2)
            return (",".join([str(a) for a in Aa_f]),",".join([str(b) for b in Ba_f]))
        else:
            print(",".join([str(a) for a in A]))
            print(ind1)
            print("\n")
            print(",".join([str(b) for b in B_r]))
            print(ind2)
            return (",".join([str(a) for a in Aa_r]),",".join([str(b) for b in Ba_r]))

    def pairwise_circle_align(self,iso1,iso2):

        matchScore = 10
        mismatchScore = -matchScore
        gapOpen = -2
        gapExtend = -1

        ind1 = np.where(self.isolateNames == iso1)[0][0]
        ind2 = np.where(self.isolateNames == iso2)[0][0]

        A = self.grabGeneOrder(ind1)
        B_f = self.grabGeneOrder(ind2) 
        L_min = min(len(B_f),len(A))
        print(L_min)

        d_f = np.zeros(L_min)
        for m in range(L_min):
            if (L_min == len(A)):
                Aroll = np.roll(A,m,axis=0)
                self.align = utils.seqAlign(A=Aroll,B=B_f)
            else:
                Broll = np.roll(B_f,m,axis=0)
                self.align = utils.seqAlign(A=A,B=Broll)
            d_f[m] = self.align.localAlign(matchScore,mismatchScore,gapOpen,gapExtend,return_align=False)

        B_r = self.grabGeneOrder(ind2,flip=True)
        d_r = np.zeros(L_min)
        for m in range(L_min):
            if (L_min == len(A)):
                Aroll = np.roll(A,m,axis=0)
                self.align = utils.seqAlign(A=Aroll,B=B_r)
            else:
                Broll = np.roll(B_r,m,axis=0)
                self.align = utils.seqAlign(A=A,B=Broll)
            d_r[m] = self.align.localAlign(matchScore,mismatchScore,gapOpen,gapExtend,return_align=False)
        
        Df = np.max(d_f)
        Dr = np.max(d_r)
        Df /=  matchScore * L_min
        Dr /=  matchScore * L_min
        
        return d_f,d_r

    def clustering_entropy(self,L=range(10,110,10),cutoff=[.25]):
        import scipy.cluster.hierarchy as hcl
        import scipy.spatial.distance as ssd
        from collections import Counter 

        entropy = np.zeros((len(L),len(cutoff)))
        for n,l in enumerate(L):
            self.all_to_all_align(L=l)
            sq = ssd.squareform(1-self.ed)
            Z = hcl.linkage(sq,method='median')
            for m,x in enumerate(cutoff):
                cluster = hcl.fcluster(Z,t=x,criterion='distance')
                numbers = Counter(cluster)
                if (m == 0):
                    total = sum(numbers.values())
                entropy[n,m] = - sum( (1.*val / total) * math.log(1.*val/total) for val in numbers.values() )

        return entropy


    def plot_ed(self):
        import pandas as pd 
        import seaborn as sns 
        import matplotlib.pylab as plt

        dat = pd.DataFrame(self.ed,index=self.isolateNames,columns=self.isolateNames)
        sns.clustermap(dat,center=.5,cmap="vlag", linewidths=.75,row_cluster=True)

        plt.show()
        # plt.savefig(family + "synteny_clusters.pdf")
        # plt.close()
