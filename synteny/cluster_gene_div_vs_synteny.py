import dendropy
import glob
import json
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict,Counter,OrderedDict
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy as hcl
import seaborn as sns 
import matplotlib.pylab as plt
import align_utils

topFolder = "/scicore/home/neher/GROUP/data/Carbapenemases/blaTrees/molecules/"
familyFolders = glob.glob(topFolder + "bla*/")

def alignGeneOrder(geneOrder,i,j,flip=False):
    Si = geneOrder[i].values()
    Sj = geneOrder[j].values()
    if (flip):
        Sj = Sj[::-1]
    Si_a,Sj_a,eD = align_utils.seqToSeq(Si,Sj)
    prc = eD/(1.*len(Si_a))
    return np.array(Si_a),np.array(Sj_a),prc

for family in familyFolders:
    if (not os.path.exists(family+'syn_vs_div.pdf')):
        tree_path = glob.glob(family+'nc/*bestTree*newick')
        clusterFile = family + 'protein_faa/diamond_matches/allclusters.tsv'

        if ( len(tree_path) == 0):
            tree_path = glob.glob(family + 'nc/*newick')
        if (len(tree_path) > 0):
            bF = tree_path[0]
            t = dendropy.Tree.get(path=bF,schema="newick")
            pdm = t.phylogenetic_distance_matrix()
            distances,nf = pdm._get_distance_matrix_and_normalization_factor(
                                is_weighted_edge_distances=True,
                                is_normalize_by_tree_size=False,
                                )
            taxonPairs = pdm._all_distinct_mapped_taxa_pairs

            names = set([t1.label for t1,t2 in taxonPairs]).union(set([t2.label for t1,t2 in taxonPairs]))
            names = list(names)
            name_lut = { name:n for n,name in enumerate(names) }
            DTree = np.zeros( (len(names),len(names)) )

            for t1, t2 in taxonPairs:
                n = name_lut[t1.label]
                m = name_lut[t2.label]
                DTree[n,m] = distances[t1][t2] / nf
                DTree[m,n] = DTree[n,m]
            
            avD = np.mean(DTree[:])
            S = squareform(DTree)
            LTree = hcl.linkage(S)
            
            # Align gene order. Build a dictionary for quick lookup.
            # plasmids = glob.glob(family + "/input_GenBank/*.gbk")
            gene_names = [name.split('|')[0].replace(' ', "_") for name in names]
            plasmids = [ family + "/input_GenBank/" + name + ".gbk" for name in gene_names]
            plasmid_id = { x.split('/')[-1].rstrip('.gbk'):n for n,x in enumerate(gene_names) }

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

            if ('blaTEM' in family):
                geneClusters = pd.read_csv(clusterFile,sep='\t')
            else:
                geneClusters = pd.read_csv(clusterFile,sep='\t',header=None)

            # Read in cluster file into dictionary. (Only strings. nan if empty cluster)
            cluster_dict = { row[0]:np.array( [ np.array(member.split('|')) for member in row[1:] if type(member)==str ]) \
                            for row in geneClusters.itertuples() }

            for cluster,isolates in cluster_dict.items():
                for isolate in isolates:
                    geneOrder[plasmid_id[isolate[0]]][isolate[1]] = cluster

            Dsyn = np.zeros( (len(plasmid_id),len(plasmid_id)) )
            for n in xrange(Dsyn.shape[0]):
                for nn in xrange(Dsyn.shape[1]):
                    Sn,Snn,eD = alignGeneOrder(geneOrder,n,nn)
                    SnR,SnnR,eDR = alignGeneOrder(geneOrder,n,nn,flip=True)
                    Dsyn[n,nn] = min(eD,eDR)
                    Dsyn[nn,n] = min(eD,eDR)

            S = squareform(Dsyn)
            Lsyn = hcl.linkage(S)

            # Glue distances together.
            Dsyn = avD * Dsyn
            # for n in xrange(DTree.shape[0]-1):
            #     for nn in xrange(n+1,DTree.shape[1]):
            #         DTree[n,nn] = Dsyn[n,nn]
            # S = .5*(S + S.T)

            maskU = np.zeros_like(DTree)
            maskU[np.triu_indices_from(maskU)] = True
            datT = pd.DataFrame(DTree,index=gene_names,columns=gene_names)
            datS = pd.DataFrame(DSyn,index=gene_names,columns=gene_names)

            cm = sns.clustermap(datT,cmap="vlag", linewidths=.75,mask=maskU)
            # cm.fig.showfig()
            # cm.fig.suptitle('Average Nucleotide Diversity %5.5f'%avD)
            # cm.savefig(family+'diversity.pdf')
            break