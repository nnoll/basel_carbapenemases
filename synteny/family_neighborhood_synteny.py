import pandas as pd
import numpy as np
import glob 
import itertools
from Bio import SeqIO
from collections import defaultdict,OrderedDict
import matplotlib.pylab as plt
import align_utils
import scipy.cluster.hierarchy as hcl
from scipy.spatial.distance import squareform
# import sys
import os
import seaborn as sns 

def alignGeneOrder(geneOrder,i,j,flip=False):
    Si = geneOrder[i].values()
    Sj = geneOrder[j].values()
    if (flip):
        Sj = Sj[::-1]
    Si_a,Sj_a,eD = align_utils.seqToSeq(Si,Sj)
    prc = eD/(1.*len(Si_a))
    return np.array(Si_a),np.array(Sj_a),prc

# gene = sys.argv[1]
topFolder = "/scicore/home/neher/GROUP/data/Carbapenemases/blaTrees/molecules/"
familyFolders = glob.glob(topFolder + "bla*/")

for family in familyFolders:

    if ( not os.path.exists(family + "synteny_clusters.pdf") and os.path.exists(family + 'protein_faa/diamond_matches/allclusters.tsv')):
        clusterFile = family + 'protein_faa/diamond_matches/allclusters.tsv'
        print family
        if ('blaTEM' in family):
            geneClusters = pd.read_csv(clusterFile,sep='\t')
        else:
            geneClusters = pd.read_csv(clusterFile,sep='\t',header=None)

        # Read in cluster file into dictionary. (Only strings. nan if empty cluster)
        cluster_dict = { row[0]:np.array( [ np.array(member.split('|')) for member in row[1:] if type(member)==str ]) for row in geneClusters.itertuples() }

        # Prepare data structure to store all annotations associated to a gene cluster.
        cluster_metaData = [ {'annotation':set([]),'length':[]} for n in xrange(len(cluster_dict))]

        # Get list of genbank files
        plasmids = glob.glob(family + "/input_GenBank/*.gbk")

        # Build a dictionary for quick lookup.
        plasmid_id = { x.split('/')[-1].rstrip('.gbk'):n for n,x in enumerate(plasmids) }

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

        for cluster,isolates in cluster_dict.iteritems():
            for isolate in isolates:
                geneOrder[plasmid_id[isolate[0]]][isolate[1]] = cluster
                cluster_metaData[cluster]['annotation'].add( isolate_metaData[plasmid_id[isolate[0]]][isolate[1]]['annotation'] )

        S = np.zeros( (len(plasmid_id),len(plasmid_id)) )
        for n in xrange(S.shape[0]):
            for nn in xrange(S.shape[1]):
                Sn,Snn,eD = alignGeneOrder(geneOrder,n,nn)
                SnR,SnnR,eDR = alignGeneOrder(geneOrder,n,nn,flip=True)
                S[n,nn] = min(eD,eDR)

        S = .5*(S + S.T)
        labels = [ x.split('/')[-1].rstrip('.gbk') for x in plasmids]

        S = squareform(S)
        Lsyn = hcl.linkage(S)
        
        # dat = pd.DataFrame(1-S,index=labels,columns=labels)
        # sns.clustermap(dat,center=.5,cmap="vlag", linewidths=.75,row_cluster=True)
        # plt.savefig(family + "synteny_clusters.pdf")
        # plt.close()
