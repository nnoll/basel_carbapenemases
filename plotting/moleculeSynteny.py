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


def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

def alignGeneOrder(geneOrder,i,j,flip=False):
    Si = geneOrder[i].values()
    Sj = geneOrder[j].values()
    if (flip):
        Sj = Sj[::-1]
    Si_a,Sj_a,eD = align_utils.seqToSeq(Si,Sj)
    prc = eD/(1.*len(Si_a))
    return np.array(Si_a),np.array(Sj_a),prc

if __name__ == '__main__':
    # gene = sys.argv[1]
    topFolder = "/home/nolln/clusterData/mnt/neherGroup/Carbapenemases/blaTrees/molecules/"
    familyFolders = glob.glob(topFolder + "bla*/")

    for family in familyFolders:

        if ( os.path.exists(family + 'protein_faa/diamond_matches/allclusters.tsv')):
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
            geneNames = []
            for plasmid in plasmids:
                genePartition = OrderedDict()
                geneMetaData = defaultdict(lambda: defaultdict(list))
                with open(plasmid,'r') as gbkFile:
                    for contig in SeqIO.parse(gbkFile,'genbank'):
                        for feature in contig.features:
                            if (feature.type == 'CDS'):
                                geneName = feature.qualifiers['locus_tag'][0]
                                # geneNames.append(geneName)
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
            Z = hcl.linkage(S)
            
            tree = hcl.to_tree(Z,False)
            nwk = getNewick(tree,"",tree.dist,labels)
            with open(family + "syntenyCluster.newick",'w+') as newick:
                newick.write(nwk)

            # dat = pd.DataFrame(1-S,index=labels,columns=labels)
            # sns.clustermap(dat,center=.5,cmap="vlag", linewidths=.75,row_cluster=True)
            # plt.savefig(family + "synteny_clusters.pdf")
            # plt.close()
