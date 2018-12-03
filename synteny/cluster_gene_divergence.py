import dendropy
import glob
import json
import os
import pandas as pd
import numpy as np
from collections import defaultdict,Counter
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy as hcl
import seaborn as sns 
import matplotlib.pylab as plt

topFolder = "/scicore/home/neher/GROUP/data/Carbapenemases/blaTrees/molecules/"
familyFolders = glob.glob(topFolder + "bla*/")

for family in familyFolders:
    if (not os.path.exists(family+'diversity.pdf')):
        tree_path = glob.glob(family+'nc/*bestTree*newick')
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
            D = np.zeros( (len(names),len(names)) )

            for t1, t2 in taxonPairs:
                n = name_lut[t1.label]
                m = name_lut[t2.label]
                D[n,m] = distances[t1][t2] / nf
                D[m,n] = D[n,m]
            
            avD = np.mean(D[:])
            if (np.max(D[:]) > 0):
                D = 1 - D / np.max(D[:])
            names = [name.split('|')[0] for name in names]
            
            dat = pd.DataFrame(D,index=names,columns=names)
            cm = sns.clustermap(dat,cmap="vlag", linewidths=.75,row_cluster=True)
            cm.fig.suptitle('Average Nucleotide Diversity %5.5f'%avD)
            cm.savefig(family+'diversity.pdf')
