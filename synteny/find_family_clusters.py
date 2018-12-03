import dendropy
import glob
import json
import os
import pandas as pd
import numpy as np
from collections import defaultdict,Counter
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy as hcl

# from cStringIO import StringIO

blaFamilies = glob.glob("aln_geneFamily/*newick")
geneFamily = defaultdict( lambda: defaultdict( lambda: defaultdict(list) ) )

# nIter = 0
for bF in blaFamilies:
    # dist_matrix = StringIO()
    if (os.stat(bF).st_size>0):
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

        D = squareform(D)
        L = hcl.linkage(D)
        cluster = hcl.fcluster(L,t=.15,criterion='distance')

        clusterMembers = [ [] for i in xrange(len(set(cluster))) ]
        annotations = [ [] for i in xrange(len(set(cluster))) ]

        for c,name in zip(cluster,names):
            splitName = name.split('|')
            clusterMembers[c-1].append((splitName[0],splitName[1]))
            annotations[c-1].append(splitName[2])

        familyName = [Counter(f).most_common(1)[0][0] for f in annotations]  

        # Deal with name collisions here.
        extantNames = defaultdict(int)
        for n,f in enumerate(familyName):
            if (f in extantNames):
                extantNames[f] += 1
                familyName[n] += " v" + str(extantNames[f])
            else:
                extantNames[f] = 1

        genre = bF.replace("_ft.newick","")
        for n,f in enumerate(familyName):
            geneFamily[genre][f]['members'] = clusterMembers[n]
            geneFamily[genre][f]['annotate'] = annotations[n]
        
        for n,c in enumerate(clusterMembers):
            for member in c:
                geneFamily[genre][member[0]][member[1]] = familyName[n]

with open("defined_clusters.json","w+") as cluster_dict:
    json.dump(geneFamily,cluster_dict)

        # if (nIter > 1):
        #     print annotations
        #     print familyName
        #     break
        # nIter += 1
