import json
import utils
import glob
import matplotlib.pylab as plt
import numpy as np
from collections import defaultdict

dirName = '/scicore/home/neher/GROUP/data/Carbapenemases/'
assemblyType = 'canu_nano'

with open('/scicore/home/neher/nolln/libraries/annotationDB/card_db/aro.json','r') as db:
    cardDB = json.load(db)
cardDB = {c['accession']:c['name'] for c in cardDB}

folders = glob.glob(dirName+'/*/'+assemblyType)

strainRes = []
Ls = []
Lcontig = []
for strain in folders:

    strainName = strain
    if (strainName.endswith('/')):
        strainName = strainName[:-1]

    gbFile = glob.glob(strainName + "/prokka/*.gbf")
    if (len(gbFile) == 1):
        name = strainName.split(dirName)[1].split('/'+assemblyType)[0]
        gbFile = gbFile[0]
        resLocs,L = utils.resLocations(gbFile)
        strainRes.append(resLocs)
        Ls.append(sum(L))
        Lcontig.append(L)

nRes = []
nContigs = []
filterList = ['KPC','OXA','NDM','VIM','FOX','IMP','beta-lactamase']

geneSet = defaultdict(list)
for n in xrange(len(strainRes)):
    totNum = 0
    nContigs.append(len(strainRes[n]))
    for m in xrange(len(strainRes[n])):
        if (filterList is not None):
            for k in xrange(len(strainRes[n][m])):
                filterAgree = np.array([name in cardDB[strainRes[n][m][k][2][0]] for name in filterList])
                if (any(filterAgree)):
                    totNum += 1
                    geneSet[strainRes[n][m][k][2][0]].append(Lcontig[n][m])
        else:
            totNum += len(strainRes[n][m])
    nRes.append(totNum)


Ls = np.array(Ls)
nRes = np.array(nRes)
nContigs = np.array(nContigs)

plt.ion()
plt.figure(1)
plt.scatter(Ls, nRes)
plt.scatter(Ls[np.where(nContigs<50)], nRes[np.where(nContigs<50)])
plt.xlabel('Assembly Size')
plt.ylabel('Number of Bla Resistance Genes')

plt.figure(2)
plt.scatter(Ls, nContigs)
plt.xlabel('Assembly Size')
plt.ylabel('Number of Contigs')
