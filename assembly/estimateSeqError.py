"""
Created on Tues: January 30 2018

@author: nolln 

Purpose of this file is to either estimate the binomial error factor associated 
to sequencing (using the sequenced reference) OR estimate a p-value for a given 
variation.
"""


if __name__ == '__main__':

    import sys
    import json
    import numpy as np
    from statsmodels.distributions.empirical_distribution import ECDF
    import matplotlib.pylab as plt

    if (len(sys.argv) == 2):
        refDir = sys.argv[1]
        
        if (not refDir.endswith('/')):
            refDir += '/'

        vc = refDir + 'variants.json'
        pu = refDir + 'pileUp.json'

        with open(vc,'r') as vcf, open(pu,'r') as puf:
            pileUp = json.load(puf)[0]
            varint = json.load(vcf)[0]

            totalSamples = sum([sum(c) for c in pileUp['cov']])
            numVariants = 0
            ratio = np.zeros(len(varint.keys()))
            cov = np.zeros(len(varint.keys()))
            diff = np.zeros(len(varint.keys()))

            for nuc,pos in enumerate(varint.keys()):
                maxInd = varint[pos]['cov'].index(max(varint[pos]['cov']))
                ratio[nuc] = 1. * varint[pos]['cov'][maxInd] / np.sum(varint[pos]['cov'])
                cov[nuc] = np.sum(varint[pos]['cov'])
                diff[nuc] = cov[nuc] - varint[pos]['cov'][maxInd]
                numVariants += sum([varint[pos]['cov'][n] for n in xrange(len(varint[pos]['cov'])) 
                                if n != maxInd])

            keys = varint.keys()
            ecdf = ECDF(np.where(diff>50)[0])
            plt.step(ecdf.x, ecdf.y)
            # plt.scatter(cov, diff)
            plt.show()
            print 1. * numVariants / totalSamples
            print np.median(ratio)
    else:
        from scipy.stats import binom_test

        Dir = sys.argv[1]
        if (not Dir.endswith('/')):
            Dir += '/'

        p = sys.argv[2]
        vc = Dir + 'variants.json'
        pval = []
        with open(vc,'r') as vcf:
            varint = json.load(vcf)[0]
            for pos in varint.keys():
                localCov = sum(varint[pos]['cov'])
                numVariants = localCov - max(varint[pos]['cov'])
                pval.append(binom_test(numVariants,localCov,p))







