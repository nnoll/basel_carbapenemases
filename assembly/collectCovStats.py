def getCovStats(ill_pileUp,ill_txt,nan_pileUp,nan_txt):
    import json
    from collections import Counter
    # import itertools
    import numpy as np
    from scipy.special import factorial

    assembly_cov = {}
    assembly_cov['ill'] = {}
    assembly_cov['nan'] = {}

    p = 1e-3
    with open(ill_pileUp,'r') as iUp, open(ill_txt,'r') as iTxt, \
         open(nan_pileUp,'r') as nUp, open(nan_txt,'r') as nTxt:
        
        assembly_cov['ill']['frac_unmapped'] = float(iTxt.read())
        # assembly_cov[isoName]['nan']['frac_unmapped'] = float(nTxt.read()) 

        ### ILLUMINA PILEUP ###
        cov = json.load(iUp)
        uncoveredLength = sum([ sum(1 for x in contig['cov'] if np.sum(x)==0) for contig in cov ])
        totalLength = sum([ len(contig['cov']) for contig in cov ])
        # In order to serialize, you have to define a dictionary with string keys
        assembly_cov['ill']['minorX'] = []
        assembly_cov['ill']['anomPos'] = []
        for n in xrange(len(cov)):
            # minorVariantFraction = Counter([ (np.max(x),np.sum(x)) if len(x) > 0 else (0,0) for x in cov[n]['cov'] ])
            C = [ (np.max(x),np.sum(x)) if len(x) > 0 else (0,0) for x in cov[n]['cov'] ]
            P = [n for n,c in enumerate(C) if ( (p**(c[1]-c[0]) * np.exp(c[0]-c[1])) / (1.*factorial(c[1]-c[0])) ) < 1e-14 ]
            C = Counter(C)
            Cstringified = {str(k):v for k,v in C.iteritems()}
            assembly_cov['ill']['minorX'].append(Cstringified)
            assembly_cov['ill']['anomPos'].append(P)
        # assembly_cov['ill']['minorX'] = Counter(itertools.chain.from_iterable(minorVariantFraction))
        assembly_cov['ill']['uncovered'] = uncoveredLength / (1. * totalLength)

        ### NANOPORE PILEUP ###
        cov = json.load(nUp)
        nanoCov = [ [ np.sum(x) for x in contig['cov'] ] for contig in cov ]
        assembly_cov['nan']['covX'] = []
        for n in xrange(len(nanoCov)):
            assembly_cov['nan']['covX'].append( Counter(nanoCov[n]) )
        
        # if (len(nanoCov) == 1):
        #     assembly_cov['nan']['covX'] = Counter(itertools.chain.from_iterable(nanoCov))
        #     assembly_cov['nan']['covP'] = Counter([])
        # else:
        #     assembly_cov['nan']['covX'] = Counter(nanoCov[0])
        #     assembly_cov['nan']['covP'] = Counter(itertools.chain.from_iterable(nanoCov[1:]))
    return assembly_cov

if __name__== '__main__':
    import sys
    import json

    outIllDir = sys.argv[1]
    outNanoDir = sys.argv[2]
    outDir = sys.argv[3]
    species = sys.argv[4]

    print outDir

    ill_pileUp = outIllDir + "pileUp.json"
    ill_txt = outIllDir + "unmapped.txt"

    nan_pileUp = outNanoDir + "pileUp.json"
    nan_txt = outNanoDir + "unmapped.txt"

    out_file = outDir + "covStats.json"
    assembly_cov = getCovStats(ill_pileUp,ill_txt,nan_pileUp,nan_txt)
    assembly_cov['species'] = species

    with open(out_file,'w+') as out:
        json.dump(assembly_cov,out)
