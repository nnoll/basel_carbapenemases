import scipy.special as ss
import numpy as np
from collections import defaultdict,Counter 

def beta_binom(n,alpha,beta):
    F = np.ones(n+1)
    for k in range(n+1):
        F[k] = ss.comb(n,k)*( 1.*ss.beta(k+alpha,n-k+beta) / ss.beta(alpha,beta) )
    return F

def est_beta_params(m1,m2,n):
    alpha = (n*m1 - m2) / (n *(m2/m1 - m1 -1) + m1)
    beta = (n - m1)*(n-m2/m1) / (n *(m2/m1 - m1 -1) + m1)
    return alpha,beta

def compute_moments(p):
    m1 = sum( i*p[i] for i in range(len(p)))
    m2 = sum( (i**2)*p[i] for i in range(len(p)))
    return m1,m2

def compute_score(w, entropy=False):
    if (entropy):
        return -sum( (v/(1.*len(w))) * np.log(v/(1.*len(w))) for v in Counter(w).values() )
    else:
        score = 0
#         for n in range(1,len(w)-1):
#             score += (1./n) * sum( w[m]==w[m+n] for m in range(len(w)-n) )
        blockLen=1
        for n,c in enumerate(w):
            if (n > 0):
                if (w[n-1] == c):
                    blockLen+=1
                else:
                    score += np.exp((blockLen - 1))
                    blockLen=1
                    
        score += np.exp((blockLen - 1))
        score /= np.exp( len(w) - 1)
        return score
    
def compute_total_kmer_plot(Dd,n=5):
    binned_err = defaultdict(lambda: 0)

    for k in Dd.keys():
        if "-" not in k:
            for kk in Dd[k].keys():
                mask = ""
                for m in range(n):
                    if (k[m] == kk[m]):
                        mask += "a"
                    elif (kk[m] == "-"):
                        mask += "-"
                    else:
                        mask += "*"
                binned_err[mask] += Dd[k][kk]

    errs = [ (k,binned_err[k]) for k in sorted(binned_err.keys())]
    mask,errs = zip(*errs)
    return np.array(errs),np.array(mask)


if __name__== '__main__':
    import pickle 

    db = '/scicore/home/neher/GROUP/data/Carbapenemases/nan_kmer.pkl'
    with open(db,'rb') as data:
        D = pickle.load(data)

    print( len(D) )
    print( np.mean( [len(v) for v in D.values()] ))

    saveDir = "/scicore/home/neher/GROUP/data/Carbapenemases/carbapenamase_seq_runs/"
    kmer_errors = np.array([ [ sum(D[k][kk] for kk in D[k].keys() if k != kk), sum(D[k][kk] for kk in D[k].keys() if k == kk) ] for k in D.keys() ])
    kmers = np.array(list(D.keys()))

    np.savez(saveDir+"kmer_errors.npy",kmer_errors,kmers)

    bin_dist = []
    alpha = []
    beta = []

    for k in kmers:
        ind_dist = np.zeros(7)
        if ("-" not in k):
            for kk in D[k].keys():
                str_dist = sum(c!=cc for c,cc in zip(k,kk) )
                ind_dist[str_dist] += D[k][kk]
        if (np.sum(ind_dist) > 0):
            ind_dist /= np.sum(ind_dist)
            m1,m2 = compute_moments(ind_dist)
            a,b = est_beta_params(m1,m2,6)
        
            alpha.append(a)
            beta.append(b)
            bin_dist.append(ind_dist)

    alpha = np.array(alpha)
    beta = np.array(beta)
    bin_dist = np.array(bin_dist)
    np.savez(saveDir+"kmer_parameters.npy",alpha,beta,bin_dist)

    all_error_dist = [ (kmers[n],x[0]/(1.*sum(x))) for n,x in enumerate(kmer_errors) if "-" not in kmers[n] ]
    eKmer,errD = zip(*all_error_dist)
    eKmer = np.array(eKmer)
    errD = np.array(errD)
    sKmer = np.array([ compute_score(k,entropy=False) for k in eKmer ]) 

    np.savez(saveDir+"selected_kmers.npy",eKmer,errD,sKmer)

    BEs,masks = compute_total_kmer_plot(D,n=6)
    BEs = BEs / np.sum(BEs)

    BEs = np.array(BEs[:-1])
    masks = np.array(masks[:-1])

    gD = np.zeros((7,7))

    deltaI = np.zeros(7)
    deltaS = np.zeros(7)
    for n,m in enumerate(masks):
        numI = sum( c == "-" for c in m )
        numS = sum( c == "*" for c in m )
        
        gD[numI,numS] += BEs[n]
        if (numI == 2):
            distI = np.abs(np.diff( [i for i,c in enumerate(m) if c == "-"]))
            deltaI[distI] += BEs[n]
            
        if (numS == 2):
            distS = np.abs(np.diff( [i for i,c in enumerate(m) if c == "*"]))
            deltaS[distS] += BEs[n]
        

    deltaI = deltaI / np.sum(deltaI)
    deltaS = deltaS / np.sum(deltaS)

    np.savez(saveDir+"error_types.npy",BEs,masks,gD,deltaI,deltaS)


