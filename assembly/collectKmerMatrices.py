if __name__== '__main__':
    import glob
    import pickle 
    from collections import defaultdict

    topFolder = '/scicore/home/neher/GROUP/data/Carbapenemases/'
    kmerMatrices = glob.glob( topFolder + "carb???/errors/nan_kmer.pkl" )

    for n,matrix in enumerate(kmerMatrices):
        print(n)
        with open(matrix,'rb') as kmer:
            if (n == 0):
                D = pickle.load(kmer)
            else:
                Dtmp = pickle.load(kmer)
                for key in Dtmp.keys():
                    if (key in D):
                        for subkey in Dtmp[key].keys():
                            if (subkey in D[key].keys()):
                                D[key][subkey] += Dtmp[key][subkey]
                            else:
                                D[key][subkey] = Dtmp[key][subkey]
                    else:
                        D[key] = Dtmp[key]

    with open(topFolder + "nan_kmer.pkl",'wb+') as final_mat:
        pickle.dump(D,final_mat)
                
