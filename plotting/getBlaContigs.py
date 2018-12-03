if __name__ == "__main__":
    from Bio import SeqIO
    import glob
    import numpy as np

    topDir = "/home/nicholas/clusterData/mnt/neherGroup/Carbapenemases/"
    GOI = ["blaOXA","blaTEM","blaKPC","blaNDM","blaVIM","blaIMP","blaCTX","blaSHV"]

    longAsmbly = glob.glob(topDir + "carb???/final/pilon/prokka/entero.gbk")
    shrtAsmbly = glob.glob(topDir + "carb???/unicycler_short/prokka/*.gbf")

    # print(len(longAsmbly))
    # print(len(shrtAsmbly))

    longContigs = []
    for la in longAsmbly:
        print(la)
        for contig in SeqIO.parse(la,'genbank'):
            for segment in contig.features:
                if (segment.type == "CDS" and 'bla' in segment.qualifiers['product'][0] and any( g in segment.qualifiers['product'][0] for g in GOI ) ):
                    longContigs.append(len(contig.seq))

    longContigs = np.array(longContigs)

    shortContigs = []
    for sa in shrtAsmbly:
        print(sa)
        for contig in SeqIO.parse(sa,'genbank'):
            for segment in contig.features:
                if (segment.type == "CDS" and 'bla' in segment.qualifiers['product'][0] and any( g in segment.qualifiers['product'][0]for g in GOI ) ):
                    shortContigs.append(len(contig.seq))

    shortContigs = np.array(shortContigs)

    np.savez("blaContigSizes.npz",longContigs,shortContigs)


