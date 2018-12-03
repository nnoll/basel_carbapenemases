def parse_genbank(gbk):
    from Bio import SeqIO

    contig_sizes = []
    with open(gbk,'r') as geneF:
        for contig in SeqIO.parse(geneF,'genbank'):
            for segment in contig.features:
                if (segment.type == 'CDS' and 'bla' in segment.qualifiers['product'][0]):
                    contig_sizes.append(len(contig.seq))
    return contig_sizes

if __name__== '__main__':
    import glob
    import numpy as np

    topFolder = '/scicore/home/neher/GROUP/data/Carbapenemases/'
    short_annotations = glob.glob( topFolder + "carb???/unicycler_short/prokka/*.gbf" )
    long_annotations = glob.glob( topFolder + "carb???/final/pilon/prokka/*.gbk" )

    long_contigs = []
    for long_gbk in long_annotations:
        print(long_gbk)
        long_contigs.extend(parse_genbank(long_gbk))

    short_contigs = []
    for short_gbk in short_annotations:
        print(short_gbk)
        short_contigs.extend(parse_genbank(short_gbk))
    
    np.savez(topFolder+"carbapenamase_seq_runs/bla_contigs.npz",np.array(short_contigs),np.array(long_contigs))