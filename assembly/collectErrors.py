def collectErrors(sortFile,assembleName,baseFolder):
    import pysam
    from Bio import SeqIO
    import numpy as np
    from collections import defaultdict

    bamFile = pysam.AlignmentFile(sortFile,'rb')
    refGenome = [np.array(list(str(seq.seq))) for seq in SeqIO.parse(assembleName,'fasta')]
    errors = [ defaultdict(lambda: defaultdict(int)), defaultdict(lambda: defaultdict(int)) ]

    for read in bamFile:
        if (not read.is_secondary and not read.is_supplementary and not read.is_unmapped):
            alignment = read.get_aligned_pairs()

            qual = np.fromstring(read.query_qualities,np.int8)
            readSeq = (read.seq).upper()
            # print readSeq
            # refSeq = read.get_reference_sequence()
            # print alignment
            for nuc in alignment:
                # if (nuc[1] is not None and nuc[1] > maxPos):
                #     maxPos = nuc[1]  
                if (nuc[0] is None):
                    errors[0][1]['-'] += 1
                elif (nuc[1] is None):
                    q = qual[nuc[0]]
                    errors[0][q][readSeq[nuc[0]]+'-'] += 1
                else:
                    refS = refGenome[read.reference_id][nuc[1]]
                    readS = readSeq[nuc[0]]
                    q = qual[nuc[0]]

                    if (refS == readS):
                        errors[1][q][refS] += 1
                    else:
                        errors[0][q][readS+refS] += 1

    for k in errors[0].keys():
        if (k not in errors[1].keys()):
            errors[1][k] = {}
    for k in errors[1].keys():
        if (k not in errors[0].keys()):
            errors[0][k] = {}

    np.save("%s/errors.npy"%baseFolder,[dict(errors[0]),dict(errors[1])])

if __name__== '__main__':
    import sys
    sortFile = sys.argv[1]
    assembleName = sys.argv[2]
    baseFolder = sys.argv[3]

    collectErrors(sortFile,assembleName,baseFolder)
