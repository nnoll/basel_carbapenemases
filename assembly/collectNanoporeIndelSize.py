import pysam
import numpy as np
import glob
import os 

def parse_nanopore_pileup(sam):
    gapSize = []
    for read in sam.fetch():
        alignment = read.get_aligned_pairs()
        if (not read.is_secondary and not read.is_supplementary and not read.is_unmapped):
            gap = 0
            for pos in alignment:

                if (pos[0] is None):
                    gap += 1
                else:
                    if (gap > 0):
                        gapSize.append(gap)
                        gap = 0
    return gapSize
    
if __name__== '__main__':
    import sys
    from collections import defaultdict
    import pickle
    import json
    import glob

    fnPU = 'final/pilon/minimap2/mappedNanopore.sam'

    # isolates = glob.glob("/scicore/home/neher/GROUP/data/Carbapenemases/carb???/")
    # for iso in isolates:

    # print(iso)
    iso = sys.argv[1]
    fn_pu = iso+fnPU

    if (os.path.exists(fn_pu)):

        # NANOPORE MAPPED TO POLISHED ASSEMBLY
        fnBAM = pysam.AlignmentFile(fn_pu,"r")
        gapSize = parse_nanopore_pileup(fnBAM)
        np.save(iso + "/errors/nanop_gap_size.npy",gapSize)



