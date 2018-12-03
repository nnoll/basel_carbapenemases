import pysam
import numpy as np
import glob
import os 

def inc_matrix(mat,refID,c,alphabet):
    if (c == "." or c == ","):
        mat[refID,refID] += 1
    else:
        mat[refID,alphabet[c]] += 1

def parse_illumina_pileup(mat,bam,fa,refSeq,alphabet,refName):

    for col in bam.pileup(contig=refName,fastafile=fa):
        refID = alphabet[ refSeq[col.reference_pos] ]

        for seq in col.get_query_sequences(mark_matches=True,add_indels=True):

            if ("+" in seq): # INSERTION
                seq_split = seq.split("+")
                inc_matrix(mat,refID,seq_split[0],alphabet)
                for c in seq_split[1]:
                    if (not c.isdigit()):
                        inc_matrix(mat,4,c,alphabet)

            elif ("-" in seq): # DELETION
                seq_split = seq.split("-")
                inc_matrix(mat,refID,seq_split[0],alphabet)
                for c in seq_split[1]:
                    if (not c.isdigit()):
                        inc_matrix(mat,alphabet[c],'-',alphabet)

            else: # ALIGN
                inc_matrix(mat,refID,seq,alphabet)

        # if (index > 5e4):
        #     break

def check_quality(bam,refSeq,errors,paired=True):
    for read in bam.fetch():
        if (not read.is_secondary and not read.is_supplementary and not read.is_unmapped and (not paired or read.is_proper_pair)): 
            if (paired):
                alignment = read.get_aligned_pairs(with_seq=True)
            else:
                alignment = read.get_aligned_pairs()
                
            qual = read.query_qualities
            
            for nuc in alignment:

                if (nuc[1] is None and not paired):
                    q = qual[nuc[0]]
                    errors[q][1] += 1
                elif (not nuc[0] is None and not nuc[1] is None):
                    q = qual[nuc[0]]
                    if (paired):
                        if (not nuc[2].islower()): #refS == readS):
                            errors[q][0] += 1
                        else:
                            errors[q][1] += 1
                    else:
                        refS = refSeq[nuc[1]]
                        readS = read.query_sequence[nuc[0]].upper()
                        if (refS == readS):
                            errors[q][0] += 1
                        else:
                            errors[q][1] += 1

#     for read in bam.fetch():
    
#         if (not read.is_secondary and not read.is_supplementary and not read.is_unmapped):
        
#             alignment = read.get_aligned_pairs()
#             qual = read.query_qualities
# #             readSeq = read.get_forward_sequence().upper()
        
#             for nuc in alignment:
            
#                 # if (nuc[0] is None):
#                 #     errors[1][1] += 1
                
#                 # elif (nuc[1] is None):
#                 #     q = qual[nuc[0]]
#                 #     errors[q][1] += 1
#                 if (not nuc[0] is None and not nuc[1] is None):
#                     refS = refSeq[nuc[1]]
#                     readS = read.query_sequence[nuc[0]]
#                     q = qual[nuc[0]]
                
#                     if (refS == readS):
#                         errors[q][0] += 1
#                     else:
#                         errors[q][1] += 1
                        
def illumina_insert_size(bam):
    return np.array( [ abs(read.template_length) for read in bam.fetch() ] )
              
def parse_nanopore_pileup(mat,sam,refSeq,alphabet):

    for read in sam.fetch():
        alignment = read.get_aligned_pairs()
        if (not read.is_secondary and not read.is_supplementary and not read.is_unmapped and not read.is_reverse):
#             print(read.qname)
#             print(read.get_forward_sequence())
            for pos in alignment:

                if (pos[1] is None):
                    refInd = 4
                else:
                    refInd = alphabet[refSeq[pos[1]]]

                if (pos[0] is None):
                    qryInd = 4
                else:
                    qryInd = alphabet[read.query_sequence[pos[0]]]

                mat[refInd,qryInd] += 1
            
def nanopore_kmer_pileup(D,sam,refSeq,k=5):

    for read in sam.fetch():
        if (not read.is_secondary and not read.is_supplementary and not read.is_unmapped):
            alignment = read.get_aligned_pairs()

            refStr = ''
            qryStr = ''
            for pos in alignment:

                if (pos[1] is None):
                    refStr += "-"
                else:
                    refStr += refSeq[pos[1]]

                if (pos[0] is None):
                    qryStr += "-"
                else:
                    qryStr += read.query_sequence[pos[0]]
            
            for locus in range( len(qryStr)-k ):
                D[ refStr[locus:locus+k] ][ qryStr[locus:locus+k]  ] += 1
        

if __name__== '__main__':
    import sys
    from collections import defaultdict
    import pickle
    import json

    iAssemble = 'canu_nano/assembly.contigs.fasta'
    iiPU = 'canu_nano/bwa/mappedIllumina.sorted.bam'
    inPU = 'canu_nano/minimap2/mappedNanopore.sam'

    fAssemble = 'final/pilon/assembly.fasta'
    fiPU = 'final/pilon/bwa/mappedIllumina.sorted.bam'
    fnPU = 'final/pilon/minimap2/mappedNanopore.sam'

    alphabet = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4, "*":4,
                "a": 0, "c": 1, "g": 2, "t": 3, "N": 4, "n":4}

    iso = sys.argv[1]

    fn_pu = iso+fnPU
    in_pu = iso+fiPU

    fi_pu = iso+fiPU
    ii_pu = iso+iiPU
    
    ia = iso+iAssemble
    fa = iso+fAssemble

    if (os.path.exists(ia) and os.path.exists(fa) and os.path.exists(fn_pu) and os.path.exists(in_pu) and os.path.exists(fi_pu) and os.path.exists(ii_pu)):

        iFA = pysam.FastaFile(ia)
        fFA = pysam.FastaFile(fa)

        iRef = iFA.references
        print(iRef)
        iRefSeq = "".join( iFA.fetch(reference=ref) for ref in iRef )
        iRefSeq = iRefSeq.upper()

        fRef = fFA.references
        fRefSeq = "".join( fFA.fetch(reference=ref) for ref in fRef )
        fRefSeq = fRefSeq.upper()

        # ILLUMINA MAPPED TO UNPOLISHED ASSEMBLY
        # iiBAM = pysam.AlignmentFile(ii_pu,"rb")
        # print(iiBAM.references)
        # iErrors = np.zeros((5,5))
        # parse_illumina_pileup(iErrors,iiBAM,iFA,iRefSeq,alphabet,iRef[0])
        # np.save(iso + "/errors/ilmna_to_unpolish.npy",iErrors)

        # NANOPORE MAPPED TO UNPOLISHED ASSEMBLY
        inBAM = pysam.AlignmentFile(in_pu,"r")
        iErrors = np.zeros((5,5))
        parse_nanopore_pileup(iErrors,inBAM,iRefSeq,alphabet)
        np.save(iso + "/errors/nanop_to_unpolish.npy",iErrors)

        # ILLUMINA MAPPED TO POLISHED ASSEMBLY
        # fiBAM = pysam.AlignmentFile(fi_pu,"rb")
        # fErrors = np.zeros((5,5))
        # parse_illumina_pileup(fErrors,fiBAM,fFA,fRefSeq,alphabet,fRef[0])
        # np.save(iso + "/errors/ilmna_to_polish.npy",fErrors)

        # ILLUMINA INSERT SIZES
        # fiBAM = pysam.AlignmentFile(fi_pu,"rb")
        # inserts = illumina_insert_size(fiBAM)
        # np.save(iso + "/errors/ilmna_inserts.npy",inserts)

        # ILLUMINA QUALITY SCORES
        # fiBAM = pysam.AlignmentFile(fi_pu,"rb")
        # errors = defaultdict(lambda: np.zeros(2))
        # check_quality(fiBAM,fRefSeq,errors,paired=True)
        # with open(iso + "/errors/ill_errors.pkl" ,'wb+') as ill_errs:
        #     pickle.dump(dict(errors),ill_errs,pickle.HIGHEST_PROTOCOL)

        # NANOPORE MAPPED TO POLISHED ASSEMBLY
        fnBAM = pysam.AlignmentFile(fn_pu,"r")
        fErrors = np.zeros((5,5))
        parse_nanopore_pileup(fErrors,fnBAM,fRefSeq,alphabet)
        np.save(iso + "/errors/nanop_to_polish.npy",fErrors)

        # NANOPORE KMER DISTRIBTUION
        # fnBAM = pysam.AlignmentFile(fn_pu,"r")
        # D = defaultdict(lambda: defaultdict(lambda: 0) )
        # nanopore_kmer_pileup(D,fnBAM,fRefSeq,k=6)
        # with open(iso + "/errors/nan_kmer.pkl" ,'wb+') as nan_kmer:
        #     tmp_dict = json.loads(json.dumps(D))
        #     # json.dump(tmp_dict,nan_kmer)
        #     pickle.dump(tmp_dict,nan_kmer,pickle.HIGHEST_PROTOCOL)

        # NANOPORE QUALITY SCORES
        # fnBAM = pysam.AlignmentFile(fn_pu,"r")
        # errors = defaultdict(lambda: np.zeros(2))
        # check_quality(fnBAM,fRefSeq,errors,paired=False)
        # with open(iso + "/errors/nan_errors.pkl" ,'wb+') as nan_errs:
        #     tmp_dict = dict(errors)
        #     # json.dump(tmp_dict,nan_errs)
        #     pickle.dump(tmp_dict,nan_errs,pickle.HIGHEST_PROTOCOL)



