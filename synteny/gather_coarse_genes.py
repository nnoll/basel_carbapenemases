from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import defaultdict
import glob
import numpy as np
from itertools import chain,combinations
import os

def all_subsets(ss):
	return chain(*map(lambda x: combinations(ss,x), range(0,len(ss)+1)))

topDir = '/scicore/home/neher/GROUP/data/Carbapenemases/'
gbkWC = 'carb???/final/pilon/prokka/*gbk'

gbkFiles = glob.glob(topDir+gbkWC)
bla_genes = defaultdict(lambda: defaultdict(list))

for gbk in gbkFiles:
	isolate_name = "carb" + gbk.split('carb')[1].split('/')[0]
	print isolate_name
	with open(gbk,'r') as geneF:
		for contig in SeqIO.parse(geneF,'genbank'):
			for segment in contig.features:
				if (segment.type == 'CDS' and 'bla' in segment.qualifiers['product'][0]):
					prod = segment.qualifiers['product'][0]
					key = prod.split('-')[0].split('_')[0]
					locus_tag = segment.qualifiers['locus_tag'][0]
					seq_id = isolate_name+"|"+locus_tag+"|"+prod

					aa_chain = SeqRecord( Seq(segment.qualifiers['translation'][0],IUPAC.ExtendedIUPACProtein), 
								id=seq_id, name=seq_id )
					bla_genes[key]['isolate'].append( (isolate_name,locus_tag) )
					bla_genes[key]['aa'].append( aa_chain )

# Get list of keys
# geneTypes = bla_genes.keys()
# geneID = {k:n for n,k in enumerate(geneTypes)}
# isolate_cooccur = np.zeros( (len(geneTypes),len(geneTypes)) )
# molecule_cooccur = np.zeros( (len(geneTypes),len(geneTypes)) )

# nMolecules = 0
# nIsolates = 0
# for gbk in gbkFiles:
# 	isolate_name = "carb" + gbk.split('carb')[1].split('/')[0]

# 	isoSet = set([])
# 	# isoSet = []
# 	with open(gbk,'r') as geneF:
# 		for contig in SeqIO.parse(geneF,'genbank'):
# 			contigSet = set([])
# 			# contigSet = []
# 			for segment in contig.features:
# 				if (segment.type == 'CDS' and 'bla' in segment.qualifiers['product'][0]):
# 					prod = segment.qualifiers['product'][0]
# 					key = prod.split('-')[0].split('_')[0]
# 					# contigSet.append(geneID[key])
# 					contigSet.add(geneID[key])
# 			contigSet = list(contigSet)
# 			for n in xrange(len(contigSet)-1):
# 				for nn in xrange(n,len(contigSet)):
# 					molecule_cooccur[contigSet[n],contigSet[nn]] += 1
# 					if (contigSet[n] != contigSet[nn]):
# 						molecule_cooccur[contigSet[nn],contigSet[n]] += 1
# 			nMolecules += 1
# 			isoSet.update(set(contigSet))
# 			# isoSet.extend(contigSet)

# 	isoSet = list(isoSet)
# 	for n in xrange(len(isoSet)-1):
# 		for nn in xrange(n,len(isoSet)):
# 			isolate_cooccur[isoSet[n],isoSet[nn]] += 1
# 			if (isoSet[nn] != isoSet[n]):
# 				isolate_cooccur[isoSet[nn],isoSet[n]] += 1
# 	nIsolates += 1

# isolate_cooccur += np.diag(np.sum(isolate_cooccur,axis=1))
# molecule_cooccur += np.diag(np.sum(molecule_cooccur,axis=1))

# isolate_cooccur /= (1. * nIsolates)
# molecule_cooccur /= (1. * nMolecules)

for key,value in bla_genes.iteritems():
	if (not os.path.exists(key+".faa")):
		with open(key+".faa",'w+') as out_fasta:
			for aa in value['aa']:
				SeqIO.write(aa,out_fasta,'fasta')
