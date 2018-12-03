from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import defaultdict
import glob
import numpy as np
from itertools import chain,combinations
import os
import json

def all_subsets(ss):
	return chain(*map(lambda x: combinations(ss,x), range(0,len(ss)+1)))

topDir = '/scicore/home/neher/GROUP/data/Carbapenemases/'
gbkWC = 'carb???/final/pilon/prokka/*gbk'

gbkFiles = glob.glob(topDir+gbkWC)
with open("geneClusters.json","r") as gc:
	defined_clusters = json.load(gc)

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

						if (key in defined_clusters):
							family_name = defined_clusters[key][isolate_name][locus_tag.replace("_"," ")]
						else:
							family_name = key

						aa_chain = SeqRecord( Seq(segment.qualifiers['translation'][0],IUPAC.ExtendedIUPACProtein), 
									id=seq_id, name=seq_id )
						nuc_seq  = SeqRecord( segment.extract(contig.seq),id=seq_id, name=seq_id )

						family_name = family_name.replace(" ","_")
						bla_genes[family_name]['isolate'].append( (isolate_name,locus_tag) )
						bla_genes[family_name]['aa'].append( aa_chain )
						bla_genes[family_name]['nc'].append( nuc_seq )

	for key,value in bla_genes.iteritems():
		if (not os.path.exists(key+".faa")):
			with open("fasta/"+key+".faa",'w+') as out_faa, open("fasta/"+key+".fna",'w+') as out_fna:
				for aa in value['aa']:
					SeqIO.write(aa,out_faa,'fasta')
				for nuc in value['nc']:
					SeqIO.write(nuc,out_fna,'fasta')
