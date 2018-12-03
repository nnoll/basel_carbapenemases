from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from shutil import copyfile
import sys,os
import json
from collections import defaultdict

topFolder = "/scicore/home/neher/GROUP/data/Carbapenemases/"
gbkLoc = "/final/pilon/prokka/entero.gbk"
clusterFolder = "/scicore/home/neher/GROUP/data/Carbapenemases/blaTrees/molecules/"

with open("geneClusters.json","r") as gc:
	defined_clusters = json.load(gc)
	gene_families = defined_clusters.keys()

	for gf in gene_families:
		# print gf
		subFamilies = [k for k in defined_clusters[gf].keys() if 'carb' not in k]
		for sub_gf in subFamilies:
			locus_tags = defaultdict(list)
			# print len(defined_clusters[gf][sub_gf]['members'])
			for member in defined_clusters[gf][sub_gf]['members']:
				locus_tags[member[0]].append(member[1])

			genomePos = defaultdict(list)
			for isolate,locus in locus_tags.iteritems():
				locusNames = set([l.replace(" ","_") for l in locus])
				# print locusNames
				with open(topFolder + isolate + gbkLoc,'r') as gbk:
					for n,contig in enumerate(SeqIO.parse(gbk,'genbank')):
						for m,segment in enumerate(contig.features):
								if (segment.type == 'CDS' and segment.qualifiers['locus_tag'][0] in locusNames):
									genomePos[isolate].append((n,m))

			if (len(genomePos) > 2):
				sub_gf_folder = clusterFolder + sub_gf.replace(" ","_") +"/"
				ncFolder = sub_gf_folder + "nc/"
				aaFolder = sub_gf_folder + "aa/"

				# sub_gf_folder = clusterFolder + sub_gf.replace(" ","_") +"/"
				os.mkdir(sub_gf_folder)
				os.mkdir(ncFolder)
				os.mkdir(aaFolder)

				for isolate,pos in genomePos.iteritems():
					for iso,x in enumerate(pos):
						contigNum = x[0]
						geneNum = x[1]
						with open(topFolder + isolate + gbkLoc,'r') as gbk, open(sub_gf_folder+isolate+"_%02d.gbk"%iso,"w+") as new_gbk, \
							open(ncFolder+"nucleotide.fna","a+") as fna, open(aaFolder+"protein.faa","a+") as faa:
							for n,contig in enumerate(SeqIO.parse(gbk,'genbank')):
								if (n == contigNum and n == 0):
									subContig = SeqRecord(contig.seq,contig.id,contig.name)
									for m, segment in enumerate(contig.features):
										if (abs(m-geneNum) <= 50):
											subContig.features.append(segment)
											if (m == geneNum): # Save gene sequence
												seq_id = isolate + "_%02d|"%iso + "|" + segment.qualifiers['product'][0]
												aa_chain = SeqRecord( Seq(segment.qualifiers['translation'][0],
												IUPAC.ExtendedIUPACProtein),id=seq_id, name=seq_id )
												nuc_seq  = SeqRecord( segment.extract(contig.seq),id=seq_id, name=seq_id )
												SeqIO.write(aa_chain,faa,'fasta')
												SeqIO.write(nuc_seq,fna,'fasta')
									break
								elif (n == contigNum and n > 0):
									subContig = contig
									for m, segment in enumerate(contig.features):
										if (m == geneNum): # Save gene sequence
											# print segment.qualifiers['product'][0]
											seq_id = isolate + "_%02d|"%iso + "|" + segment.qualifiers['product'][0]
											aa_chain = SeqRecord( Seq(segment.qualifiers['translation'][0],
											IUPAC.ExtendedIUPACProtein),id=seq_id, name=seq_id )
											nuc_seq  = SeqRecord( segment.extract(contig.seq),id=seq_id, name=seq_id )
											SeqIO.write(aa_chain,faa,'fasta')
											SeqIO.write(nuc_seq,fna,'fasta')
									break
							SeqIO.write(subContig,new_gbk,'genbank')
# gene = sys.argv[1]
# if (os.path.exists("./"+gene)):
#     input_dir = '../input_GenBank/'
#     output_dir = "./" + str(gene) + "/"
#     fasta = False

#     if (fasta):
#         with open(output_dir+gene+'_plasmids.txt','r') as kpcP, open(output_dir+gene+"_plasmids.faa",'w+') as faa, open(output_dir+gene+"_plasmids.faa",'w+') as fna:
#             fileNames = set([line.split(".gbk")[0] for line in kpcP.readlines()])
#             for name in fileNames:
#                 gbkName = input_dir + name + '.gbk'
#                 print gbkName
#                 with open(gbkName,'r') as gbk:
#                     for plasmid in SeqIO.parse(gbk,"genbank"):
#                         SeqIO.write(plasmid,fna,"fasta")
#                         plasmidProtein = SeqRecord(plasmid.seq.translate(), plasmid.id,plasmid.name)
#                         SeqIO.write(plasmidProtein,faa,"fasta")
#     else:
#         with open(output_dir+gene+'_plasmids.txt','r') as kpcP:
#             fileNames = set([line.split(".gbk")[0] for line in kpcP.readlines()])
#             for name in fileNames:
#                 name = name.replace("input_GenBank/","")
#                 gbkName = input_dir + name + '.gbk'
#                 new_gbkName = output_dir + name + '.gbk'
#                 copyfile(gbkName,new_gbkName)
#                 # print gbkName
# else:
#     raise ValueError("Not a recognized gene")