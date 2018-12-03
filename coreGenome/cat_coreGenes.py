import numpy as np
from scipy.special import loggamma

from collections import defaultdict,Counter 

import os 
import gzip 

from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq 
from Bio.Alphabet import generic_dna
from Bio import SeqIO 

# from treetime import TreeAnc 
import json

class mga(object):
    
    def __init__(self, panDir, overwrite=False, ref_strain=None, partition=False):

        panDir = panDir.rstrip("/")

        # Use predictability of folder structure output from panX
        isolateDir = panDir + "/input_GenBank/"
        clusterFile = panDir + "/vis/geneCluster.json"
        coreTree = panDir + "/vis/strain_tree.nwk"

        aln_dir = panDir + "/vis/geneCluster/"
        self.core_aln = panDir + "/core.fna"

        store_file = panDir + "/rec_anc.json"
        meta_file = panDir + "/meta_info.npz"

        if (not os.path.exists(self.core_aln) or overwrite):
            # Loading step. This simply organizes the data in a convienent manner.
            print("STARTING STEP 1")
            isolates = parseIsolates(isolateDir)
            isolates,clusters = parseJSON(clusterFile,isolates)

            self.isolates = list(isolates.keys())
            if (ref_strain is None or ref_strain not in self.isolates):
                self.ref_iso = self.isolates[0]
            else:
                self.ref_iso = ref_strain

            core_seqs = { iso:Seq('',alphabet=generic_dna) for iso in self.isolates }

            core_pos = []
            for locus,gene in isolates[self.ref_iso].items():
                c = gene['cluster']
                if ( c is not None and clusters[c]["dup"] == "no" and clusters[c]["count"] == len(self.isolates) ):
                    aln_file = aln_dir + clusters[c]["msa"] + "_na_aln.fa.gz"
                    if (os.path.exists(aln_file)):
                        with gzip.open( aln_file, "rt" ) as fna:
                            for segment in SeqIO.parse( fna, "fasta" ):
                                isoName = segment.name.split("-")[0]
                                core_seqs[isoName] += segment.seq
                        breaks = gene['start'] + gene['end']
                        core_pos.extend( np.linspace(min(breaks),max(breaks),len(segment.seq)).tolist() )
                    else:
                        print(aln_file)
            print(len(core_pos))
            np.savez(meta_file,np.array(core_pos),np.array([self.ref_iso]))

            with open(self.core_aln,"w+") as coreFasta:
                SeqIO.write( [SeqRecord(seq,iso,iso) for iso,seq in core_seqs.items()], coreFasta, 'fasta')

        print("FINISHED STEP 1")
        print("STARTING STEP 2")

        if (not os.path.exists(store_file) or overwrite):
            myTree = TreeAnc(gtr='Jukes-Cantor', tree=coreTree, aln=self.core_aln, verbose=0)
            myTree.infer_ancestral_sequences(infer_gtr=True, marginal=False)

            to_json = store_subtree(myTree.tree.clade)
            with open(store_file,"w+") as outF:
                json.dump(to_json,outF)

        print("FINISHED STEP 2")
        self.homoplasy = homoplasy(store_file)

        print("FINISHED STEP 3")
        if (partition):
            print("STARTING STEP 4")
            self.partitionDir = panDir + "/tree_partition/"
            if (not(os.path.exists(self.partitionDir))):

                os.mkdir(self.partitionDir)
                self.partition_alignment()
                print("Alt tree construction")
                call_fastTree(self.partitionDir+"alt.fna",self.partitionDir+"alt.nwk",submitJob=False)
                print("Ref tree construction")
                call_fastTree(self.partitionDir+"new.fna",self.partitionDir+"new.nwk",submitJob=False)

            print("FINISHED STEP 4")
        
            print("STARTING STEP 5")
            alt_store = self.partitionDir+"alt.json"
            if (not os.path.exists(alt_store)):
                myTree = TreeAnc(gtr='Jukes-Cantor', tree=self.partitionDir+"alt.nwk", aln=self.partitionDir+"alt.fna", verbose=0)
                myTree.infer_ancestral_sequences(infer_gtr=True, marginal=False)
                to_json = store_subtree(myTree.tree.clade)
                with open(alt_store,"w+") as outF:
                    json.dump(to_json,outF)
            print("FINISHED STEP 5")

            print("STARTING STEP 6")
            ref_store = self.partitionDir+"ref.json"
            if (not os.path.exists(ref_store)):
                myTree = TreeAnc(gtr='Jukes-Cantor', tree=self.partitionDir+"new.nwk", aln=self.partitionDir+"new.fna", verbose=0)
                myTree.infer_ancestral_sequences(infer_gtr=True, marginal=False)
                to_json = store_subtree(myTree.tree.clade)
                with open(alt_store,"w+") as outF:
                    json.dump(to_json,outF)
            print("FINISHED STEP 6")

    def partition_alignment(self,cutoff=.02):
        # import matplotlib.pylab as plt 
        from Bio import AlignIO 

        window = 10000
        loci,h,_ = self.homoplasy.scan(remove_singletons=True)
        nh = np.array( [ 1 for elem in h if len(elem) >= 1 ] )
        xs,nh_avg = self.homoplasy.window_average(loci,nh,W=window)

        sites = np.array( sorted(set(np.floor(xs[nh_avg > cutoff]))) )
        breakpoints = np.where( np.abs(np.diff(sites)) > window)[0]

        # breakpoints = np.asarray( [1] + breakpoints.tolist() + [len(sites)-1], dtype=int )
        sites = np.asarray(sites, dtype=int)
        
        start_ints = np.array( [sites[0]] + sites[breakpoints+1].tolist() )
        end_ints = np.array( sites[breakpoints].tolist() +  [sites[-1]] )

        np.savez(self.partitionDir + "alt_sites.npz",start_ints,end_ints)

        msa = AlignIO.read(self.core_aln,format="fasta")
        L = msa.get_alignment_length()

        for n,(s,e) in enumerate(zip(start_ints,end_ints)):
            if (n == 0):
                alt_msa = msa[:,s:e]
            else:
                alt_msa += msa[:,s:e]
        
        end_ints = np.array( [0] + end_ints.tolist() )
        start_ints = np.array( start_ints.tolist() + [L] )

        np.savez(self.partitionDir + "sites.npz",end_ints,start_ints)

        for n,(s,e) in enumerate(zip(end_ints,start_ints)):
            if (n == 0):
                new_msa = msa[:,s:e]
            else:
                new_msa += msa[:,s:e]

        with open(self.partitionDir+"alt.fna",'w+') as alt, open(self.partitionDir+"new.fna",'w+') as new:
            AlignIO.write(alt_msa,alt,"fasta")
            AlignIO.write(new_msa,new,"fasta")

        # alt_sites = []
        # for s,e in zip(start_ints,end_ints):
        #     alt_sites.extend( list( range(s,e) ) )
        # alt_sites = np.array(alt_sites)

        # m1 = np.zeros( max(loci) )
        # m1[alt_sites] = cutoff

        # plt.plot(xs,nh_avg)
        # plt.plot(m1)


        # plt.show()

def good_elem(elem):
    return elem[1] > 1 or elem[0][-1] == "P"

class homoplasy(object):

    def __init__(self, json_file):

        self.mut_pileup = defaultdict(lambda: [])
        self.tree_length = 0

        with open(json_file,'r') as inF:
            tree = json.load(inF)
            self.__store_subtree__(tree)

    def __store_subtree__(self, node):
        if ("children" in node.keys() and len(node['children'])>0):
            for child in node['children']:
                self.__store_subtree__(child)

        if ('mutation_length' not in node.keys()):
            self.tree_length += node["branch_length"] #node['attr']["clock_length"]
        else:
            self.tree_length += node['mut_length']
        
        # if ("children" in node.keys() and len(node['children'])>0):
        for mut in node['muts']:
            locus = int( mut[1:-1] )
            if ("children" in node.keys() and len(node['children'])>0):
                suffix = "P"
            else:
                suffix = "F"
            self.mut_pileup[locus].append( mut[0]+mut[-1]+suffix ) 

    def mutation_density(self,snp=True):
        
        if (snp):
            muts = [ (k, sum(1 for val in vs if check_mutation(val,0) ) ) for k,vs in self.mut_pileup.items() ]
        else:
            muts = [ (k, sum(1 for val in vs if check_mutation(val,1) ) ) for k,vs in self.mut_pileup.items() ]
        
        locus,mu = zip(*muts)
        locus = np.array(locus)
        mu = np.array(mu)

        indices = np.argsort(locus)
        locus = locus[indices]
        mu = mu[indices]

        locus = locus[mu > 0]
        mu = mu[mu > 0]

        return locus,mu

    def scan(self, remove_singletons=True):

        if (remove_singletons):
            muts = [ (k, list(Counter([ val for val in vs if check_mutation(val,0) ]).items()) ) for k,vs in self.mut_pileup.items() ]
        else:
            muts = [ (k, list(Counter([ val for val in vs if check_mutation(val,0) ]).items()) ) for k,vs in self.mut_pileup.items() ]
        loci,hplasy = zip(*muts)

        if (remove_singletons):
            loci = np.array([locus for n,locus in enumerate(loci) if is_homoplastic(hplasy[n]) ])
            hplasy = np.array([ [elem for elem in h if elem[1] > 1] for h in hplasy if is_homoplastic(h) ])
        else:
            loci = np.array([locus for n,locus in enumerate(loci) if len(hplasy[n]) ])
            hplasy = [ [ elem for elem in h if good_elem(elem) ] for h in hplasy if len(h) ]
            good_indices = [ n for n in range(len(hplasy)) if (len(hplasy[n])) > 0]
            hplasy = np.array(hplasy)
            loci = loci[good_indices]
            hplasy = hplasy[good_indices]
            print(hplasy.shape)
            print(loci.shape)
        indices = np.argsort(loci)
        loci = loci[indices]
        hplasy = hplasy[indices]

        hist = Counter()
        for h in hplasy:
            hist.update( Counter([ elem[1] for elem in h]) )

        return loci,hplasy,hist
    
    def window_average(self, x, y, xsample=None, W=5000):
        x = np.array(x)
        y = np.array(y)

        if (len(x) != len(y)):
            raise ValueError("Check dimensions of input arrays!")

        if (xsample is None):
            xsample = np.linspace(0,np.max(x),5000)
            
        y_avg = np.array( [ np.sum( y[ np.abs(x - xs) <= W ] ) / (2. * W + 1.) for xs in xsample ] )
        return xsample,y_avg

    def kl_divergence(self, W=5000, cutoff=None):
        from scipy.stats import poisson 
        from scipy.stats import entropy 

        mLoc,num_muts = self.mutation_density()
        num_muts[num_muts > 0] = 1 # Just count sites
        mLoc,mu = self.window_average(mLoc,num_muts)

        hloc,h,_ = self.scan(remove_singletons=False)

        empirical_dist = []
        for loc in mLoc:
            hInd = np.where( np.abs(hloc - loc) <= W )[0]
            C = Counter()
            for n in hInd:
                C.update( Counter( [elem[1] for elem in h[n]] ) )
            empirical_dist.append(count_to_dist(C, 2*W + 1) )

        KL = np.array( [ entropy(p,poisson.pmf(np.arange(len(p)),mu[n])) if p is not None and mu[n]>0 else 0 for n,p in enumerate(empirical_dist) ] )
        return KL, mu, mLoc

def parseJSON(homologs,isolates):
    """
    Parses the gene cluster JSON file output from PanX and assigns each gene
    the assigned cluster.
    """
    import json
    
    with open(homologs) as clusterFile:
        geneClusters = json.load(clusterFile)
        clusterData = []

        for n,cluster in enumerate(geneClusters):
            leaves = cluster['locus'].split(' ')
            meta_info = {"dup":cluster['dupli'],"count":cluster["count"],"msa":cluster["msa"],"loci":[]}

            for leaf in leaves:
                IDs = str(leaf).split('_')
                isolateName = IDs[0] 
                geneName = '_'.join( IDs[n] for n in range(1,len(IDs)) )
                meta_info["loci"].append((isolateName,geneName))

                if (isolateName in isolates):
                    isolates[isolateName][geneName]['cluster'] = n

            clusterData.append(meta_info)

    return isolates,clusterData

def parseIsolates(isolateDir):
    """
    Parses a directory containing gen bank files that contain isolate genomes
    used to build a pan-genome from the PanX pipeline.
    """
    from Bio import SeqIO
    import glob
    from collections import OrderedDict
    import scipy.stats
    
    fileList = glob.glob(isolateDir + "*.gbk")
    isolates = defaultdict()

    for gbk in fileList:
        isolateName = gbk.split('.gbk')[0].split(isolateDir)[1]
        genePartition = OrderedDict()
        numContig = 0

        for contig in SeqIO.parse(gbk,'genbank'):
            L = len(contig.seq)
            for feature in contig.features:

                if (feature.type == 'CDS'): # If it recognizes a gene, store its name and position on the reference.
                    geneName = feature.qualifiers['locus_tag'][0]
                    # Recover position of gene. Must deal with fuzzy boundaries and compound maps (multiple disjoint intervals)
                    pos = defaultdict(list)                    
                    for interval in feature.location.parts: 
                        pos['start'].append(interval.nofuzzy_start)
                        pos['end'].append(interval.nofuzzy_end)
                        pos['strand'].append(interval.strand)

                    strand = scipy.stats.mode(pos['strand'])
                    pos['strand'] = strand.mode[0]
                    pos['start'],pos['end'] = checkGenePosition(pos['start'],pos['end'],L)
                    pos['cluster'] = None
                    pos['contig'] = numContig
                    genePartition[geneName] = pos

            numContig += 1

        isolates[isolateName] = genePartition

    return isolates

def checkGenePosition(start,end,L):
    """
    Deals with the edge case that the gene straddles the contig start. If so,
    it swaps the order so that the first start point is to the left of the 
    successive start points
    """
  
    if (len(start) > 1):
        delta = np.diff(start)
        bigDiff = np.where(delta > L/2.)[0]
        if (len(bigDiff) > 0):
            rollAmt = -(bigDiff[0]+1)
            return np.roll(start,rollAmt).tolist(),np.roll(end,rollAmt).tolist()
        else:
            return start,end
    else:
        return start,end

def returnGeneList(n,isolates):
    names = isolates.keys()

    cogs = [isolates[names[n]][geneName]['cluster'] for geneName in isolates[names[n]].keys() if isolates[names[n]][geneName]['cluster'] is not None]
    dirs = [isolates[names[n]][geneName]['strand'] for geneName in isolates[names[n]].keys() if isolates[names[n]][geneName]['cluster'] is not None]
    gene = [geneName for geneName in isolates[names[n]].keys() if isolates[names[n]][geneName]['cluster'] is not None]
    start = [isolates[names[n]][geneName]['start'] for geneName in isolates[names[n]].keys() if isolates[names[n]][geneName]['cluster'] is not None]
    end = [isolates[names[n]][geneName]['end'] for geneName in isolates[names[n]].keys() if isolates[names[n]][geneName]['cluster'] is not None]

    return cogs,dirs,gene,start,end

def store_subtree(node):
    node_dict = {}

    node_dict["strain"] = node.name
    node_dict["branch_length"] = node.branch_length
    node_dict["branch_time"] = node.mutation_length
    node_dict["muts"] = [ m[0] + str(m[1]) + m[2] for m in node.mutations ]

    if ( len(node.clades) > 0 ):
        node_dict["children"] = []
        for child in node.clades:
            node_dict["children"].append( store_subtree(child) )

    return node_dict


def check_mutation(transition, case):
    if (case == 0):
        return "N" not in transition and "-" not in transition
    elif (case == 1):
        return "N" not in transition
    elif (case == 2):
        return "N" not in transition and "-" not in transition and "F" not in transition

def is_homoplastic(muts):
    if (len(muts) == 0):
        return False
    else:
        homoplastic = False
        for m in muts:
            if (m[1] > 1):
                homoplastic = True
                break
        return homoplastic
    
def count_to_dist(C,W):
    if (len(C)):
        numEvents = sum(C.values())
        numNonEvents = W-numEvents
        C[0] = numNonEvents
        kMax = max(C.keys())
        P = np.array( [ C[k]/(1.*W) if k in C.keys() else 0 for k in range(kMax+1) ] )
        return P
    else:
        return None

def entropy(P):
    return np.sum( P*np.array( [ np.log(p) if p > 0 else 0 for p in P] ) )

def mean_dist(P):
    return np.sum( P * np.arange(len(P)) )

def log_factorial_dist(P):
    return np.sum(P * np.real(loggamma(np.arange(len(P)) + 1)) )

def call_fastTree(inputFile,outputFile,submitJob=False):
    import subprocess
    from socket import gethostname

    workSpace = gethostname()
    fastTreeExist = True
    
    if (workSpace == 'nicholasWorkStation'):
        fastTreeBin = '~/Documents/postdoc/'
    elif (workSpace.endswith('.cluster.bc2.ch')):
        fastTreeBin = '/scicore/home/neher/GROUP/bin/'
    else:
        fastTreeExist = False # Tree won't assemble
        raise ValueError("Fast Tree not found") 

    if (fastTreeExist):  
        # Reads aligned sequences from a fasta file.
        
        if (workSpace.endswith('.cluster.bc2.ch')):
            # fastTreeCMD = "sbatch --export=input='%s',output='%s',tmp_output='%s' /scicore/home/neher/nolln/genomeAssembly/submitScripts/fastTree.sh"%(inputFile,outputFile)
            if (submitJob):
                fastTreeCMD = "sbatch --export=input='%s',output='%s' /scicore/home/neher/nolln/genomeAssembly/submitScripts/fastTree.sh"%(inputFile,outputFile)
            else:
                fastTreeCMD = fastTreeBin +"FastTree -quiet -nt " + inputFile + " > " + outputFile
        else:
            fastTreeCMD = fastTreeBin + "FastTree -quiet -nt " + inputFile + " > " + outputFile

        # print(fastTreeCMD)
        process = subprocess.Popen(fastTreeCMD,shell=True)
        out, err = process.communicate()