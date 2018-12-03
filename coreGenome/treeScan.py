from Bio import AlignIO
from Bio import SeqIO
import numpy as np

import os 
from collections import defaultdict
import networkx as nx 
import ete3 

import tables
import tempfile 
import glob


class keyedList(object):
    def __init__(self,l,key):
        self.l = l
        self.key = key
    def __len__(self):
        return len(self.l)
    def __getitem__(self,index):
        return self.key(self.l[index])

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


class scanner(object):

    def __init__(self, topDir):
        
        # Prepare h5 file
        self.dir = topDir + "/"

        self.h5name = topDir + "/src.h5"
        self.aln = topDir + "/core.fna"
        self.scanName = "scan"
        self.blockSize = 3e3
        self.tempTree = topDir + "/tmp.nwk"
        self.strainTree = topDir + "/vis/strain_tree.nwk"
        self.branch_cutoff = 5e-4

        self.isolate_names = [ x.rstrip(".gbk").split('/')[-1] for x in glob.glob(topDir+"/input_GenBank/*") ]
        self.numStrains = len(self.isolate_names)

        if (not os.path.exists(self.h5name)):
            db = tables.open_file(self.h5name,mode="w",title="home")
            db.create_group("/",self.scanName,"Reconstructed trees from 10 kB blocks")
            db.close()
 
    def scan(self,overwrite=False):
        db = tables.open_file(self.h5name,mode="a")
        treeScan = db.get_node("/",self.scanName)

        if (not db.__contains__("/" + self.scanName + "/trees") or overwrite):
            full_msa = AlignIO.read(self.aln,"fasta")
            L = full_msa.get_alignment_length()

            numBlocks = np.floor(L/self.blockSize)
            breakpoints = np.linspace(0,L,num=numBlocks,dtype=np.int)
            intervals = [( breakpoints[n], breakpoints[n+1]) for n in range(len(breakpoints)-1)]
            if (L - breakpoints[-1] > 1000):
                intervals += [ (breakpoints[-1],L) ]

            trees = []
            true_intvals = []

            for interval in intervals:
                sub_msa = full_msa[:,interval[0]:interval[1]]
                frac_gaps = sum( record.seq.count("-")/(interval[1]-interval[0]) for record in sub_msa ) / len(sub_msa)

                if (frac_gaps < 1):  
                    tmp_alnF = tempfile.NamedTemporaryFile(mode="w")
                    AlignIO.write(sub_msa,tmp_alnF,'fasta')
                    tmp_alnF.flush()
                    call_fastTree(tmp_alnF.name,self.tempTree,submitJob=False)

                    trees.append( self.__newick__() )
                    true_intvals.append(interval)

            tree = np.array(trees)
            breakpoints = np.array(intervals)

            if (overwrite):
                del db["/" + self.scanName + "trees"]
                del db["/" + self.scanName + "ivals"]

            db.create_array(treeScan,"trees",tree)
            db.create_array(treeScan,"ivals",true_intvals)

        db.flush()
        db.close()

        if (os.path.exists(self.tempTree)):
            os.remove(self.tempTree)

    def grab_trees(self,db=None,collapse=True):

        if (db is None):
            db = tables.open_file(self.h5name,mode="r")
            close = True
        else:
            close = False
        trees = [ete3.Tree(tree.decode('utf-8')) for tree in db.get_node("/" + self.scanName + "/trees",classname="Array")]
        intvals = [(x[0],x[1]) for x in db.get_node("/" + self.scanName + "/ivals",classname="Array")]

        if (close):
            db.close()

        if (collapse):
            for tree in trees:
                for node in tree.get_descendants():
                    if (not node.is_leaf() and node.dist < self.branch_cutoff):
                        node.delete()

        return trees, intvals

    def compute_dist_matrix(self, dendropy=False, weighted=False, resolve=True, overwrite=False):

        import dendropy 
        from dendropy.calculate import treecompare

        db = tables.open_file(self.h5name,mode="a")
        trees,intvals = self.grab_trees(db)

        if (not db.__contains__("/" + "dist_matrix") or overwrite):
            D = np.zeros( (len(trees), len(trees)) )
            if (not dendropy):
                for n in range(len(trees)-1):
                    for nn in range(n+1,len(trees)):
                        D[n,nn] = self.compare_trees(trees[n],trees[nn])
                        D[nn,n] = D[n,nn]
            else:
                T = dendropy.TreeList([dendropy.Tree.get(data=t.write(),schema='newick') for t in trees])
                for n in range(len(trees)-1):
                    for nn in range(n+1,len(trees)):
                        if (weighted):
                            w_rf = treecompare.euclidean_distance(T[n],T[nn]) #weighted_robinson_foulds_distance(T[n],T[nn])
                        else:
                            w_rf = treecompare.symmetric_difference(T[n],T[nn]) #weighted_robinson_foulds_distance(T[n],T[nn])
                        D[n,nn] = w_rf
                        D[nn,n] = w_rf
            if (overwrite):
                del db["/dist_matrix"]
            db.create_array("/","dist_matrix",D)
        else:
            D = np.array([np.array(row) for row in db.get_node("/dist_matrix",classname="Array")])

        db.flush()
        db.close()
        return D

    def compare_to_strain_tree(self,rf=True):
        import dendropy 
        from dendropy.calculate import treecompare
        from ete3 import Tree 

        if (rf):
            with open(self.strainTree,'r') as ref:
                Ref = dendropy.Tree.get(file=ref,schema="newick")

            T = dendropy.TreeList( [Ref] + [dendropy.Tree.get(data=t.write(),schema='newick') for t in self.grab_trees()[0]] )
            d = np.zeros(len(T)-1)
            for n in range(1,len(T)):
                d[n-1] = treecompare.symmetric_difference(T[0],T[n])
        else:
            Ref = Tree(self.strainTree,format=1)
            trees = self.grab_trees()[0]
            d = np.zeros(len(trees))
            for n,tree in enumerate(trees):
                d[n] = self.compare_trees(Ref,tree)
                
        return d

    def compare_trees(self,A,B):
        childrenA = A.get_cached_content()
        childrenB = B.get_cached_content()

        leafA = set([node.name for node in childrenA[A]])
        leafB = set([node.name for node in childrenB[B]])

        for node in A.traverse():
            if (not node.is_root()):
                tmpA = set( [node.name for node in childrenA[node]] )
                childrenA[node] = (tmpA,leafA.difference(tmpA))
        
        for node in B.traverse():
            if (not node.is_root()):
                tmpB = set( [node.name for node in childrenB[node]] )
                childrenB[node] = (tmpB,leafB.difference(tmpB))

        D = 0
        for node_a in A.traverse():
            if (not node_a.is_root()):
                for node_b in B.traverse():
                    if (not node_b.is_root()):
                        D += ( (not childrenA[node_a][0].isdisjoint(childrenB[node_b][0])) and \
                               (not childrenA[node_a][1].isdisjoint(childrenB[node_b][1])) and \
                               (not childrenA[node_a][0].isdisjoint(childrenB[node_b][1])) and \
                               (not childrenA[node_a][1].isdisjoint(childrenB[node_b][0])) )
        return D
    
    def cluster_trees(self,numTrees=8,cluster_cutoff=5,dist_cutoff=5,method='centroid'):
        import scipy.cluster.hierarchy as hcl
        import scipy.spatial.distance as ssd 
        from collections import Counter

        D = self.compute_dist_matrix()
        trees,intvals = self.grab_trees()

        S = ssd.squareform(D)
        Z = hcl.linkage(S,method=method)

        clusters = hcl.fcluster(Z,cluster_cutoff,criterion='distance')
        multiplicity = Counter(clusters)

        isolate_types = Counter([name.split('_')[0] if "nan" not in name else name for name in self.isolate_names])
        tree_clusters = [ common_clus[0] for common_clus in multiplicity.most_common(numTrees) ]

        clusters = np.array(clusters)
        globedClusters = np.zeros_like(clusters) 
        minDist = np.zeros_like(clusters)
        DMax = np.max(D)

        for n,subCluster in enumerate(multiplicity.most_common()):
            tree_indx = np.where(clusters==subCluster[0])[0]
            if (n < numTrees):
                globedClusters[tree_indx] = n
            else:
                clust_dist = np.array([ np.median(np.take(np.take(D,tree_indx,axis=0),np.where(clusters==base_cluster)[0],axis=1).flatten()) for base_cluster in tree_clusters ])
                min_ind = np.argmin(clust_dist)
                minDist[tree_indx] = clust_dist[min_ind]

                if (clust_dist[min_ind] < dist_cutoff):
                    globedClusters[tree_indx] = min_ind
                else:
                    globedClusters[tree_indx] = -1
        
        tree_intvals = {}
        for n in range(numTrees):
            tree_intvals[n] = [x for m,x in enumerate(intvals) if globedClusters[m] == n ]
            
        return tree_intvals,minDist

    def partition_coreFA(self,numTrees=8,cluster_cutoff=13,dist_cutoff=5,method="centroid"):
        import os, shutil

        from Bio.Alphabet import generic_dna
        from Bio import AlignIO
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        partDir = self.dir + "/partitions"
        if (os.path.exists(partDir)):
            shutil.rmtree(partDir)
        os.mkdir(partDir)

        coreFNA = AlignIO.read(open(self.aln,'r'),'fasta') 
        names = [ rec.id for rec in coreFNA ]
        coreFNA = np.array( [ list(rec) for rec in coreFNA], np.character ) 

        tree_intvals, _ = self.cluster_trees(numTrees,cluster_cutoff,dist_cutoff,method)
        for n,val in tree_intvals.items():
            sub_msa = np.concatenate( [ coreFNA[:,v[0]:v[1]] for v in val ], axis=1)
            # print(sub_msa)
            write_msa = MultipleSeqAlignment( [SeqRecord( Seq("".join(c.decode('utf-8') for c in sub_msa[m,:]), generic_dna), id=names[m]) for m in range(sub_msa.shape[0]) ] )
            AlignIO.write(write_msa, partDir + "/core_partition%02d.fna"%n,"fasta")

    def build_partition_trees(self):
        import os, glob

        partDir = self.dir + "/partitions"
        if (not os.path.exists(partDir)):
            raise ValueError("First partition the core genome!")
        else:
           files = glob.glob(partDir + "/*.fna")
           print(files)
           for f in files:
               outName = f.rstrip(".fna") + ".nwk"
               call_fastTree(f,outName,submitJob=True)

    def plot_cluster_trees(self, cutoff=1, repNum=10, W=500, H=500, M=20, f=.25, TreeM=50, numTrees=4, numRows=1, method='centroid'):
        import scipy.cluster.hierarchy as hcl
        import scipy.spatial.distance as ssd 
        import pickle
        from collections import Counter
        import toyplot
        from ete3 import Tree 
        # import toyplot.svg 

        D = self.compute_dist_matrix()
        trees,intvals = self.grab_trees()

        S = ssd.squareform(D)
        Z = hcl.linkage(S,method=method)

        clusters = hcl.fcluster(Z,cutoff,criterion='distance')
        multiplicity = Counter(clusters)

        # isolate_colors = toyplot.color.brewer.map("Spectral")
        isolate_colors = pickle.load( open(self.dir + "isoColors.pkl","rb") )
        isolate_types = Counter([name.split('_')[0] if "nan" not in name else name for name in self.isolate_names])
        # colors = toyplot.color.broadcast(isolate_colors,shape=len(isolate_types))

        # isolate_colors = {}
        # isolate_shapes = {"o":[],"^":[]}

        # nColor = 0
        # for iso_type,num_types in isolate_types.items():
        #     for n in range(num_types):
        #         if ('nan' in iso_type):
        #             name = iso_type
        #         else:
        #             name = iso_type + "_i" + str(n)
        #         isolate_colors[ name ] = colors[nColor]
                #
                # if (n == 0):
                #     isolate_shapes["o"].append(name)
                # else:
                #     isolate_shapes["^"].append(name)

            # nColor += 1

        # isolate_shapes["o"] = set(isolate_shapes["o"])
        # isolate_shapes["^"] = set(isolate_shapes["^"])

        partitions = np.matlib.repmat(clusters,repNum,1)

        base_palette = toyplot.color.brewer.palette("Set1")
        tree_clusters = [ common_clus[0] for common_clus in multiplicity.most_common(numTrees) ]

        clusters = np.array(clusters)
        cluster_palette = [None] * (np.max(clusters)+1)
        
        DMax = np.max(D)
        for n,subCluster in enumerate(multiplicity.most_common()):
            if (n < numTrees):
                cluster_palette[subCluster[0]] = base_palette.color(n)
            else:

                clust_dist = np.array([ np.median(np.take(np.take(D,np.where(clusters==subCluster[0])[0],axis=0),np.where(clusters==base_cluster)[0],axis=1).flatten()) \
                                      for base_cluster in tree_clusters ])
                clust_dist = DMax - clust_dist
                clust_dist /= np.sum(np.power( clust_dist, 1))
                cL = sum( c * toyplot.color.to_lab(base_palette.color(m))[0] for m,c in enumerate(clust_dist) )
                cA = sum( c * toyplot.color.to_lab(base_palette.color(m))[1] for m,c in enumerate(clust_dist) )
                cB = sum( c * toyplot.color.to_lab(base_palette.color(m))[2] for m,c in enumerate(clust_dist) )
                cluster_palette[subCluster[0]] = toyplot.color.lab(cL,cA,cB)

        for n,color in enumerate(cluster_palette):
            if color is None:
                cluster_palette[n] = toyplot.color.black

        palette = toyplot.color.Palette(colors=np.array(cluster_palette))
                
        colormap = toyplot.color.CategoricalMap(palette)

        ## Plot with toyplot into canvas 
        canvas = toyplot.Canvas(width=W,height=H)
        table = canvas.matrix((partitions,colormap),bounds=(M,(W-M),M,f*(H-M)),margin=0,bshow=False,tshow=False,rshow=False,lshow=False)
        table.body.gaps.columns[...] = 0
        table.body.gaps.rows[...] = 0
        
        trees_per_row = int(numTrees / numRows)

        yMin = f*(H-M)
        deltaY = ( (1-f)*(H-M) ) / numRows
        yMax = yMin + deltaY - 1.2*TreeM/2

        deltaX = (W-M) / (1. * trees_per_row)

        xMin = M
        xMax = xMin + deltaX - TreeM

        for n in range(numTrees):
            if (n % trees_per_row == 0 and n > 0):
                yMin = yMax + 1.2*TreeM/2
                yMax += deltaY

                xMin = M
                xMax = xMin + deltaX - TreeM

            # xMin = M + n*(W-M) / (1. * numTrees)
            # xMax = M + (n+1)*(W-M) / (1. * numTrees) - TreeM
            tree_axes = canvas.cartesian( bounds=(xMin,xMax,yMin,yMax) )

            xMin = xMax + TreeM
            xMax += deltaX

            # tree_indx = next(i for i in range(len(clusters)) if clusters[i]==tree_clusters[n])
            trees_in_cluster = [(i,len(tree.get_edges())) for i,tree in enumerate(trees) if clusters[i] == tree_clusters[n]]
            tree_indx = max(list(range(len(trees_in_cluster))),key = lambda x: trees_in_cluster[x][1] )
            tree_indx = trees_in_cluster[tree_indx][0]
            # if (n == 0):
                # isolate_colors = self.plot_unrooted_tree(trees[tree_indx],tree_axes,None,isolate_shapes)
            # else:
            curT = Tree( self.dir + "partitions/core_partition%02d.nwk"%n )
            self.plot_unrooted_tree(curT,tree_axes,isolate_colors) #,isolate_shapes)

            tree_axes.label.text = "Partition " + str(n+1)
            tree_axes.label.style = {"fill":palette.css(tree_clusters[n]),"font-weight":"bold","text-anchor":"middle"} #colors[tree_clusters[n]]

            tree_axes.padding = 20
            tree_axes.show = True
            tree_axes.x.ticks.show = False
            tree_axes.y.ticks.show = False
            if (n % trees_per_row == 0):
                tree_axes.y.show = False
            tree_axes.x.show = False
            tree_axes.x.ticks.labels.show = False
            tree_axes.y.ticks.labels.show = False

        return (Z,clusters,partitions,len(multiplicity)),canvas

    def plot_unrooted_tree(self, tree, axes, tip_colors=None, tip_shapes=None):
        import toyplot

        # Count the number of children under each node
        for node in tree.traverse('postorder'):
            if (node.is_leaf()):
                node.add_feature('leafCount',1)
            else:
                numLeafs = 0
                for child in node.children:
                    numLeafs += child.leafCount
                node.add_feature('leafCount',numLeafs)

        # Initialize the root at the origin.       
        tree.add_features(x=0,y=0,tau=0,w=2*np.pi)
        self.__place_unrooted_subtree__(tree)

        # Plot as graph in toyplot
        for n,node in enumerate(tree.traverse('preorder')):
            if (node.name == ""):
                node.name = "INT" + str(n)

        name_id = {node.name:n for n,node in enumerate(tree.traverse('preorder'))}
        edges = np.array([ [name_id[node.name],name_id[node.up.name] ] for node in tree.traverse('preorder') if not node.is_root() ])

        # Get calculate layout 
        x = [ node.x for node in tree.traverse('preorder') ]
        y = [ node.y for node in tree.traverse('preorder') ]

        vcoords = np.array([x,y]).T
        axes.graph(edges,
                   vcoordinates=vcoords,
                   vlshow=False,
                   vcolor=None,
                   vopacity=0,
                   vsize=0,
                   ecolor='black',
                   eopacity=.8
                  )
        axes.show = False

        # Scatter plot for tip markers.
        ctX = [] 
        ctY = []
        ctColor = []
        cName = []

        if (tip_colors is None):
            angular_tip_colors = defaultdict(list)

        for node in tree.traverse('preorder'):
            if (node.name in self.isolate_names):
                ctX.append(node.x)
                ctY.append(node.y)
                if (tip_colors is not None):
                    ctColor.append(tip_colors[node.name])
                else:
                    angular_tip_colors[node.name].append(np.arctan2(node.y,node.x))
                    # ctColor.append("black")
                cName.append(node.name)

        if (tip_colors is None):
            avg_angle = [ (key,np.mean(value)) for key,value in angular_tip_colors.items() ]
            keys,values = zip(*avg_angle)
            index = np.argsort(np.array(values))

            isolate_colors = toyplot.color.brewer.map("Spectral")
            colors = toyplot.color.broadcast(isolate_colors,shape=len(index))
            angular_tip_colors = {keys[ii]:colors[ii] for ii in index}

            for node in tree.traverse('preorder'):
                if (node.name in self.isolate_names):
                    ctColor.append(angular_tip_colors[node.name])


        ctX = np.array(ctX)
        ctY = np.array(ctY)

        darkColors = [ toyplot.color.rgba(.65*x.tolist()[0],.65*x.tolist()[1],.65*x.tolist()[2],.65*x.tolist()[3]) for x in ctColor ]

        axes.scatterplot(ctX,ctY,color=darkColors,size=8,title=cName)
        axes.scatterplot(ctX,ctY,color=ctColor,size=6,title=cName)

        if (tip_colors is None):
            return angular_tip_colors
        else:
            return None

    def __place_unrooted_subtree__(self,node):
        if (not node.is_root()):
            px = node.up.x
            py = node.up.y 
            angle = node.tau + node.w*.5 
            Dl = 1e-4*np.log10( (1.0*node.dist)/1e-4 + 1)
            node.add_features(x= px + Dl*np.cos(angle),y= py + Dl*np.sin(angle) )

        eta = node.tau
        for child in node.children:
            child.add_features(w = 2*np.pi*child.leafCount / self.numStrains, tau=eta )
            eta += child.w 
            self.__place_unrooted_subtree__(child)

    def __newick__(self):
        nwk = ""
        with open(self.tempTree,'r') as tmpTree:
            nwk = tmpTree.read().replace('\n','')
        return nwk


        
def discrete_cmap(N, base_cmap=None):
    import matplotlib.pylab as plt

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0,1,N))
    cmap_name = base.name + str(N)

    return base.from_list(cmap_name,color_list,N)


def plotAccessoryDensity(W=500,H=100,M=5,strain="ecoli",ref="carb013"):
    import synteny.analysis as analysis
    import toyplot 

    exp = analysis.strains(strain + "/")
    V = exp.synGraph.coreEnvVar(ref.encode('utf-8'))
    V = np.matlib.repmat(V,2,1)

    canvas = toyplot.Canvas(width=W,height=H)
    dat = V
    colormap = toyplot.color.linear.map("Blackbody",domain_min=np.min(dat), domain_max=np.max(dat))
    table = canvas.matrix( (dat,colormap) ,bounds=(M,(W-M),M,(H-M)),margin=0,bshow=False,tshow=False,rshow=False,lshow=False)
    return canvas,table



