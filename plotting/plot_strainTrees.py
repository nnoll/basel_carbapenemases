import json 
import os 
from collections import defaultdict
from genomePlotter import plotter
import glob 

import numpy as np
import matplotlib.pylab as plt 

import toytree 
import toyplot 

inFile = "mnt/neherGroup/Carbapenemases/blaTrees/geneClusters.json"
outFile = "mnt/neherGroup/Carbapenemases/panGenome/pa.json"

def gather_PA(overwrite=False):
    if (overwrite or not os.path.exists(outFile)):
        pa_db = defaultdict(lambda: list())
        pa_db['categories'] = set([])
        with open(inFile,'r') as gc:
            db = json.load(gc)
            for k in db.keys():
                for kk in db[k].keys():
                    if ('carb' in kk):
                        genes = list(db[k][kk].values())
                        pa_db[kk] += genes
                        pa_db['categories'].update(genes)

        pa_db['categories'] = list(pa_db['categories'])
        with open(outFile,'w+') as pa:
            json.dump(pa_db,pa)

def plot_tree(species,headSize=5,frac=.3,W=700,H=500,ydiff=.1,xdiff=.045):
    from ete3 import Tree
    import pandas as pd
    import pickle
    import matplotlib

    if (not os.path.exists(outFile)):
        gather_PA()
    
    with open(outFile,'r') as in_db:
        db = json.load(in_db)

    folder = "mnt/neherGroup/Carbapenemases/panGenome/proc/" + species.rstrip("/") + "/"
    treeFile = folder + "vis/strain_tree_ml.nwk"

    mlst = {}
    tbl = pd.read_csv("mnt/neherGroup/Carbapenemases/carbapenamase_seq_runs/ST_types.tsv",sep = "\t")
    for index, row in tbl.iterrows():
        isolate = row['isolate']
        mlst[isolate] = row['ST']

    mlsts_id = { st:n for n,st in enumerate(set(mlst.values()))}
    mlst_db = { iso: mlsts_id[st] for iso,st in mlst.items() }

    tmp_tree = Tree(treeFile,format=1)
    tmp_tree.set_outgroup( tmp_tree.get_midpoint_outgroup() )

    tree = toytree.tree(tmp_tree.write(format=0),format=0)

    isolates = tree.get_tip_labels()
    reduced_cats = {"CTX-M-24" : "CTX-M", "CTX-M-15" : "CTX-M", "KPC-3":"KPC", "NDM-4":"NDM", "NDM-6":"NDM", "SHV-11":"SHV-11", "SHV-143":"SHV-11", "OXA-320":"OXA-48", "OXA-48":"OXA-48" }
    all_cats = list( set( [ reduced_cats[k.replace("bla","").split(" ")[0]] if k.replace("bla","").split(" ")[0] in reduced_cats else k.replace("bla","").split(" ")[0] for k in db['categories'] ] ) )
    
    pa_matrix = np.zeros( (len(isolates),len(all_cats)) )
    cat_id = { cat:n for n,cat in enumerate(all_cats)}

    tip_labels = [ "%03d"%int(mlst[isolate].replace("ST","")) if mlst[isolate] is not "-" else "---" for isolate in isolates ]
    for n,isolate in enumerate(isolates):
        for g in db[isolate]:
            G = g.replace("bla","").split(" ")[0]
            if (G in reduced_cats):
                G = reduced_cats[G]
            pa_matrix[n,cat_id[G]] = 1
        # pa_matrix[n,0] = mlst[isolate]

    labels = np.array(all_cats) #[""]*pa_matrix.shape[1]

    index = ( np.sum(pa_matrix,axis=0) > 1 ) & ( np.sum(pa_matrix,axis=0) < pa_matrix.shape[0] )
    # index = np.array(index.tolist() )
    pa_matrix = pa_matrix[:,index]
    labels = labels[index]
    sort_index = np.argsort(np.sum(pa_matrix,axis=0))
    pa_matrix = pa_matrix[:,sort_index[::-1]]
    labels = labels[sort_index[::-1]]
    # labels = np.array( labels[index[1:]].tolist() )

    binaryMap = toyplot.color.Palette( [toyplot.color.css("SlateGray"), toyplot.color.css("Wheat") ] )

    MLSTS = np.array( sorted( set( [ mlst_db[iso] for iso in isolates ]) ))
    cmap = matplotlib.cm.get_cmap('jet')
    tip_colors = { iso:cmap( (1.0* np.where(MLSTS == mlst_db[iso])[0][0] )/ (len(MLSTS)-1.0) ) for iso in isolates }
    tip_colors = { k:toyplot.color.rgba(v[0],v[1],v[2],v[3]) for k,v in tip_colors.items() }
    with open(folder + 'isoColors.pkl','wb+') as outPKL:
        pickle.dump(tip_colors,outPKL)

    cmap = [binaryMap]*(pa_matrix.shape[1]) 
    canvas = plotTree_w_attributes(tree,pa_matrix,tlabel=labels,tip_labels=tip_labels,
            tip_colors=tip_colors,colors=cmap,W=W,H=H,frac=frac,colGap=2,headSize=headSize,xdiff=xdiff,ydiff=ydiff)

    return canvas 

def plotTree_w_attributes(tree,matrix,tlabel,tip_labels,tip_colors,colors,W=800,H=400,M=50,frac=.23,rowH=.9,colGap=20,headSize=9,xdiff=.045,ydiff=.1):

    canvas = toyplot.Canvas(width=W,height=H)
    tree_axes = canvas.cartesian(bounds=(M,frac*(W-M),M,H-M),xlabel="Divergence")
    for node in tree.tree.traverse():
        if node.is_leaf():
            node.color= 'Grey'
        else:
            node.color= 'WhiteSmoke'
    nodeColors = tree.get_node_values('color', show_root=1, show_tips=1)
    tips = set(tree.get_tip_labels())
    tipColors = [  tip_colors[k] for k in tree.get_tip_labels() ] 

    _,tree_axes = tree.draw(axes=tree_axes,
                            node_labels=None,
                            node_size=None,
                            node_color=nodeColors,
                            node_style={"stroke":"black"},
                            use_edge_lengths=True,
                            tip_labels_align=True,
                            tip_labels=tip_labels,
                            padding=None,
                            show_tips=True,
                            show_root=False,
                            )

    tree_axes.y.show = False
    # tree_axes.x.show = False

    verts = tree.verts

        # Subsample verts for just leaves.
    # label = tree.get_node_values(feature='name', show_root=True, show_tips=True)
    # indx = np.array( [x in tips for x in label] )
    verts = verts[:len(tips),:]

    tree_axes.scatterplot(np.max(verts[:,0])*np.ones_like(verts[:,1])+xdiff,verts[:,1][::-1]+ydiff,color=tipColors,size=6)

    tree_ymin = tree_axes.__dict__['_ymin_range'] #+ tree_axes.__dict__['_padding']/2
    tree_ymax = tree_axes.__dict__['_ymax_range'] #- tree_axes.__dict__['_padding']/2

    # Linearly interpolate tree vertex position into tree domain
    vert_max = np.max(verts[:,1])
    vert_min = np.min(verts[:,1])
    
    vertY = (((tree_ymax - tree_ymin) / (vert_max - vert_min) ) * ( sorted(verts[:,1]) - vert_min ) ) + tree_ymin 
    rowH = rowH * np.mean( np.abs(np.diff(vertY)))
    topRow = vertY - rowH/2
    botRow = vertY + rowH/2
    
    tableMin = np.min(topRow) 
    tableMax = np.max(botRow) + headSize*rowH

    # Enlarge matrix.
    matrix = toyplot.require.scalar_matrix(matrix)

    # Build up the attribute table incrementally
    xmin_range, xmax_range, _, _ = toyplot.layout.region( 0, canvas._width, 0, canvas._height, 
            bounds=(frac*(W-M),W-M,M,H-M), rect=None, corner=None, grid=None, margin=None)
    xmin_range += 15
    table = toyplot.coordinates.Table(
                            rows=2*matrix.shape[0]-1,
                            columns=matrix.shape[1],
                            annotation=False,
                            brows=1,
                            trows=0,
                            rcolumns=0,
                            lcolumns=0,
                            xmax_range=xmax_range,
                            xmin_range=xmin_range,
                            ymax_range=tableMax,
                            ymin_range=tableMin,
                            parent=canvas,
                            label=None,
                            filename=None
                        )
    
    for n,bRow in enumerate(table.bottom.column):
        bRow.data = tlabel[n]
        bRow.height = headSize*rowH 
        bRow.angle = 90
        bRow.lstyle = {"font-size":10}

    # Color each cell.
    for i in np.arange(matrix.shape[0]):
        for j in np.arange(matrix.shape[1]):
            cell = table.body.cell[2*i, j]
            cell.style = {"stroke": "black", "stroke-width":.5, "fill": toyplot.color.to_css(colors[j].color(matrix[i,j]))}
            cell.title = matrix[i, j]

    # Change each cell's geometry.     
    for i in np.arange(2*matrix.shape[0]-1):
        for j in np.arange(matrix.shape[1]):
            cell = table.body.cell[i, j]
            if (i % 2 == 0):
                cell.height = rowH
            else:
                cell.height = topRow[int(i/2)+1] - botRow[int(i/2)] 

    table.body.gaps.columns[...] = colGap
    canvas._children.append(table)

    return canvas
