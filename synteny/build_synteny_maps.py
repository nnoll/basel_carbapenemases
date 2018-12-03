import glob
from geneSynteny import wordBlock

if __name__ == '__main__':
    geneFamilies = glob.glob("./molecules/*/")
    for gene in geneFamilies:
        syn = wordBlock(gene=gene,globS=False,overwrite=False)
        syn.all_to_all_circle_align(save=True,overwrite=False)
        syn.synteny_tree(overwrite=True)