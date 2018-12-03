import glob, os
import subprocess

fasttree_aa_cmd = "FastTree "
fasttree_nc_cmd = "FastTree -gtr -nt "
raxml_aa_cmd = "raxml -m PROTGAMMAWAG -p12345 -T 8 -f d"
raxml_nc_cmd = "raxml -m GTRCAT -p12345 -T 8 -f d"

topFolder = "/scicore/home/neher/GROUP/data/Carbapenemases/blaTrees/molecules/"
familyFolders = glob.glob(topFolder + "bla*/")

for family in familyFolders:
    name = family.lstrip(topFolder).strip("/")
    if (not os.path.exists(family + "input_GenBank") or not os.path.exists(family + "protein_faa/diamond_matches/allclusters.tsv")):
        # print family
        sbatchCMD = "sbatch --export=folder='%s',name=%s panX_geneFamily.sh"%(family,name)
        # print sbatchCMD
        os.system(sbatchCMD)
        # break
