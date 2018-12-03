import glob
import subprocess

topFolder = "/scicore/home/neher/GROUP/data/Carbapenemases/blaTrees/molecules/"
familyFolders = glob.glob(topFolder + "bla*/")
for family in familyFolders:
    faaFiles = glob.glob(family+ "aa/protein.faa")
    mafft_cmd = "mafft-linsi --amino "
    for fa in faaFiles:
        isoName = fa.rstrip(".faa")
        cmd = mafft_cmd + fa + " > " + isoName + "_aln.faa"
        process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        process.wait()
        print process.returncode

    fnaFiles = glob.glob(family + "nc/*.fna")
    mafft_cmd = "mafft-linsi --nuc "
    for fn in fnaFiles:
        isoName = fn.rstrip(".fna")
        cmd = mafft_cmd + fn + " > " + isoName + "_aln.fna"
        process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        process.wait()
        print process.returncode