import glob
import subprocess

faaFiles = glob.glob("fasta/*.faa")
mafft_cmd = "mafft-linsi --amino "
for fa in faaFiles:
    isoName = fa.rstrip(".faa").replace("fasta/","aln/")
    cmd = mafft_cmd + fa + " > " + isoName + "_aln.faa"
    process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    process.wait()
    print process.returncode

fnaFiles = glob.glob("fasta/*.fna")
mafft_cmd = "mafft-linsi --nuc "
for fn in fnaFiles:
    isoName = fn.rstrip(".fna").replace("fasta/","aln/")
    cmd = mafft_cmd + fn + " > " + isoName + "_aln.fna"
    process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    process.wait()
    print process.returncode