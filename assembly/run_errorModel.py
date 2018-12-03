import glob 
import os 
from subprocess import call

topFolder = '/scicore/home/neher/GROUP/data/Carbapenemases/'
isoFolders = glob.glob( topFolder + "carb???/" )
submitScript = "/scicore/home/neher/nolln/genomeAssembly/submitScripts/errorModel.sh"

for iso in isoFolders:

    polishedRef = iso + "final/pilon/assembly.fasta"
    unpolishedRef = iso + "canu_nano/assembly.contigs.fasta"
    unpolishedFolder = iso + "canu_nano/"
    polishedFolder = iso + "final/pilon/"

    nano = iso + "nanopore_reads.fastq.gz"
    read1 = iso + "illumina_r1.fq.gz"
    read2 = iso + "illumina_r2.fq.gz"

    if (os.path.exists(polishedRef) and os.path.exists(unpolishedRef) and os.path.exists(nano) and os.path.exists(read1) and os.path.exists(read2)):
        print(iso)
        sbatchCall = "sbatch --export=polishedRef='%s',unpolishedRef='%s',nano='%s',read1='%s',read2='%s',polishedFolder='%s',unpolishedFolder='%s',isoFolder='%s' %s"%(polishedRef,unpolishedRef,nano,read1,read2,polishedFolder,unpolishedFolder,iso,submitScript)
        call(sbatchCall, shell=True)
        # break
