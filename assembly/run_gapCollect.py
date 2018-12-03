import glob 
import os 
from subprocess import call

topFolder = '/scicore/home/neher/GROUP/data/Carbapenemases/'
isoFolders = glob.glob( topFolder + "carb???/" )
submitScript = "/scicore/home/neher/nolln/genomeAssembly/submitScripts/nanoporeGaps.sh"

for iso in isoFolders:

    sam = iso + "final/pilon/minimap2/mappedNanopore.sam"

    if (os.path.exists(sam)):
        print(iso)
        sbatchCall = "sbatch --export=isoFolder='%s' %s"%(iso,submitScript)
        call(sbatchCall, shell=True)
        # break
