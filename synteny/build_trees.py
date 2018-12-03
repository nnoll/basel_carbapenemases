import glob, os
import subprocess

aa_alnFiles = glob.glob("aln/*_aln.faa")
nc_alnFiles = glob.glob("aln/*_aln.fna")

fasttree_aa_cmd = "FastTree "
fasttree_nc_cmd = "FastTree -gtr -nt "
raxml_aa_cmd = "raxml -m PROTGAMMAWAG -w /scicore/home/neher/GROUP/data/Carbapenemases/blaTrees/trees -p12345 -T 8 -f d"
raxml_nc_cmd = "raxml -m GTRCAT -w /scicore/home/neher/GROUP/data/Carbapenemases/blaTrees/trees -p12345 -T 8 -f d"

for aln in aa_alnFiles:
    if (os.stat(aln).st_size > 0):

        nEntries = 0
        with open(aln,"r") as aln_file:
            for line in aln_file.readlines():
                if (line[0] == ">"):
                    nEntries += 1

        isoName = aln.rstrip("_aln.faa")
        if (nEntries > 3):
            rax_name = isoName.replace("aln/","") + "_aa_rax.newick"
            cmd2 = raxml_aa_cmd + " -s " + aln + " -n " + rax_name + " -c 25"
            # print cmd2
            process2 = subprocess.Popen(cmd2,shell=True,stdout=subprocess.PIPE)
            process2.wait()
            # break
        else:
            ft_name = isoName.replace("aln/","") + "_aa_ft.newick"
            cmd1 = fasttree_aa_cmd + aln + " > " + ft_name
            process = subprocess.Popen(cmd1,shell=True,stdout=subprocess.PIPE)
            process.wait()

    # break
    # print process.returncode
    # ft_name = isoName + "_ft.newick"
    # cmd1 = fasttree_cmd + aln + " > " + ft_name
    # process = subprocess.Popen(cmd1,shell=True,stdout=subprocess.PIPE)
    # process.wait()

for aln in nc_alnFiles:
    if (os.stat(aln).st_size > 0):
        
        nEntries = 0
        with open(aln,"r") as aln_file:
            for line in aln_file.readlines():
                if (line[0] == ">"):
                    nEntries += 1

        isoName = aln.rstrip("_aln.fna")
        if (nEntries > 3):
            rax_name = isoName.replace("aln/","") + "_nc_rax.newick"
            cmd2 = raxml_nc_cmd + " -s " + aln + " -n " + rax_name + " -c 25"
            # print cmd2
            process2 = subprocess.Popen(cmd2,shell=True,stdout=subprocess.PIPE)
            process2.wait()
        else:
            ft_name = isoName.replace("aln/","") + "_nc_ft.newick"
            cmd1 = fasttree_nc_cmd + aln + " > " + ft_name
            process = subprocess.Popen(cmd1,shell=True,stdout=subprocess.PIPE)
            process.wait()
