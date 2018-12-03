import glob, os
import subprocess

fasttree_aa_cmd = "FastTree "
fasttree_nc_cmd = "FastTree -gtr -nt "
raxml_aa_cmd = "raxml -m PROTGAMMAWAG -p12345 -T 8 -f d"
raxml_nc_cmd = "raxml -m GTRCAT -p12345 -T 8 -f d"

topFolder = "/scicore/home/neher/GROUP/data/Carbapenemases/blaTrees/molecules/"
familyFolders = glob.glob(topFolder + "bla*/")

for family in familyFolders:
    aaF = glob.glob(family+"aa/*_aln.faa")[0]
    ncF = glob.glob(family+"nc/*_aln.fna")[0]

    raxCMD = raxml_aa_cmd + " -w " + family + "aa/"
    if (os.stat(aaF).st_size > 0):
        nEntries = 0
        with open(aaF,"r") as aln_file:
            for line in aln_file.readlines():
                if (line[0] == ">"):
                    nEntries += 1

        if (nEntries > 3):
            rax_name = "aa_rax.newick"
            cmd = raxCMD + " -s " + aaF + " -n " + rax_name + " -c 25"
            print cmd
            process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
            process.wait()
        else:
            ft_name = "aa_ft.newick"
            cmd = fasttree_aa_cmd + aaF + " > " + ft_name
            print cmd
            process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
            process.wait()
            os.rename(ft_name, family + "aa/" + ft_name)

    raxCMD = raxml_nc_cmd + " -w " + family + "nc/"
    if (os.stat(ncF).st_size > 0):
        nEntries = 0
        with open(ncF,"r") as aln_file:
            for line in aln_file.readlines():
                if (line[0] == ">"):
                    nEntries += 1

        if (nEntries > 3):
            rax_name = "nc_rax.newick"
            cmd = raxCMD + " -s " + ncF + " -n " + rax_name + " -c 25"
            print cmd
            process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
            process.wait()
        else:
            ft_name = "nc_ft.newick"
            cmd = fasttree_nc_cmd + ncF + " > " + ft_name
            process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
            process.wait()
            os.rename(ft_name, family + "nc/" + ft_name)

