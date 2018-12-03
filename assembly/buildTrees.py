import os, sys, time, glob, subprocess, shutil
from ete2 import Tree

def resolve_polytomies(infileName, outfileName):
#    newickString=open(infileName, 'rb').readline().rstrip().replace('[&R] ', '')
    tree = Tree(infileName);
    tree.resolve_polytomy(recursive=True)
    with open(outfileName, 'wb') as outfile:
        outfile.write(tree.write(format=1))

def midpointRooting(infileName, outfileName):
    """ using ete2 for mid-point rooting """
    newickString=open(infileName, 'rb').readline().rstrip().replace('[&R] ', '')
    tree = Tree(newickString);
    
    if tree.get_midpoint_outgroup()!=None:
        tree.set_outgroup( tree.get_midpoint_outgroup() )
    tree.ladderize()
    with open(outfileName, 'wb') as outfile:
        outfile.write(tree.write(format=1))

def aln_to_Newick(alnPath, timeLimit, threads):
    """ function: build tree using SNP alignment
    """

    name = alnPath.split('.fna')[0]
    ft = 'FastTree' 
    ## run fasttree
    ftProcess = subprocess.Popen(ft +' -gtr -nt -gamma -mlacc 2 -slownni '
                                +alnPath+' > initial_tree.newick0 2> ' + name + '_ft.log',shell=True) ;
    ftProcess.communicate()
    resolve_polytomies('initial_tree.newick0','initial_tree.newick')

    ## run raxml
    outName = name + '.nwk'
    if (timeLimit>0):
        print '%s%d%s'%('RAxML tree optimization within the timelimit of ',timeLimit, ' minutes')
        # exec for killing process
        end_time = time.time() + int(timeLimit*60) #

        raxml_program = 'raxml' 
        process = subprocess.Popen('exec '+ raxml_program +' -f d -T '+ str(threads) +
                                  ' -j -s '+alnPath+' -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick > '
                                  +name+'_raxml.log', shell=True)
#        process.communicate()
        while (time.time() < end_time):
            if os.path.isfile('RAxML_result.topology'):
                break
            time.sleep(10)
        process.terminate()

        checkpoint_files = glob.glob('RAxML_checkpoint*')
        if os.path.isfile('RAxML_result.topology'):
            checkpoint_files.append('RAxML_result.topology')
        if (len(checkpoint_files) > 0):
            last_tree_file = checkpoint_files[-1]
            shutil.copy(last_tree_file, outName)
        else:
            shutil.copy('initial_tree.newick', outName)
    else:
        shutil.copy('initial_tree.newick', outName)

if (__name__ == "__main__") : 
    
    folder = sys.argv[1]
    if (not folder.endswith('/')):
        folder += '/'
    nThreads = sys.argv[2]
    raxTime = int(sys.argv[3])
    
    cwd = os.getcwd()
    os.chdir(folder)
    alnFiles = glob.glob('./*.fna')
    for aln in alnFiles:
        aln_to_Newick(aln,raxTime,nThreads)
        rmvFiles = glob.glob('*topology*')
        for rF in rmvFiles:
            os.remove(rF)
        if (os.path.exists(aln+'.reduced')):
            os.remove(aln+'.reduced')
        
    os.remove('initial_tree.newick')
    os.remove('initial_tree.newick0')
    os.chdir(cwd)
        