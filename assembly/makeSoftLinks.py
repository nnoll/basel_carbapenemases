import pandas
import glob
from collections import defaultdict
import os

speciesLUT = pandas.read_csv('../isolateAssemblies/speciesTable.csv')
assemblyFolders = glob.glob('../assemblies/*')

isolateNames = [name.split('.fasta')[0] for name in speciesLUT['genomeName']]
isolateSpecies = [name[0:4] for name in speciesLUT['organismName']]

isolate = defaultdict(list)
for n,name in enumerate(isolateNames):
    spName = isolateSpecies[n]
    if (spName == 'Shig'):
        spName = 'Esch'
    isolate[spName].append(name)

dirName = '..'    
for key in isolate.keys():
    strainFolder = dirName + '/panGenomes/' + key + '/'
    if (not os.path.exists(strainFolder)):
        os.mkdir(strainFolder)
        
    for n,individual in enumerate(isolate[key]):
        fname = os.path.abspath('../assemblies/' + individual) + '/prokka/PROKKA_11282017.gbk'
        target = os.path.abspath(strainFolder) + '/isolate_' + str(individual).replace('/','_').replace('-','_').replace(' ','_').replace('.','_') + '.gbk' # + '_%03d.gbk'%n
        cmd = "ln -s %s %s"%(fname, target)
        os.system(cmd)
#        print fname
#        print target