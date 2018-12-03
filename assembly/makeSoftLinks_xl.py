import openpyxl as xlsx
import glob
import sys
import os        
from collections import defaultdict
import numpy as np

dirName = sys.argv[1]
if (dirName.endswith('/')):
    dirName = dirName[:-1]

assemblyName = dirName+'/assemblies'
folders = glob.glob(assemblyName +'/*')

strains = xlsx.load_workbook('20171114 Carba Liste.xlsx')
wb = strains.worksheets[0]
datCol = wb.dimensions.split(':')[1][0]
numOfStrains = int(wb.dimensions.split(':'+datCol)[1])

species = []
name = []
for n in xrange(numOfStrains):
    speciesName = str(wb['A'+str(n+1)].value.lower())
    
    # Take care of typos in the database file.
    if (speciesName == 'klpne'):
        speciesName = 'klepne'
    elif (speciesName == 'klepnc'):
        speciesName = 'klepne'
    elif (speciesName == 'klepnec'):
        speciesName = 'klepne'
    elif (speciesName == 'klepnee'):
        speciesName = 'klepne'
    elif (speciesName == 'esccole'):
        speciesName = 'esccol'
    elif (speciesName == 'klpene'):
        speciesName = 'klepne'
    elif (speciesName == 'escol'):
        speciesName = 'esccol'   
    elif (speciesName == 'pesaer'):
        speciesName= 'pseaer'
    elif (speciesName == 'acineto'):
        speciesName = 'acibau'
        
    species.append(speciesName)
    tmpName = str(wb['B'+str(n+1)].value) #.split('RV ')
    name.append(tmpName)
    
isolates = defaultdict(list)
unPaired = []
for ID,folder in enumerate(folders):
    isolateName = folder.split(assemblyName +'/')[1].split('_S')[0].split('-')
    if (len(isolateName) == 1):
        isolateName = isolateName[0]
    else:
        isolateName = isolateName[0] + '-' + isolateName[1]
    
    contained = [ isolateName in n for n in name ]
    if (any(contained)):
        ind = np.where(contained)[0][0]
        isolates[species[ind]].append([ID,name[ind]])
    else:
        # Deal with inconsistent number of dashes in names within the database
        isolateName = folder.split(assemblyName +'/')[1].split('_S')[0].split('-')
        isolateName = isolateName[0]
        contained = [ isolateName in n for n in name ]
        if (any(contained)):
            ind = np.where(contained)[0][0]
            isolates[species[ind]].append([ID,name[ind]])
                    
for key in isolates.keys():
    strainFolder = dirName + '/panGenomes/' + key + '/'
    if (not os.path.exists(strainFolder)):
        os.mkdir(strainFolder)
        
    for n,individual in enumerate(isolates[key]):
        fname = os.path.abspath(folders[individual[0]]) + '/prokka/PROKKA_11272017.gbk'
        target = os.path.abspath(strainFolder) + '/isolate_' + str(individual[1]).replace('/','_').replace('-','_').replace(' ','_').replace('.','_') + '.gbk' # + '_%03d.gbk'%n
        cmd = "ln -s %s %s"%(fname, target)
        os.system(cmd)
    