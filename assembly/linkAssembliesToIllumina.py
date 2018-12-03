#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 17:26:29 2017

@author: nolln

Submits an assembly script for all isolates.

"""

import os
import sys
import utils

dirName = sys.argv[1]
dirName = os.path.abspath(dirName)
if (not dirName.endswith('/')):
    dirName += '/'

dir2Name = os.path.abspath(dirName+'/../')
if (not dir2Name.endswith('/')):
    dir2Name += '/'
    
readNames = utils.obtainReadNames(dirName,True)
for n,readName in enumerate(readNames):
        strainName = readName[0].split(dirName)[1].split('_R1_')[0]
        fname1 = readName[0]
        fname2 = readName[1]
        
        target1 = dir2Name + 'assemblies/' + strainName + '/illumina_r1.fq.gz'
        target2 = dir2Name + 'assemblies/' + strainName + '/illumina_r2.fq.gz'

        cmd1 = "ln -s %s %s"%(fname1, target1)
        cmd2 = "ln -s %s %s"%(fname2, target2)

        os.system(cmd1)
        os.system(cmd2)