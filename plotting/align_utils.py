
"""
Created on Tue Nov  7 16:37:11 2017

@author: nolln

Contains utility functions to be used in aligning various features of the 
synteny graph. Useful for personal knowledge.

"""

import numpy as np

def seqToSeqEditDistance(x,y,scalar=False):
    D = np.zeros((len(x)+1,len(y)+1))
    D[0,:] = range(len(y)+1)
    D[:,0] = range(len(x)+1)
    
    for n in xrange(1,len(x)+1):
        for m in xrange(1,len(y)+1):
            distHor = D[n,m-1] + 1
            distVer = D[n-1,m] + 1
            if (x[n-1] == y[m-1]):
                distDiag = D[n-1,m-1]
            else: 
                distDiag = D[n-1,m-1] + 1
            D[n,m] = min(distHor,distVer,distDiag)
            
    if (scalar):
        return D[len(x),len(y)]
    else:
        return D

def seqToSeq(x,y):

    D = seqToSeqEditDistance(x,y)

    x = ['-'] + x
    y = ['-'] + y

    # Trace back our path
    N = len(x)-1
    M = len(y)-1
    eD = D[N,M]
    Yalign = []
    Xalign = []
    # nn = 0
    while (N > 0 or M > 0):

        delta = 1 - 1.0*(x[N] == y[M])

        if (N >= 1 and M >= 1):
            left = D[N-1,M]
            up = D[N,M-1]
            diag = D[N-1,M-1]
            # print delta
            # print left
            # print up
            # print diag
            if (D[N,M] - delta == diag):
                Yalign = [y[M]] + Yalign 
                Xalign = [x[N]] + Xalign 
                N -= 1
                M -= 1
            elif (D[N,M] - 1 == left):
                Yalign = ['-'] + Yalign 
                Xalign = [x[N]] + Xalign 
                N -= 1
            elif (D[N,M] - 1 == up):
                Yalign = [y[M]] + Yalign 
                Xalign = ['-'] + Xalign 
                M -= 1
        elif (M >= 1 and N == 0):
            Yalign = [y[M]] + Yalign 
            Xalign = ['-'] + Xalign 
            M -= 1
        else:
            Yalign = ['-'] + Yalign 
            Xalign = [x[N]] + Xalign
            N -= 1
        # print N
        # print M
        # nn += 1
        # if (nn > 1):
        #     break
    if (Xalign[0] == '-' and Yalign[0]=='-'):
        Xalign = Xalign[1:]
        Yalign = Yalign[1:]

    return Xalign,Yalign,eD

def loopDistance(A,B):
    D = ( seqToSeqEditDistance(A,B,scalar=True), seqToSeqEditDistance(A,B[::-1],scalar=True) )
    sigma = np.argmin(D)
    D = D[sigma]
    return 1. * D / max(len(A),len(B)), sigma

def matchLoops(loops):
    from munkres import Munkres   
    
    numLoops = [ len(loop) for loop in loops ] 
    maxInd = np.argmax(numLoops)
    
    matches = []
    sigma = []
    M = Munkres()
    for n in xrange(len(loops)):
        if (n != maxInd):
            eD = np.array( [ [loopDistance(loopRef,loop) for loop in loops[n] ] for loopRef in loops[maxInd] ] )
            matchIndex = M.compute(np.max(np.max(eD[:,:,0])) - eD[:,:,0])
            matches.append(np.array(matchIndex))
            sigma.append(np.array([eD[i[0],i[1],1] for i in matchIndex]))
        else:
            matches.append(np.array([(i,i) for i in xrange(len(loops[n])) ]))
            sigma.append(np.ones(numLoops[maxInd]))
    return matches,sigma,maxInd