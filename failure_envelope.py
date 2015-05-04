import numpy as np
import scipy.optimize as opti

def loadedLaminate(a,b,analysis):
    N = a * np.matrix([[10.0e3],[5.0e3],[2.5e3]])
    M = b * np.matrix([[500],[700],[100]])
    NM = np.matrix(np.zeros([6,1]))
    NM[0:3,0]=N
    NM[3:,0]=M
    analysis.calculatePlyStressStrain(NM)

def abFailureIndex(a,b,analysis,failureType,kind='any'):
    loadedLaminate(a,b,analysis)
    index = failureType(analysis.laminate, kind)
    return 1-max(index)
    #print('a={a} b={b} i1={1} '.format(a=a,b=b,i=max(index)))

def baFailureIndex(b,a,analysis,failureType,kind='any'):
    # This just reverses the a and b input order. Sometimes useful to
    # fool the solver into solving in the other direction.
    loadedLaminate(a,b,analysis)
    index = failureType(analysis.laminate, kind)
    return 1-max(index)
    #print('a={a} b={b} i1={1} '.format(a=a,b=b,i=max(index)))

def makeLinVarEnvelope(analysis, failureType, kind='any'):
    aMin = opti.bisect(abFailureIndex,0,-100,args=(0,analysis,failureType,kind))
    aMax = opti.bisect(abFailureIndex,0,100,args=(0,analysis,failureType,kind))
    bMin = opti.bisect(baFailureIndex,0,-100,args=(0,analysis,failureType,kind))
    bMax = opti.bisect(baFailureIndex,0,100,args=(0,analysis,failureType,kind))


    print('{amin:.5g} <= a0 <= {amax:.5g}'.format(amin=aMin, amax=aMax))
    print('{bmin:.5g} <= b0 <= {bmax:.5g}'.format(bmin=bMin, bmax=bMax))

    aaRangeLeft = np.arange(2*aMin,0,0.01).tolist()
    aaRangeRight = np.arange(0,2*aMax,0.01).tolist()
    bUpper = list()
    bLower = list()
    aaRange = list()

    for a in aaRangeLeft:
        try:
            indexDown = opti.brentq(baFailureIndex,0,-100,args=(a,analysis,failureType,kind))
            indexUp = opti.brentq(baFailureIndex,indexDown*0.99,100,args=(a,analysis,failureType,kind))
            bLower.append(indexDown)
            bUpper.append(indexUp)
            aaRange.append(a)
            checkUp = abFailureIndex(a,indexUp,analysis,failureType,kind)
            checkDown = abFailureIndex(a,indexDown,analysis,failureType,kind)
            print('a={aa:.3g} {bLow:.3g} ({cLow:.3g}) < b < {bUp:.3g} ({cUp:.3g})'.format(aa=a,bLow=indexDown,bUp=indexUp,cLow=checkDown,cUp=checkUp))
        except ValueError:
            print('a={aa:.3g} No convergence'.format(aa=a))
    for a in aaRangeRight:
        try:
            indexUp = opti.brentq(baFailureIndex,0,100,args=(a,analysis,failureType,kind))
            indexDown = opti.brentq(baFailureIndex,-100,indexUp*0.99,args=(a,analysis,failureType,kind))
            bLower.append(indexDown)
            bUpper.append(indexUp)
            aaRange.append(a)
            checkUp = abFailureIndex(a,indexUp,analysis,failureType,kind)
            checkDown = abFailureIndex(a,indexDown,analysis,failureType,kind)
            print('a={aa:.3g} {bLow:.3g} ({cLow:.3g}) < b < {bUp:.3g} ({cUp:.3g})'.format(aa=a,bLow=indexDown,bUp=indexUp,cLow=checkDown,cUp=checkUp))
        except ValueError:
            print('a={aa:.3g} No convergence'.format(aa=a))

    dictOut = {'aMin':aMin, 'aMax':aMax, 'bMin':bMin, 'bMax':bMax,\
               'aaRange':aaRange, 'bUpper':bUpper, 'bLower':bLower}
    return dictOut
