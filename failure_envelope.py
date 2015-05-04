import numpy as np
import scipy.optimize as opti

def loadedLaminate(a,b,N,M,analysis):
    Na = N * a
    Mb = M * b
    NM = np.matrix(np.zeros([6,1]))
    NM[0:3,0] = Na
    NM[3:,0] = Mb
    analysis.calculatePlyStressStrain(NM)

def abFailureIndex(a,b,N,M,analysis,failureType,kind='any'):
    loadedLaminate(a,b,N,M,analysis)
    index = failureType(analysis.laminate, kind)
    return 1-max(index)
    #print('a={a} b={b} i1={1} '.format(a=a,b=b,i=max(index)))

def baFailureIndex(b,a,N,M,analysis,failureType,kind='any'):
    # This just reverses the a and b input order. Sometimes useful to
    # fool the solver into solving in the other direction.
    loadedLaminate(a,b,N,M,analysis)
    index = failureType(analysis.laminate, kind)
    return 1-max(index)
    #print('a={a} b={b} i1={1} '.format(a=a,b=b,i=max(index)))

def makeNMLinVarEnvelope(N,M,analysis, failureType, kind='any'):
    solveArgs = (0,N,M,analysis,failureType,kind)
    aMin = opti.brentq(abFailureIndex,0,-100,args=solveArgs,xtol=1e-16)
    aMax = opti.brentq(abFailureIndex,0,100,args=solveArgs,xtol=1e-16)
    bMin = opti.brentq(baFailureIndex,0,-100,args=solveArgs,xtol=1e-16)
    bMax = opti.brentq(baFailureIndex,0,100,args=solveArgs,xtol=1e-16)


    print('{amin:.5g} <= a0 <= {amax:.5g}'.format(amin=aMin, amax=aMax))
    print('{bmin:.5g} <= b0 <= {bmax:.5g}'.format(bmin=bMin, bmax=bMax))

    aaRangeLeft = np.arange(2*aMin,0,0.01).tolist()
    aaRangeRight = np.arange(0,2*aMax,0.01).tolist()
    bUpper = list()
    bLower = list()
    aaRange = list()

    for a in aaRangeLeft:
        try:
            solveArgs = (a,N,M,analysis,failureType,kind)
            indexDown = opti.brentq(baFailureIndex,0,-100,args=solveArgs,xtol=1e-16)
            indexUp = opti.brentq(baFailureIndex,indexDown*0.99,100,args=solveArgs,xtol=1e-16)
            bLower.append(indexDown)
            bUpper.append(indexUp)
            aaRange.append(a)
            checkUp = abFailureIndex(a,indexUp,N,M,analysis,failureType,kind)
            checkDown = abFailureIndex(a,indexDown,N,M,analysis,failureType,kind)
            print('a={aa:.3g} {bLow:.3g} ({cLow:.3g}) < b < {bUp:.3g} ({cUp:.3g})'.format(aa=a,bLow=indexDown,bUp=indexUp,cLow=checkDown,cUp=checkUp))
        except ValueError:
            print('a={aa:.3g} No convergence'.format(aa=a))
    for a in aaRangeRight:
        try:
            solveArgs = (a,N,M,analysis,failureType,kind)
            indexUp = opti.brentq(baFailureIndex,0,100,args=solveArgs,xtol=1e-16)
            indexDown = opti.brentq(baFailureIndex,-100,indexUp*0.99,args=solveArgs,xtol=1e-16)
            bLower.append(indexDown)
            bUpper.append(indexUp)
            aaRange.append(a)
            checkUp = abFailureIndex(a,indexUp,N,M,analysis,failureType,kind)
            checkDown = abFailureIndex(a,indexDown,N,M,analysis,failureType,kind)
            print('a={aa:.3g} {bLow:.3g} ({cLow:.3g}) < b < {bUp:.3g} ({cUp:.3g})'.format(aa=a,bLow=indexDown,bUp=indexUp,cLow=checkDown,cUp=checkUp))
        except ValueError:
            print('a={aa:.3g} No convergence'.format(aa=a))

    aPlot = np.array([aMin] + aaRange + [aMax] + aaRange[::-1] + [aMin])
    bPlot = np.array([0] + bUpper + [0] + bLower[::-1] + [0])

    dictOut = {'aMin':aMin, 'aMax':aMax, 'bMin':bMin, 'bMax':bMax,\
               'aPlot':aPlot, 'bPlot':bPlot}
    return dictOut
