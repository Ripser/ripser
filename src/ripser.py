"""
Programmer: Chris Tralie (ctralie@alumni.princeton.edu)
Purpose: To wrap around Ripser to compute persistence diagrams
"""
import subprocess
import os
import numpy as np
import time
import matplotlib.pyplot as plt

def plotDGM(dgm, color = 'b', sz = 20, label = 'dgm', axcolor = np.array([0.0, 0.0, 0.0]), marker = None):
    if dgm.size == 0:
        return
    # Create Lists
    # set axis values
    axMin = np.min(dgm)
    axMax = np.max(dgm)
    axRange = axMax-axMin
    a = max(axMin - axRange/5, 0)
    b = axMax+axRange/5
    # plot line
    plt.plot([a, b], [a, b], c = axcolor, label = 'none')
    plt.hold(True)
    # plot points
    if marker:
        H = plt.scatter(dgm[:, 0], dgm[:, 1], sz, color, marker, label=label, edgecolor = 'none')
    else:
        H = plt.scatter(dgm[:, 0], dgm[:, 1], sz, color, label=label, edgecolor = 'none')
    # add labels
    plt.xlabel('Time of Birth')
    plt.ylabel('Time of Death')
    return H

def doRipsFiltrationDM(D, maxHomDim, thresh=-1, coeff=2):
    """
    Wrapper around Uli Bauer's Ripser code
    :param D: An NxN pairwise distance matrix
    :param maxHomDim: The dimension up to which to compute persistent homology
    :param thresh: Threshold up to which to add edges.  If not specified, add all
        edges up to the full clique
    :param coeff: A prime to use as the field coefficients for the PH computation

    :return: PDs (array of all persistence diagrams from 0D up to maxHomDim).
        Each persistence diagram is a numpy array
    """
    from pyRipser import doRipsFiltrationDM as DRFDM
    #import _ripser
    N = D.shape[0]
    if thresh == -1:
        thresh = np.max(D)*2
    [I, J] = np.meshgrid(np.arange(N), np.arange(N))
    DParam = np.array(D[I > J], dtype=np.float32)
    threshParam = np.array(thresh, dtype=np.float32)
    #res = _ripser.ripser(DParam, coeff, maxHomDim, threshParam)
    res = DRFDM(DParam, maxHomDim, thresh, coeff)
    PDs = []
    istart = 0
    for dim in range(maxHomDim+1):
        N = int(res[-1 - (maxHomDim-dim)])
        if dim > 0:
            N -= int(res[-2 - (maxHomDim-dim)])
        I = np.array(res[istart*2:(istart+N)*2])
        PDs.append(np.reshape(I, (N, 2)))
        istart += N
    return PDs

def doRipsFiltrationDMTextfile(D, maxHomDim, thresh=-1, coeff=2, getCocycles=False):
    """
    FOR DEBUGGING: Wrapper around the ripser executable with a
        text file and subprocess
    :param D: An NxN pairwise distance matrix
    :param maxHomDim: The dimension up to which to compute persistent homology
    :param thresh: Threshold up to which to add edges.  If not specified, add all
        edges up to the full clique
    :param coeff: A prime to use as the field coefficients for the PH computation
    :param getCocycles: True if cocycles should be computed and returned

    :return: PDs (array of all persistence diagrams from 0D up to maxHomDim).
        Each persistence diagram is a numpy array
        OR
        tuple (PDs, Cocycles) if returning cocycles
    """
    N = D.shape[0]
    #Step 1: Extract and output distance matrix
    fout = open("D.txt", "w")
    for i in range(0, N):
        for j in range(0, N):
            fout.write("%g "%D[i, j])
        if i < N-1:
            fout.write("\n")
    fout.close()

    #Step 2: Call ripser
    callThresh = 2*np.max(D)
    if thresh > 0:
        callThresh = thresh
    if getCocycles:
        proc = subprocess.Popen(["./ripser-representatives", "--format", "distance", "--dim", "%i"%maxHomDim, "--threshold", "%g"%callThresh, "--modulus", "%i"%coeff, "D.txt"], stdout=subprocess.PIPE)
    elif coeff > 2:
        proc = subprocess.Popen(["./ripser-coeff", "--format", "distance", "--dim", "%i"%maxHomDim, "--threshold", "%g"%callThresh, "--modulus", "%i"%coeff, "D.txt"], stdout=subprocess.PIPE)
    else:
        proc = subprocess.Popen(["./ripser", "--format", "distance", "--dim", "%i"%maxHomDim, "--threshold", "%g"%callThresh, "D.txt"], stdout=subprocess.PIPE)
    #stdout = proc.communicate()[0]
    PDs = []
    AllCocycles = []
    while True:
        output=proc.stdout.readline()
        if (output == b'' or output == '') and proc.poll() is not None:
            break
        if output:
            s = output.strip()
            if output[0:4] == b"dist":
                continue
            elif output[0:4] == b"valu":
                continue
            elif output[0:4] == b"pers":
                if len(PDs) > 0:
                    PDs[-1] = np.array(PDs[-1])
                PDs.append([])
                AllCocycles.append([])
            else:
                if getCocycles:
                    s = s.split(": ")
                    if len(s) > 1:
                        [s, s1] = s
                        c = parseCocycle(s1)
                        AllCocycles[-1].append(c)
                    else:
                        s = s[0]
                s = s.replace(b"[", b"")
                s = s.replace(b"]", b"")
                s = s.replace(b"(", b"")
                s = s.replace(b")", b"")
                s = s.replace(b" ", b"")
                fields = s.split(b",")
                b = float(fields[0])
                d = -1
                if len(fields[1]) > 0:
                    d = float(fields[1])
                PDs[-1].append([b, d])
        rc = proc.poll()
    PDs[-1] = np.array(PDs[-1])
    if getCocycles:
        return (PDs, AllCocycles)
    return PDs


def getSSM(X):
    """
    Given a set of Euclidean vectors, return a pairwise distance matrix
    :param X: An Nxd array of N Euclidean vectors in d dimensions
    :returns D: An NxN array of all pairwise distances
    """
    XSqr = np.sum(X**2, 1)
    D = XSqr[:, None] + XSqr[None, :] - 2 * X.dot(X.T)
    D[D < 0] = 0  # Numerical precision
    D = np.sqrt(D)
    return D


def doRipsFiltration(X, maxHomDim, thresh=-1, coeff=2):
    """
    Run ripser assuming Euclidean distance of a point cloud X
    :param X: An N x d dimensional point cloud
    :param maxHomDim: The dimension up to which to compute persistent homology
    :param thresh: Threshold up to which to add edges.  If not specified, add all
        edges up to the full clique
    :param coeff: A prime to use as the field coefficients for the PH computation

    :return: PDs (array of all persistence diagrams from 0D up to maxHomDim).
        Each persistence diagram is a numpy array
    """
    D = getSSM(X)
    return doRipsFiltrationDM(D, maxHomDim, thresh, coeff)


if __name__ == '__main__':
    """
    Test to make sure the numpy wrapper matches the output
    of the ripser exectuable
    """
    np.random.seed(10)
    X = np.random.randn(300, 2)
    X = X / np.sqrt(np.sum(X**2, 1)[:, None])
    D = getSSM(X)

    tic = time.time()
    PDs1 = doRipsFiltrationDMTextfile(D, 2, coeff=3)
    print("Elapsed Time Text File: %g"%(time.time() - tic))
    tic = time.time()
    PDs2 = doRipsFiltrationDM(D, 2, coeff=3)
    print("Elapsed Time Numpy Wrapper: %g"%(time.time() - tic))
    plt.subplot(141)
    plt.plot(X[:, 0], X[:, 1], '.')
    plt.subplot(142)
    plotDGM(PDs1[0][0:-1, :])
    plotDGM(PDs2[0][0:-1, :], color = 'r', marker = 'x')
    plt.title("H0")
    plt.subplot(143)
    plotDGM(PDs1[1])
    plotDGM(PDs2[1], color = 'r', marker = 'x')
    plt.title("H1")
    plt.subplot(144)
    plotDGM(PDs1[2])
    plotDGM(PDs2[2], color = 'r', marker = '*')
    plt.title("H2")
    plt.show()
