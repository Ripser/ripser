import numpy as np
import matplotlib.pyplot as plt
from ripser import Rips

def drawLineColored(X, C):
    plt.hold(True)
    for i in range(X.shape[0]-1):
        plt.plot(X[i:i+2, 0], X[i:i+2, 1], c=C[i, :], lineWidth = 3)

def plotCocycle2D(r, X, cocycle, thresh):
    """
    Given a 2D point cloud X, display a cocycle projected
    onto edges under a given threshold "thresh"
    """
    D = r.dm_
    plt.hold(True)
    #Plot all edges under the threshold
    N = X.shape[0]
    t = np.linspace(0, 1, 10)
    c = plt.get_cmap('Greys')
    C = c(np.array(np.round(np.linspace(0, 255, len(t))), dtype=np.int32))
    C = C[:, 0:3]

    for i in range(N):
        for j in range(N):
            if D[i, j] <= thresh:
                Y = np.zeros((len(t), 2))
                Y[:, 0] = X[i, 0] + t*(X[j, 0] - X[i, 0])
                Y[:, 1] = X[i, 1] + t*(X[j, 1] - X[i, 1])
                drawLineColored(Y, C)
    #Plot cocycle projected to edges under the chosen threshold
    for k in range(cocycle.shape[0]):
        [i, j, val] = cocycle[k, :]
        if D[i, j] <= thresh:
            [i, j] = [min(i, j), max(i, j)]
            a = 0.5*(X[i, :] + X[j, :])
            plt.text(a[0], a[1], '%g'%val, color='b')
    #Plot vertex labels
    for i in range(N):
        plt.text(X[i, 0], X[i, 1], '%i'%i, color='r')

if __name__ == '__main__':
    #Initialize circle
    plt.figure(figsize=(18, 6))

    """
    np.random.seed(5)
    n = 12
    """

    np.random.seed(6)
    n = 6
    
    t = np.linspace(0, 2*np.pi, n+1)[0:n]
    x = np.zeros((n, 2))
    x[:, 0] = np.cos(t)
    x[:, 1] = np.sin(t)
    x += 0.1*np.random.randn(n, 2)
    
    r = Rips(coeff=17, do_cocycles=True)
    diagrams = r.fit_transform(x)
    plt.subplot(131)
    r.plot(diagrams, show = False)
    
    #Now find the index of the maximum persistence point
    #in the diagram and highlight that point
    dgm1 = diagrams[1]
    idx = np.argmax(dgm1[:, 1] - dgm1[:, 0])
    plt.scatter(dgm1[idx, 0], dgm1[idx, 1], 60, 'k')
    plt.title("Max 1D birth = %.3g, death = %.3g"%(dgm1[idx, 0], dgm1[idx, 1]))
    
    #Now plot the representative cocycle associated with this point
    cocycle = r.cocycles_[1][idx]
    print(cocycle)
    plt.subplot(132)    
    thresh = dgm1[idx, 1]
    plotCocycle2D(r, x, cocycle, thresh)
    plt.title("1-Form Thresh=%.3g"%thresh)

    #Lower the threshold slightly below the cocycle birth time (cycle death time)
    plt.subplot(133)    
    #thresh = dgm1[idx, 0]+0.01
    thresh = dgm1[idx, 1]-0.02
    plotCocycle2D(r, x, cocycle, thresh)
    plt.title("1-Form Thresh=%.3g"%thresh)
    
    plt.savefig("cocycle.svg")
