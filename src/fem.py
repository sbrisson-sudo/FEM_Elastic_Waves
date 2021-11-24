# IMPORTATIONS

import numpy as np
import matplotlib.pyplot as plt

# INTERPOLATIONS

class Interpolators:

    implemented = [
        "P1Quad",
    ]

    def __init__(self, itp):

        if itp == "P1Quad":
            self.N = 4
            self.phi = [
                    lambda u,v : 0.25*(1-u)*(1-v), # (-1,-1)
                    lambda u,v : 0.25*(1+u)*(1-v), # (+1,-1)
                    lambda u,v : 0.25*(1+u)*(1+v), # (-1,+1)
                    lambda u,v : 0.25*(1-u)*(1+v), # (+1,+1)
            ]
            self.duphi = [
                    lambda u,v : -0.25*(1-v),   # (-1,-1)
                    lambda u,v : 0.25*(1-v),    # (+1,-1)
                    lambda u,v : 0.25**(1+v),   # (-1,+1)
                    lambda u,v : -0.25*(1+v),   # (+1,+1)
            ]
            self.dvphi = [
                    lambda u,v : -0.25*(1-u),   # (-1,-1)
                    lambda u,v : -0.25*(1+u),   # (+1,-1)
                    lambda u,v : 0.25*(1+u),    # (-1,+1)
                    lambda u,v : 0.25*(1-u),    # (+1,+1)
            ]

    def interp(self, X):
        """Return the interpolation of the field with values X on nodes"""
        return lambda u,v : sum([phi(u,v)*x for phi,x in zip(self.phi, X)])

    def plot(self):

        if self.N == 4:
            U = np.linspace(-1,1,100)

            fig, axs = plt.subplots(1,4)
            levels = np.linspace(0,1,20)

            for (i,phi) in enumerate(self.phi):
                axs[i].contourf(U,U,np.array([[phi(u,v) for u in U] for v in U]),levels = levels)
                axs[i].set_aspect("equal")
                axs[i].set_title(f"$\Phi${i+1}")

        plt.show()

# QUADRATURE

def getGLLquad(N, dir="../quadratureRules/"):
    data = np.loadtxt(f"{dir}/gll_{N:02d}.tab")
    return data[0,:], data[1,:], data[2:,:] 

def getGLquad(N, dir="../quadratureRules/"):
    data = np.loadtxt(f"{dir}/gl_{N:02d}.tab")
    return data[:,0], data[:,1] 



# CLASSES ELEMENTS (P1)

class Elements:

    type = 3            
    table = []          # elements table
    nodesTable = []     # nodes table

    def __init__(self, nodes):

        self.nodes = nodes
        self.type = len(nodes)  # element type (3: triangles, 4: quadrangles)
        self.nodesCoord = np.array([Elements.nodesTable[n-1] for n in self.nodes])

        self.rho = 1000.0
        self.mu = 1e6

        self.itp = Interpolators("P1Quad")
        self.computeJacobian()

        self.id = len(Elements.table)
        Elements.table.append(self)

    def computeJacobian(self):
        X,Y = self.nodesCoord[:,0],self.nodesCoord[:,1]
        dxdu = lambda u,v : sum([self.itp.duphi(u,v)*x for phi,x in zip(self.phi, X)])
        dxdv = lambda u,v : sum([self.itp.dvphi(u,v)*x for phi,x in zip(self.phi, X)])        
        dydu = lambda u,v : sum([self.itp.duphi(u,v)*x for phi,x in zip(self.phi, Y)])
        dydv = lambda u,v : sum([self.itp.dvphi(u,v)*x for phi,x in zip(self.phi, Y)])
        self.detJ = lambda u,v : dxdu(u,v)*dydv(u,v) - dxdv(u,v)*dydu(u,v)

    def quadratureSetup(self, ui, wi):
        return

    def __repr__(self):
        out = f"Element {self.id}:\n"
        out += f"Type : {'Quadrangle' if self.type == 4 else 'Triangle'}\n"
        out += f"Nodes : {self.nodes}\n"
        return out

def readGmsh(meshFile):
    """
    Read a gmsh file (version 2) of 2D quadrangles.
    """

    with open(meshFile, "r") as f:

        lines = f.readlines()
        nodesStart, nodesEnd, elemStart, elemEnd = 0,0,0,0

        for (i,l) in enumerate(lines):
            if l == "$Nodes\n": nodesStart = i+2
            elif l == "$EndNodes\n": nodesEnd = i
            elif l == "$Elements\n": elemStart = i+2
            elif l == "$EndElements\n": elemEnd = i

        Elements.nodesTable = [list(map(float, l.split()[1:3])) for l in lines[nodesStart:nodesEnd]]

        for l in lines[elemStart:elemEnd]:
            Elements(list(map(int, l.split()))[5:])
            

def plotMesh():
    
    fig, ax = plt.subplots()

    for e in Elements.table:
        nodesCoord = e.getNodesCoord()
        plt.fill([n[0] for n in nodesCoord], [n[1] for n in nodesCoord], color="greenyellow", alpha=0.5, ec="k")

    ax.set_aspect("equal", "box")
    plt.show()
    

if __name__ == "__main__":
    
    #-------------
    # readGmsh("../meshes/octogon.msh")
    # plotMesh()
    # print(Elements.table[0])

    #-------------
    # Elements.nodesTable = np.array([
    #     [0.0, 0.0],
    #     [1.0, 0.0],
    #     [1.0, 1.0],
    #     [0.0, 1.0],
    # ])

    # e = Elements([1,2,3,4])
    # plotMesh()

    #-------------
    # itp = Interpolators("P1Quad")
    # itp.plot()

    #-------------
    # itp = Interpolators("P1Quad")
    # nodes = np.array([
    #     [0.0, 0.0],
    #     [2.0, 1.0],
    #     [3.0, 4.0],
    #     [1.0, 2.0],
    # ]) # nodes
    # T = [0.0, 120.0, 210.0, 30.0] # field values on nodes

    # Titp = itp.interp(T)

    # U = np.linspace(-1,1,100)
    
    # plt.contourf(U,U,np.array([[Titp(u,v) for u in U] for v in U]))
    # plt.colorbar()
    # plt.gca().set_aspect("equal")
    # plt.show()

    #-------------
    N = 5 # degré 
    wi, xi = getGLquad(N)
    print(xi)
    print(wi)

    wij = np.array([[wi[i]*wi[j] for i in range(N)] for j in range(N)])
    xij = np.array([[xi[i] for i in range(N)] for j in range(N)])

    linearize = lambda x : x.reshape(x.shape[0]*x.shape[1])

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(linearize(xij), linearize(xij.T), linearize(wij))
    plt.xlim([-1,1])
    plt.ylim([-1,1])
    plt.show()






