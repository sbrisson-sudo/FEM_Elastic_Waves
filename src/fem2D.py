import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

import gmsh

# definition d'un objet element

class Node:
    id = 0
    def __init__(self, x, y, region):
        self.id = Node.id
        Node.id += 1
        self.x = x
        self.y = y
        self.region = region
        
    def __repr__(self):
        return f"Node {self.id} [{self.x},{self.y}], region = {self.region}\n"

class Element:
    id = 0
    def __init__(self, nodes, N=1):
        """N : order of element (1 -> linear)"""
        self.nodes = nodes
        self.N = N
        self.id = Element.id
        Element.id += 1

    def getCoords(self):
        # if self.N == 1:
        #     return np.array([[n.x,n.y] for n in self.nodes])
        return np.array([[n.x,n.y] for n in [self.nodes[0], 
                                             self.nodes[self.N],
                                             self.nodes[(self.N+1)**2-1],
                                             self.nodes[(self.N+1)**2-1-self.N],
                                            ]])
    
    def __repr__(self):
        return f"Element {self.id}\nNodes : {[n.id for n in self.nodes]}\nCoords : {self.getCoords()}\n"


#-----------------------------
# SHAPE FUNCTIONS

shapeFunctions = {
    "P1" : {
        "phi" : [
    lambda u,v : 0.25*(1-u)*(1-v), # (-1,-1)
    lambda u,v : 0.25*(1+u)*(1-v), # (+1,-1)
    lambda u,v : 0.25*(1+u)*(1+v), # (+1,+1)
    lambda u,v : 0.25*(1-u)*(1+v), # (-1,+1)
],
        "duphi" : [
    lambda u,v : -0.25*(1-v),   # (-1,-1)
    lambda u,v : 0.25*(1-v),    # (+1,-1)
    lambda u,v : 0.25*(1+v),   # (-1,+1)
    lambda u,v : -0.25*(1+v),   # (+1,+1)
],
        "dvphi" : [
    lambda u,v : -0.25*(1-u),   # (-1,-1)
    lambda u,v : -0.25*(1+u),   # (+1,-1)
    lambda u,v : 0.25*(1+u),    # (-1,+1)
    lambda u,v : 0.25*(1-u),    # (+1,+1)
],
    }

}   

def interp(phi, X):
    return lambda u,v : sum([f(u,v)*x for f,x in zip(phi, X)])

#-----------------------------
# LECTURE REGLES DE QUADRATURE

def readGL(N):
    data = np.loadtxt(f"../quadratureRules/gl_{N:02d}.tab")
    w, xi = data[:,0], data[:,1]
    return xi, w

def readGLL(N):
    data = np.loadtxt(f"../quadratureRules/gll_{N+1:02d}.tab")
    xi, w, l = data[0,:], data[1,:], data[2:,:]
    return xi, w, l

def getInt2D(xi,w):
    """Return a function computing the integral of a function f over [-1,1]^2"""
    return lambda f: sum([ sum([w1*w2 * f(xi1, xi2) for xi1,w1 in zip(xi, w)]) for xi2,w2 in zip(xi, w)])

#-----------------------------
# FONCTIONS DE PLOT

def plotMesh(elements, nodes, regions = [], allNodes=False):
    """
    plotMesh(elements, nodes, regions=[])
    """
    
    fig, ax = plt.subplots()
    
    regionsSymbols = ["+b", "+g", "+r", "+m"]
    
    for e in elements:
        coords = e.getCoords()
        plt.fill(coords[:,0], coords[:,1], color="greenyellow", alpha=0.5, ec="k")
    
    if allNodes:
        for n in nodes:
            plt.plot(n.x, n.y, "k+")
    
    for (i,r) in enumerate(regions):
        first = True
        for n in nodes:
            if n.region == r: 
                plt.plot(n.x, n.y, regionsSymbols[i%4],
                         label = r if first else "")
                first = False
    
    ax.legend()
    ax.set_aspect("equal", "box")
    plt.show()
    
def plotNodes(e):
    for i,n in enumerate(e.nodes):
        plt.plot(n.x, n.y, "k+")
        plt.annotate(str(i), (n.x, n.y))
    plt.gca().set_aspect("equal", "box")
    plt.show()

#-----------------------------
# LECTURE MAILLAGES GMSH


def readGmsh4(file, regions = []):
    """ Read a gmsh mesh file version 4, regions : list of tuples (dim, tag) of physical groups to identify """
    
    gmsh.initialize()
    gmsh.open(file)
    
    # Lecture des noeuds
    
        
    nodeTags, nodeCoords,_ = gmsh.model.mesh.getNodes()
    nodeCoords = nodeCoords.reshape((len(nodeTags), 3))
    
    nodes = [Node(x, y, 1) for x,y,_ in nodeCoords]
    nodesList = dict(zip(nodeTags, nodes))
    
    # Lecture des regions
    
    for dim_region, tag_region in regions:
        nodes_in_region, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim_region, tag_region)
        for n_id in nodes_in_region:
            nodesList[n_id].region = tag_region
    
    nodes = list(nodesList.values())
    
    # Lecture des éléments
    
    elements = []
    elemtypes, elemtags, nodestags = gmsh.model.mesh.getElements()
    
    elmtType2order = {3:1, 10:2, 37:4}
    
    reOrder = {
        1 : np.array([0,1,3,2]),
        2 : np.array([0, 5, 1, 7, 8, 5, 3, 6, 2]),
        4 : np.array([0, 4, 5, 6, 1, 15, 16, 17, 18, 7, 14, 23, 24, 19, 8, 13, 22, 21, 20, 9, 3, 12, 11, 10, 2])
    }

    for i in range(len(elemtypes)):
        # lecture des quadrangles
        if elemtypes[i] in elmtType2order.keys():
            
            order = elmtType2order[elemtypes[i]]
            nnodes = (order+1)**2
            
            for j,tag in enumerate(elemtags[i]):
                elements.append(Element([nodesList[n] for n in nodestags[i][nnodes*j:nnodes*(j+1)][reOrder[order]]], order))

    gmsh.finalize()
                
    return elements, nodes

def transform2GL(nodes, elements):
    """Move the inner nodes to the Gauss Lobatto points"""
    
    phi = shapeFunctions["P1"]["phi"]
    
    N = int(np.sqrt(len(elements[0].nodes))) - 1
    xi,_,_ = readGLL(N)
        
    for e in elements:
        
        coords = e.getCoords()
        X = interp(phi, coords[:,0])
        Y = interp(phi, coords[:,1])
        
        for k,n in enumerate(e.nodes):
            i = k//(N+1)
            j = k - i*(N+1)
            n.x = X(xi[j], xi[i])
            n.y = Y(xi[j], xi[i])


if __name__ == "__main__":
    
    #---------------------
    
    # file = "../meshes/octogon2.msh"
    
    # regions = [
    #     (1, 21),
    #     (1, 23),
    #     (1, 22),
    #     (1, 24),
    #     ]
    
    # elements, nodes = readGmsh4(file, regions)    
    
    #---------------------
    
    nodes = [   
        Node(0.0, 0.0, 1),
        Node(1.0, 0.0, 1),
        Node(0.0, 1.0, 1),
        Node(1.0, 1.0, 1),
    ]
    
    elements = [Element(nodes)]

    #---------------------
    
    N = int(np.sqrt(len(elements[0].nodes))) - 1
    
    print(f"N = {N}")
    
    if N > 1:
        transform2GL(nodes, elements)
    
    plotMesh(elements, nodes, [21,22,23,24], allNodes=True)