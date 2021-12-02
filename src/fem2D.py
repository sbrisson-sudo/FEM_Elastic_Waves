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
    def __init__(self, nodes):
        self.nodes = nodes
        self.id = Element.id
        Element.id += 1

    def getCoords(self):
        return np.array([[n.x,n.y] for n in self.nodes[:4]])
    
    def __repr__(self):
        return f"Element {self.id}\nNodes : {[n.id for n in self.nodes]}\nCoords : {self.getCoords()}\n"
    
class FEdomain:
    
    def __init__(self, nodes, elements):
        self.nodes = deepcopy(nodes)
        self.elements = deepcopy(elements)
        
    def addInnerPoints(self, xi):
        """
        Add inner points in each elements, defined as space along 1D on [-1,1] as the x points
        """
        
        N = len(xi)
        
        for e in self.elements:
            
            coords = e.getCoords()
            X, Y = coords[:,0], coords[:,1]
            phiP1 = shapeFunctions["quad"]["P1"]["phi"]
            interpX = interp(phiP1, X)
            interpY = interp(phiP1, Y)
            
            n1 = e.nodes[0]
            
            for i,xi_1 in enumerate(xi):
                for j,xi_2 in enumerate(xi):
                    
                    if abs(xi_1) == 1 and abs(xi_2) == 1: continue
                    
                    print(f"[{i},{j}] add node")
                    
                    n = Node(interpX(xi_2,xi_1), interpY(xi_2,xi_1), 1)
                    self.nodes.append(n)
                    
                    if xi_1 == 1.0:
                        e.nodes.insert(i*(N+1)-j, n)
                    else:
                        e.nodes.insert(i*N+j, n)

#-----------------------------
# SHAPE FUNCTIONS

shapeFunctions = {
    "quad" : {
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
    data = np.loadtxt(f"../quadratureRules/gll_{N:02d}.tab")
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

    for i in range(len(elemtypes)):
        # lecture des quadrangles
        if elemtypes[i] == 3:
            for j,tag in enumerate(elemtags[i]):
                elements.append(Element([nodesList[n] for n in nodestags[i][4*j:4*(j+1)]]))

    gmsh.finalize()
            
    return elements, nodes


if __name__ == "__main__":
    
    #---------------------
    
    # file = "../meshes/octogon.msh"
    
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
        Node(1.0, 1.0, 1),
        Node(0.0, 1.0, 1),
    ]
    
    elements = [Element(nodes)]
    
    
    V = FEdomain(nodes, elements)
    
    N = 4

    xi, w ,_ = readGLL(N)
    
    print(V.elements[0])
    print(V.nodes)
    

    V.addInnerPoints(xi)
    
    print(V.elements[0])
    print(V.nodes)

    
    plotNodes(V.elements[0])
    

    #plotMesh(elements, nodes, [21,22,23,24], allNodes=True)

    
    