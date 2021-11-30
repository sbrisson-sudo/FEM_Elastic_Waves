# Importations

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import gmsh

# Classes

# definition d'un objet element

class Mesh:
    
    def __init__(self, nodes, elements, regions):
        self.nodes = nodes
        self.elements = elements
        self.regions = regions
        

class FEequation:
    
    def __init__():
        return

class FEproblem:
    
    def __init__(self, meshNodes, meshElements, SF, quadrature, equation):
        return

class Node:
    """Node class"""
    id = 0
    def __init__(self, x, y, region):
        """
        Node(x, y, region)
        """
        self.id = Node.id
        Node.id += 1
        self.x = x
        self.y = y
        self.region = region
        
    def __repr__(self):
        return f"Node {self.id} [{self.x},{self.y}], region = {self.region}"

class Element:
    id = 0
    def __init__(self, nodes):
        """
        Element(nodes)
        """
        self.nodes = nodes
        self.nnodes = len(nodes)
        self.id = Element.id
        Element.id += 1
    
    def setShapeFunction(self, SF):
        self.SF = SF    
        
    def computeJacobian(self):
        if not(self.SF): raise("Element does not have a shape function assigned")
        
        coords = self.getCoords()
        X, Y = coords[:,0], coords[:,1]

        dxdu = ShapeFunctions.interp(self.SF.duphi, X)
        dxdv = ShapeFunctions.interp(self.SF.dvphi, X)        
        dydu = ShapeFunctions.interp(self.SF.duphi, Y)
        dydv = ShapeFunctions.interp(self.SF.dvphi, Y)
        
        self.detJ = lambda u,v : dxdu(u,v)*dydv(u,v) - dxdv(u,v)*dydu(u,v)

        self.iJ = lambda u,v :  np.array(
            [[dydv(u,v), -dydu(u,v)],
            [-dxdv(u,v), dxdu(u,v)]]
        )/self.detJ(u,v)
        
    def computeKe(self, quadrature):
        
        self.Ke = np.zeros((self.nnodes,self.nnodes))
    
        for i in range(self.nnodes):
            for j in range(self.nnodes):

                fk = lambda u,v: (self.iJ(u,v)@self.SF.gradPhi[i](u,v))@(self.iJ(u,v)@self.SF.gradPhi[j](u,v))*self.detJ(u,v)
                self.Ke[i,j] = quadrature.int2D(fk)
        
    def getCoords(self):
        return np.array([[n.x,n.y] for n in self.nodes])
    
    def __repr__(self):
        return f"Element {self.id}\nNodes : {[n.id for n in self.nodes]}\nCoords : {self.getCoords()}"
    
    
    
class ShapeFunctions:
    
    def __init__(self, nameSF, typeElement):
        
        global shapeFunctions
        SF = shapeFunctions.get(typeElement).get(nameSF)
        
        self.phi = SF.get("phi")
        self.duphi = SF.get("duphi")
        self.dvphi = SF.get("dvphi")
        self.gradPhi = [
            lambda u,v: np.array([self.duphi[0](u,v), self.dvphi[0](u,v)]),
            lambda u,v: np.array([self.duphi[1](u,v), self.dvphi[1](u,v)]),
            lambda u,v: np.array([self.duphi[2](u,v), self.dvphi[2](u,v)]),
            lambda u,v: np.array([self.duphi[3](u,v), self.dvphi[3](u,v)]),
        ]
            
def interp(phi, X):
    return lambda u,v : sum([f(u,v)*x for f,x in zip(phi, X)])
    
        
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

class QuadratureRule:
    
    def __init__(self, name, order):
        if not(name in ["GL", "GLL"]): raise(f"{name} : invalid quadrature rule.")
        
        try: data = np.loadtxt(f"../quadratureRules/gl_{order:02d}.tab")
        except FileNotFoundError: raise(f"{name} quadrature at order {order} not found.")
        
        self.w, self.xi = data[:,0], data[:,1] 
    
    def int2D(self, f):
        return sum([ sum([w1*w2 * f(xi1, xi2) for xi1,w1 in zip(self.xi, self.w)]) for xi2,w2 in zip(self.xi, self.w)])
    

def plotMesh(elements, nodes, regions = []):
    """
    plotMesh(elements, nodes, regions=[])
    """
    
    fig, ax = plt.subplots()
    
    regionsSymbols = ["+b", "+g", "+r", "+m"]
    
    for e in elements:
        coords = e.getCoords()
        plt.fill(coords[:,0], coords[:,1], color="greenyellow", alpha=0.5, ec="k")
    
    
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


def readGmsh4(file, regions = []):
    """ 
    readGmsh4(file, regions=[])
    Read a gmsh mesh file version 4.
    - `regions` : dict of keys tuples (dim, tag) of physical groups to
    identify, and values the name of the regions. If an element is part
    of several physical groups, the last specify will be retains.
     """
    
    gmsh.initialize() # régler le niveau de verbosité
    gmsh.open(file)
    
    # Lecture des noeuds
        
    nodeTags, nodeCoords,_ = gmsh.model.mesh.getNodes()
    nodeCoords = nodeCoords.reshape((len(nodeTags), 3))
    
    nodes = [Node(x, y, "main") for x,y,_ in nodeCoords]
    nodesList = dict(zip(nodeTags, nodes))
    
    # Lecture des regions
    
    for i, (dim, tag) in enumerate(regions.keys()):
        nodeTagsRegion, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
        for tag in nodeTagsRegion:
            nodesList[tag].region = list(regions.values())[i]
    
    nodes = nodesList.values()
    
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
    
    file = "../meshes/octogon.msh"
    
    regions = {
        (1, 21) : "bottom",
        (1, 23) : "top",
        (1, 22) : "right",
        (1, 24) : "left"
        }
    
    elements, nodes = readGmsh4(file, regions)    
    plotMesh(elements, nodes, regions.values())
    
