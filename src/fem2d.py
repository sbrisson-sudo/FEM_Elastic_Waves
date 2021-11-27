# Importations

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import gmsh

# Classes

# definition d'un objet element

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
        self.id = Element.id
        Element.id += 1
        
    def getCoords(self):
        return np.array([[n.x,n.y] for n in self.nodes])
    
    def __repr__(self):
        return f"Element {self.id}\nNodes : {[n.id for n in self.nodes]}\nCoords : {self.getCoords()}"
    
    

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
    
