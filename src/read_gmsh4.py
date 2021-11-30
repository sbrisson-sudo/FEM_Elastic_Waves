import matplotlib.pyplot as plt
import numpy as np

import gmsh

# definition d'un objet element

class Node:
    id = 1
    def __init__(self, x, y, region):
        self.id = Node.id
        Node.id += 1
        self.x = x
        self.y = y
        self.region = region
        
    def __repr__(self):
        return f"Node {self.id} [{self.x},{self.y}], region = {self.region}\n"

class Element:
    id = 1
    def __init__(self, nodes):
        self.nodes = nodes
        self.id = Element.id
        Element.id += 1

    def getCoords(self):
        return np.array([[n.x,n.y] for n in self.nodes])
    
    def __repr__(self):
        return f"Element {self.id}\nNodes : {[n.id for n in self.nodes]}\nCoords : {self.getCoords()}\n"

def plotMesh(elements, nodes, regions = []):
    
    fig, ax = plt.subplots()
    
    regionsSymbols = ["+b", "+g", "+r", "+y"]
    
    for e in elements:
        coords = e.getCoords()
        plt.fill(coords[:,0], coords[:,1], color="greenyellow", alpha=0.5, ec="k")
        
    for (i,r) in enumerate(regions):
        for n in nodes:
            if n.region == r: plt.plot(n.x, n.y, regionsSymbols[i])
            
    ax.set_aspect("equal", "box")
    plt.show()


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
    
    for i, (dim, tag) in enumerate(regions):
        nodeTagsRegion, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
        for tag in nodeTagsRegion:
            nodesList[tag].region = i+2
    
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
    
    file = "../meshes/octogon4.msh"
    
    regions = [
        (1, 21),
        (1, 23),
        (1, 22),
        (1, 24),
        ]
    
    elements, nodes = readGmsh4(file, regions)    
    plotMesh(elements, nodes, [2,3,4,5])
    

    
    