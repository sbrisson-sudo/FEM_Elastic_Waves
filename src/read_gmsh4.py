import matplotlib.pyplot as plt
import numpy as np

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
        return np.array([[n.x,n.y] for n in self.nodes])
    
    def __repr__(self):
        return f"Element {self.id}\nNodes : {[n.id for n in self.nodes]}\nCoords : {self.getCoords()}\n"

def plotMesh(elements, nodes, regions = [], allNodes=False):
    """
    plotMesh(elements, nodes, regions=[])
    """
    
    fig, ax = plt.subplots(figsize=(10,10))
    
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

def readGmsh4(file, regions = []):
    """ Read a gmsh mesh file version 4, regions : list of tuples (dim, tag) of physical groups to identify """
    
    gmsh.initialize()
    gmsh.open(file)
    
    # Lecture des noeuds
        
    nodeTags, nodeCoords,_ = gmsh.model.mesh.getNodes()
    nodeCoords = nodeCoords.reshape((len(nodeTags), 3))
    
    nodes = [Node(x, y, 0) for x,y,_ in nodeCoords]
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
    
    file = "../meshes/octogon.msh"
    
    regions = [
        (1, 21),
        (1, 23),
        (1, 22),
        (1, 24),
        ]
    
    elements, nodes = readGmsh4(file, regions)    
    plotMesh(elements, nodes, [21,22,23,24])