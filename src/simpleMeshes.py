import gmsh
import numpy as np
from fem2D import Node, Element, transform2GL

def squareMesh(L, dx, N = 1):
    """Build a square mesh, return the elements and nodes list, `L` square size, `dx` elements size, `N=1` elements order"""
    
    return rectangleMesh(L,L,dx,N)


def rectangleMesh(Lx, Ly, dx, N=1):
    
    if not(N in [1,2,4]): raise(Exception("Only orders 1,2 and 4 implemented."))

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 0)
    
    # defining geometry

    p1 = gmsh.model.geo.addPoint(0, 0, 0, dx)
    p2 = gmsh.model.geo.addPoint(Lx, 0, 0, dx)
    p3 = gmsh.model.geo.addPoint(Lx, Ly, 0, dx)
    p4 = gmsh.model.geo.addPoint(0, Ly, 0, dx)

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    outside = gmsh.model.geo.addCurveLoop([l1,l2,l3,l4])
    surface = gmsh.model.geo.addPlaneSurface([outside])

    gmsh.model.geo.mesh.setTransfiniteSurface(surface)
    gmsh.model.geo.mesh.setRecombine(2, surface)
    
    gmsh.option.setNumber("Mesh.ElementOrder", N)

    gmsh.model.geo.synchronize()

    bot,right,top,left = [2,3,4,5] # physical lines tag

    gmsh.model.addPhysicalGroup(1, [l1], bot)
    gmsh.model.addPhysicalGroup(1, [l2], right)
    gmsh.model.addPhysicalGroup(1, [l3], top)
    gmsh.model.addPhysicalGroup(1, [l4], left)
    
    # meshing

    gmsh.model.mesh.generate(2)
    
    # reading mesh
    
    regions = list(zip([1,1,1,1], [right,left,bot,top]))
    
    nodeTags, nodeCoords,_ = gmsh.model.mesh.getNodes()
    nodeCoords = nodeCoords.reshape((len(nodeTags), 3))
    
    nodes = [Node(x, y, 1) for x,y,_ in nodeCoords]
    nodesList = dict(zip(nodeTags, nodes))
        
    for dim_region, tag_region in regions:
        
        nodes_in_region, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim_region, tag_region)
        
        for n_id in nodes_in_region:
            nodesList[n_id].region = tag_region
    
    nodes = list(nodesList.values())
        
    elements = []
    elemtypes, elemtags, nodestags = gmsh.model.mesh.getElements()
        
    elmtType2order = {3:1, 10:2, 37:4}
    
    reOrder = {
        1 : np.array([0,1,3,2]),
        2 : np.array([0, 4, 1, 7, 8, 5, 3, 6, 2]),
        4 : np.array([0, 4, 5, 6, 1, 15, 16, 17, 18, 7, 14, 23, 24, 19, 8, 13, 22, 21, 20, 9, 3, 12, 11, 10, 2])
    }

    for i in range(len(elemtypes)):
        if elemtypes[i] in elmtType2order.keys():
            order = elmtType2order[elemtypes[i]]
            nnodes = (order+1)**2
            for j,tag in enumerate(elemtags[i]):
                elements.append(Element([nodesList[n] for n in nodestags[i][nnodes*j:nnodes*(j+1)][reOrder[order]]], order))

    gmsh.finalize()
    
    # Moving inner points to GL points
    
    N = int(np.sqrt(len(elements[0].nodes))) - 1
    if N > 1:
        transform2GL(nodes, elements)
                
    return elements, nodes

if __name__ == "__main__":
    
    from fem2D import plotMesh
    import matplotlib.pyplot as plt
    
    elements, nodes = rectangleMesh(200,100,10)
    bot,right,top,left = [2,3,4,5]
    
    plotMesh(elements, nodes, [bot,right,top,left])
    plt.show()