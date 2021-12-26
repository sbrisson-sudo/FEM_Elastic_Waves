import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from functools import reduce
from scipy.spatial import ConvexHull
import scipy.interpolate as itp 
from matplotlib.path import Path

import os

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

def eraseCurrentMesh():
    Element.id = 0
    Node.id = 0

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

# lagrangian interpolants
def lagrangian2D(N):
    
    xi_l = list(readGLL(N)[0])

    # passage des paramètres par arguments nécessaire
    h = [
        lambda x,xi1=xi1b,i=ib : reduce(lambda a,b:a*b, [(lambda xi2=xi2b : (x-xi2)/(xi1-xi2))(xi2b) for xi2b in xi_l[:i]+xi_l[i+1:]])
        for (ib,xi1b) in enumerate(xi_l)
    ]
    hij = []

    for hib in h:
        for hjb in h:
            hij.append(lambda u,v,hi=hib,hj=hjb : hi(v)*hj(u))
            
    return hij
        
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
        
    fig, ax = plt.subplots(figsize = (6,6))
    
    regionsSymbols = ["+b", "+g", "+r", "+m"]
    
    for e in elements:
        coords = e.getCoords()
        ax.fill(coords[:,0], coords[:,1], color="greenyellow", alpha=0.5, ec="k")
    
    if allNodes:
        for n in nodes:
            ax.plot(n.x, n.y, "k+")
    
    for (i,r) in enumerate(regions):
        first = True
        for n in nodes:
            if n.region == r: 
                ax.plot(n.x, n.y, regionsSymbols[i%4],
                         label = r if first else "")
                first = False
    
    ax.legend()
    ax.set_aspect("equal", "box")
    return ax
    
def plotNodes(e):
    for i,n in enumerate(e.nodes):
        plt.plot(n.x, n.y, "k+")
        plt.annotate(str(i), (n.x, n.y))
    plt.gca().set_aspect("equal", "box")
    plt.show()
    
def plotNodesValues(U, nodes, ax):
    dots = ax.scatter([n.x for n in nodes], [n.y for n in nodes], c=U, zorder=5)
    plt.colorbar(dots)
    ax.set_aspect("equal", "box")
    
def plotMeshLimits(elements, ax):
    for e in elements:
        coords = e.getCoords()
        ax.fill(coords[:,0], coords[:,1], facecolor="none", ec="k", zorder=1)
        
def plotDeformedMesh(elements, U, ax, s = 1.0):
    """plot deformed mesh, U being a 2D displacement field"""
    for e in elements:
        coords = e.getCoords()

        nodesVertices = [
            e.nodes[0], 
            e.nodes[e.N],
            e.nodes[(e.N+1)**2-1],
            e.nodes[(e.N+1)**2-1-e.N],
        ]
        u1 = np.array([U[2*n.id] for n in nodesVertices])
        u2 = np.array([U[2*n.id+1] for n in nodesVertices])
        ax.fill(coords[:,0]+s*u1, coords[:,1]+s*u2, facecolor="none", ec="k", zorder=1)

#-----------------------------
# LECTURE MAILLAGES GMSH


def readGmsh4(file, regions = []):
    """ Read a gmsh mesh file version 4, regions : list of tuples (dim, tag) of physical groups to identify """
    
    if not(os.path.exists(file)): raise Exception(f"mesh file {file} not found.")
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 0)
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
        2 : np.array([0, 4, 1, 7, 8, 5, 3, 6, 2]),
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
    
    if not(elements): raise(Exception("No quadrangles of degree 1,2 or 4 found."))
    
    # Moving inner points to GL points
    
    N = int(np.sqrt(len(elements[0].nodes))) - 1
    if N > 1:
        transform2GL(nodes, elements)
                
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
            
def getCenterNode(nodes):
    """get the node at the center of the mesh"""
    
    coords = np.array([[n.x,n.y] for n in nodes])
    meanX, meanY = coords[:,0].mean(), coords[:,1].mean()
    
    distToCenter = (coords[:,0] - meanX)**2 + (coords[:,1] - meanY)**2  
    
    return np.argmin(distToCenter)
            
#-------------------------
# Post-Processing

def outputGrid(N,nodes):
    """Return a regular spaced grid and a mask on the mesh"""
    nodesCoord = np.array([[n.x,n.y] for n in nodes])

    xmin,xmax = nodesCoord[:,0].min(),nodesCoord[:,0].max()
    ymin,ymax = nodesCoord[:,1].min(),nodesCoord[:,1].max()

    NptsGrid = 200

    gridx = np.linspace(xmin,xmax,NptsGrid)
    gridy = np.linspace(ymin,ymax,NptsGrid)

    maskMesh = np.full((NptsGrid,NptsGrid), True)

    gridxx, gridyy = np.meshgrid(gridx,gridy)
    
    hull = ConvexHull(nodesCoord)
    hull_path = Path(nodesCoord[hull.vertices])

    for i in range(NptsGrid):
        for j in range(NptsGrid):
            if hull_path.contains_point((gridxx[i,j],gridyy[i,j])): maskMesh[i,j] = False
    
    return gridx,gridy,maskMesh

def interpOnGrid(U,elements,gridx,gridy,maskMesh):
    """interpolate the value of the field U over the grid"""
    
    Npts = 50

    u = np.linspace(-1, 1, Npts)
    uu, vv = np.meshgrid(u,u)

    vmin, vmax = U.min(), U.max()

    X,Y,Z = np.zeros(Npts**2*len(elements)),np.zeros(Npts**2*len(elements)),np.zeros(Npts**2*len(elements))

    # interpolants
    phi = shapeFunctions["P1"]["phi"]
    N = int(np.sqrt(len(elements[0].nodes))) - 1
    hij = lagrangian2D(N)
    
    for i,e in enumerate(elements):
        
        coords = e.getCoords()
        
        zitp = interp(hij, [U[n.id] for n in e.nodes])
        xitp = interp(phi, coords[:,0])
        yitp = interp(phi, coords[:,1])
        
        X[i*Npts**2:(i+1)*Npts**2] = xitp(uu, vv).reshape(Npts**2)
        Y[i*Npts**2:(i+1)*Npts**2] = yitp(uu, vv).reshape(Npts**2)
        Z[i*Npts**2:(i+1)*Npts**2] = zitp(uu, vv).reshape(Npts**2)
        
    points = list(zip(X,Y))
        
    # interpolation entre ces points  
    gridxx, gridyy = np.meshgrid(gridx,gridy)
    
    Uitp = itp.griddata(points,Z, (gridxx, gridyy), method='nearest')
    UitpMask = np.ma.masked_where(maskMesh, Uitp)
    
    return UitpMask

def computeField(U,nodes,elements,N):
    """Compute the field on an NxN grid (with a mask on the mesh)"""
    gridx,gridy,maskMesh = outputGrid(N,nodes)
    Uinterp = interpOnGrid(U,elements,gridx,gridy,maskMesh)
    return gridx,gridy,Uinterp

def exportVTK(U, nodes, filename):
    """Export the nodes values as legacy vtk format"""
    
    out = "# vtk DataFile Version 2.0\nSEM 2D data : nodes values\nASCII\n"
    
    # writing nodes positions
    out += f"DATASET STRUCTURED_GRID\nDIMENSIONS {len(nodes)} {len(nodes)} 1\nPOINTS {len(nodes)} float\n"
    for n in nodes:
        out += f"{np.float32(n.x)} {np.float32(n.y)} {np.float32(0.0)}\n"
        
    # writing displacement on nodes
    out += f"POINT_DATA {len(nodes)}\nSCALARS displacement float 1\nLOOKUP_TABLE default\n"
    for n in nodes:
        out += f"{np.float32(U[n.id])}\n"
      
    with open(filename, "w", encoding="UTF8") as f:
        f.write(out)
        

if __name__ == "__main__":
    
    #---------------------
    
    #file = "../meshes/octogon4.msh"
    # regions = [
    #     (1, 21),
    #     (1, 23),
    #     (1, 22),
    #     (1, 24),
    #     ]
    
    file = "../meshes/square4.msh"
    regions = [
        (1, 31),
        (1, 33),
        (1, 32),
        (1, 34),
        ]
    
    elements, nodes = readGmsh4(file, regions)    
    
    #---------------------
    
    # nodes = [   
    #     Node(0.0, 0.0, 1),
    #     Node(1.0, 0.0, 1),
    #     Node(0.0, 1.0, 1),
    #     Node(1.0, 1.0, 1),
    # ]
    
    # elements = [Element(nodes)]

    #---------------------
 
    plotMesh(elements, nodes, [r[1] for r in regions], allNodes=True)