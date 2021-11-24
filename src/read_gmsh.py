import matplotlib.pyplot as plt

import numpy as np

def readGmsh(meshFile):
    """
    Read a gmsh file (version 2) of 2D quadrangles.
    Return an array of nodes positions and an array of the nodes composing each element.
    """

    nodes = []
    elements = []

    with open(meshFile, "r") as f:

        lines = f.readlines()
        nodesStart, nodesEnd, elemStart, elemEnd = 0,0,0,0

        for (i,l) in enumerate(lines):
            if l == "$Nodes\n": nodesStart = i+2
            elif l == "$EndNodes\n": nodesEnd = i
            elif l == "$Elements\n": elemStart = i+2
            elif l == "$EndElements\n": elemEnd = i

        print(nodesStart, nodesEnd, elemStart, elemEnd)

        nodes = [list(map(float, l.split()[1:3])) for l in lines[nodesStart:nodesEnd]]
        for l in lines[elemStart:elemEnd]:
            if len(l.split())==9:
                elements.append(list(map(int, l.split()[5:])))

    return np.array(nodes), np.array(elements)


def plotMesh(nodes, elements):

    fig, ax = plt.subplots()

    for e in elements:
        n1,n2,n3,n4 = e
        ne = [nodes[n1-1],nodes[n2-1],nodes[n3-1],nodes[n4-1]]
        plt.fill([n[0] for n in ne], [n[1] for n in ne], color="greenyellow", alpha=0.5, ec="k")

    ax.set_aspect("equal", "box")

    plt.show()


if __name__ == "__main__":

    nodes, elements = readGmsh("../gmsh_meshes/octogon.msh")
    print(f"{len(nodes)} nodes")
    print(f"{len(elements)} elements")
    plotMesh(nodes, elements)