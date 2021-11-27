#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 18:50:38 2021
@author: Sylvain Brisson
"""

import gmsh
import sys

if len(sys.argv) < 2:
    print("Usage: " + sys.argv[0] + " file.msh")
    exit(0)

gmsh.initialize()
gmsh.open(sys.argv[1])

print("Model name: " + gmsh.model.getCurrent())

# get all elementary entities in the model
entities = gmsh.model.getEntities()

for e in entities:
    print("Entity " + str(e) + " of type " + gmsh.model.getType(e[0], e[1]))
    # get the mesh nodes for each elementary entity
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(e[0], e[1])
    # get the mesh elements for each elementary entity
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[0], e[1])
    # count number of elements
    numElem = sum(len(i) for i in elemTags)
    print(" - mesh has " + str(len(nodeTags)) + " nodes and " + str(numElem) +
          " elements")
    boundary = gmsh.model.getBoundary([e])
    print(" - boundary entities " + str(boundary))
    partitions = gmsh.model.getPartitions(e[0], e[1])
    
# all mesh node coordinates
nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
x = dict(zip(nodeTags, nodeCoords[0::3]))
y = dict(zip(nodeTags, nodeCoords[1::3]))
z = dict(zip(nodeTags, nodeCoords[2::3]))



gmsh.finalize()
