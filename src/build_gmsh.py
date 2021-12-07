import gmsh

gmsh.initialize()
gmsh.model.add("t3")

### To modify global mesh size

gmsh.option.setNumber("Mesh.MeshSizeFactor", 0.18)
#gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 2)

### To modify colors
gmsh.option.setColor("Geometry.Points", 255, 0, 0)
gmsh.option.setColor("Mesh.Points", 255, 127, 0)


### To build the figure
geo = gmsh.model.geo

#Building
#h = [0, 0.2, 3, 3.2, 6, 6.2, 9, 9.2, 12, 12.2, 15, 15.2] #height
#n = len(h)
#w = [0, 0.2, 6, 6.2, 10, 10.2, 12, 12.2, 16, 16.2, 22, 22.2] #width
#m = len(w)

#h = [0, 0.2, 3, 3.2, 6, 6.2, 9, 9.2] #height
#n = len(h)
#w = [0, 0.2, 6, 6.2, 10, 10.2] #width
#m = len(w)

h = [0, 0.2, 3, 3.2, 6, 6.2] #height
n = len(h)
w = [0, 0.2, 6, 6.2] #width
m = len(w)

#Points-Lines-Surface for exterior
geo.addPoint(w[0], h[0], 0, meshSize=0.0, tag=1)
geo.addPoint(w[-1], h[0], 0, tag=2)
geo.addPoint(w[0], h[-1], 0, tag=3)
geo.addPoint(w[-1], h[-1], 0, tag=4)
geo.addLine(1, 2, 1)
geo.addLine(2, 4, 2)
geo.addLine(4, 3, 4)
geo.addLine(3, 1, 3)
geo.addCurveLoop([1, 2, 4, 3], 1)

# To set boundaries
gmsh.model.addPhysicalGroup(1, [1], 1)

#Points for holes
for i in range(n-2):
    for j in range(m-2):
        geo.addPoint(w[j+1], h[i+1], 0, meshSize=0.0, tag= 5 + (m-2)*i + j)

#Lines-Surface for holes
l = [1]
l_completed = [2,3,4]
for i in range(n//2 - 1): #nb of holes in height dim
    for j in range(m//2 -1): #nb of holes in width dim
        a = 5 +2*j +2*(m-2)*i
        b = 6 +2*j +2*(m-2)*i
        c = 6+(m-2) +2*j +2*(m-2)*i
        d = 5+(m-2) +2*j +2*(m-2)*i
        geo.addLine(a, b, a)
        geo.addLine(b, c, b)
        geo.addLine(c, d, c)
        geo.addLine(d, a, d)
        geo.addCurveLoop([a, b, c, d], a)
        l.append(a)
        l_completed.extend([a, b, c, d])

# To set boundaries
gmsh.model.addPhysicalGroup(1, l_completed, 2)

#Complete the build
geo.addPlaneSurface(l, 1)
gmsh.model.geo.synchronize()

### To generate the mesh
gmsh.model.mesh.generate(2)

### To recombine
gmsh.model.mesh.recombine()
#gmsh.model.mesh.optimize("HighOrderElastic")

### To visualize the mesh
gmsh.fltk.run()

### To save the mesh
gmsh.option.setNumber("Mesh.SaveAll", 1)
#gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("../FEM_Elastic_Waves/meshes/t3.msh")

gmsh.finalize()