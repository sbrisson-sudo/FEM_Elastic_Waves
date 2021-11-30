import gmsh

gmsh.initialize()
gmsh.model.add("t1")

### To modify global mesh size

gmsh.option.setNumber("Mesh.MeshSizeFactor", 0.1)

### To build the figure
geo = gmsh.model.geo

#Building
h = [0, 0.2, 3, 3.2, 6, 6.2, 9, 9.2, 12, 12.2, 15, 15.2] #height
n = len(h)
w = [0, 0.2, 6, 6.2, 10, 10.2, 12, 12.2, 16, 16.2, 22, 22.2] #width
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

#Points for holes
for i in range(n-2):
    for j in range(m-2):
        geo.addPoint(w[j+1], h[i+1], 0, meshSize=0.0, tag= 5 + (m-2)*i + j)

#Lines-Surface for holes
l = [1]
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

#Complete the build
geo.addPlaneSurface(l, 1)
gmsh.model.geo.synchronize()

### To generate the mesh
gmsh.model.mesh.generate(2)

### To recombine
gmsh.model.mesh.recombine()

### To visualize the mesh
gmsh.fltk.run()

### To save the mesh
gmsh.option.setNumber("Mesh.SaveAll", 1)
#gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("../meshes/t1_quad.msh")

gmsh.finalize()