import gmsh

gmsh.initialize()
gmsh.model.add("building")

# gmsh.option.setNumber("General.Verbosity", 0)

# element order
N = 2
gmsh.option.setNumber("Mesh.ElementOrder", N)

# recombine into quadrangles by default
gmsh.option.setNumber("Mesh.RecombineAll", 1)

# recombine into quadrangles by default
gmsh.option.setNumber("Geometry.AutoCoherence", 1)

# assign geometry to variable
geo = gmsh.model.geo

#Building
heights = [0, 0.5, 3, 3.5, 6, 6.5] #height
H = len(heights)
widths = [0, 0.5, 5, 5.5, 10, 10.5] #width
W = len(widths)

dx = 0.25 # taille éléments

points = [[geo.addPoint(w,h,0,dx) for w in widths] for h in heights]

surfaces = []

for i in range(0,H,2):
    pts = [points[i][0], points[i][W-1], points[i+1][W-1], points[i+1][0]]
    lines = [geo.addLine(pts[i], pts[(i+1)%4]) for i in range(4)]
    curveLine = geo.addCurveLoop(lines)
    surfaces.append(geo.addPlaneSurface([curveLine]))

for i in range(1,H-2,2):
    for j in range(0,W,2):
        pts = [points[i][j], points[i][j+1], points[i+1][j+1], points[i+1][j]]
        lines = [geo.addLine(pts[i], pts[(i+1)%4]) for i in range(4)]
        curveLine = geo.addCurveLoop(lines)
        surfaces.append(geo.addPlaneSurface([curveLine]))


# add physical groups at boundaries
cornerPoints = [
    points[0][0],
    points[0][W-1],
    points[H-1][W-1],
    points[H-1][0]
]

lines = [geo.addLine(cornerPoints[i], cornerPoints[(i+1)%4]) for i in range(4)]
bot,right,top,left = [2,3,4,5]

gmsh.model.addPhysicalGroup(1, [lines[0]], bot)
gmsh.model.addPhysicalGroup(1, [lines[1]], right)
gmsh.model.addPhysicalGroup(1, [lines[2]], top)
gmsh.model.addPhysicalGroup(1, [lines[3]], left)

for surface in surfaces:
    geo.mesh.setTransfiniteSurface(surface)
    




geo.synchronize()

# # To generate the mesh
# gmsh.model.mesh.generate(2)

### To visualize the mesh
gmsh.fltk.run()

# # To save the mesh
gmsh.option.setNumber("Mesh.SaveAll", 1)
# gmsh.write(f"../FEM_Elastic_Waves/meshes/building{N}.msh")

gmsh.finalize()