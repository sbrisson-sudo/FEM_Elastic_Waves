## Geometry commands
- ***POINTS***
  - `Point(tag) = {x, y, z};`
  - `Physical Point {"name", tag} = {point-tag};`
- ***CURVES***
  - `Line(tag) = {pt1, pt2}`
  - `Curve Loop(tag) = {curve_list_tags};`
  Autre : Bézier, Spline, Circle, Ellipse...
  - `Physical Curve("name",tag) = {curve_list_tags}`

- ***SURFACES***
  - `Plane Surface(tag) = {curve_loop_tags};`
  - `Physical Surface("name", tag) = {surfaces_list};`


## FORMAT

```

$Entities
8 8 1 0 \\ nombre de points, de lignes, de surfaces
...points...
11 0 0 0 10 0 0 1 21 2 1 -2 \\ courbe : tag,minx,miny,minz,max...,maxz,num physical tags, physical tags
$Nodes
17 107 1 107 \\ 17 entities, 107 noeuds indéxés de 1 à 107
0 1 0 1     \\ entité de dimension 0, tag 1, garbage, avec 1 élément
1           \\ Noeud n° 1
0 0 0       \\ Coordonnées
...
1 11 0 3      \\ courbe, dimension = 1
9
10
11
2.499999999995363 0 0
4.999999999992402 0 0
7.4999999999962 0 0
$EndNodes
$Elements
1 90 1 90     \\ nb entité , nombre éléments
2 111 3 90    \\ dimension entité, tag entité, type élément, nombre d'éléments
1 75 89 64 78  \\ élément 1


