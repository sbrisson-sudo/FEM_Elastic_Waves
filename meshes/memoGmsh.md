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