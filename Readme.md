# Elastic wave simulations using SEM

- `\src`
  - `sem1d.ipynb` wave equation solver (1D)
  - `fem2D.py` main module, used by all 2D FEM and SEM notebooks
  - `fem2d.ipynb` Poisson equation solver, using P1 elements (linear)
  - `fem2dDyn.ipynb` Accoustic wave equation solver using P1 elements (not working)
  - `fem2dElast.ipynb` Elastic static equation solver using P1 elements (not working)
  - `sem2d.ipynb` Poisson equation solver (SEM)
  - `sem2dDyn.ipynb` Accoustic wave equation solver (SEM)
- `\meshes` gmsh mesh files `.msh` and geometry files (for building meshes) `.geo`
- `\quadratureRules` data of GL and GLL quadrature rules
- `\references` documentation on FEM,SEM an elastic laws