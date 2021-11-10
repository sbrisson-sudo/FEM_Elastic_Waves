## FEM example in python

***GENERAL***
- before write code : construct a **test plan**

***TECHNICAL***
- *boundary-condition* and *right hand side* as general as possible 
- *quadratic elements*

## Introduction and basic implementation for finite element methods

$$
\frac{d}{dx}\left(c(x)\frac{du(x)}{dx}\right) = f(x)
$$
implique
$$
-\int_a^b\frac{d}{dx}\left(c(x)\frac{du(x)}{dx}\right) v(x)dx = \int_a^bf(x)v(x)dx
$$

Avec $u(x)$ trial functions, $v(x)$ the test function -> doivent appartenir à un espace de Sobolev.

Par intégration par partie:

$$
\int_a^b\frac{d}{dx}\left(c(x)\frac{du(x)}{dx}\right) v(x)dx = -c(b)u'(b)v(b) + c(a)u'(a)v(a) + \int_a^bcu'v'dx
$$

## Finite Element Method (MIT)

FEM approximation: passer du système *governing equation* + *boundary condition* à un système algébrique $[K]U = F$.

Pour élasticité, $[K]$ : stiffness matrix

- mesh the domain in elements
- field quantities are interpolated by polynomials over elements
- adjacent elements share the DOF at connecting nodes

**Objectif**
- obtenir équation pour chaque élément
- assembly

## FEM

**Data storage of mesh**
- matrice des coordonnées de chaque noeud
- matrice des indices des noeuds de chaque elements

## Simple implementation matlab

- triangular mesh : piecewise linear basis functions
- quadrilateral mesh : piecewise isoparametric bilinear basis

***FEM program***
1. Preprocessing : definition of data structures : mesh, material properties, BC...
    - **nodal coordinate matrix**
    - **element connectivity matrix** (nodes enumerated in counter clock-wise fashion)
    - **boundary mesh** (1D si problème 2D) : éléments connectant des noeuds 2 à 2
2. Processing : system defined (stiffness matrix, load vector), BC are enforced and system is solved
   - $Kd=f$ : $d$ nodal unknowns
    $$
    u(x) = \sum_{i=1}{nn}\Phi_i(x)d_i
    $$
    Pour $u$ scalaire. Avec $\Phi_i$ les finites elements shape functions, $nn$ le nombre de noeuds.

    Si $u$ est vectoriel, choisir un mapping associant $i$ l'indice du noeud et $j$ la composante vectoriel à une position $n$ dans $d$. (ex : $d=[u1x,u1y,...,unnx,unny]$ ou $d=[u1x,...,unnx,u1y,...,unny]$)
   - **computation of Finite Elements Operators**
    - **linear static** : stiffness matrix et load vector
    - **dynamic** : \+ mass matrix
    1. integrate operators over each elements : choose a **quadrature**
       $$
       \int_{\xi=-1}^{\xi=1}f(\xi)d\xi \simeq \sum_qf(\xi_q)w_q
       $$
       -> connaissance des points d'intégration et des poids d'intégration
    2. assembly -> éviter les boucles for
    - **enforce boundary condition**
    - 
3. Post-processing : retrieve the solution


## Implementation

***INPUTS***
- mesh : liste des position des noeuds et liste des noeuds de chaque elements
- conditions aux limites
- right hand side function (source)

$A$ : stiffness matrix : NxN avec N le nombre de noeuds (DOF)
$b$ : load vector : N

***meshing***
 - triangular 3 nodes
 - four nodes quadrilateral
 - nine nodes quadrilateral

***SRC***
- **reading** read nodes and elem vectors from gmsh file
- **plotting** plot a finite element field, a mesh
- **quadrature** return quadrature rule (points and weights)
- **lagrange_basis** shape functions and gradient of shape functions