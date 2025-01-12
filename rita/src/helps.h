/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2025 Rachid Touzani

    This file is part of rita.

    rita is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rita is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

  ==============================================================================

                               Help messages

  ==============================================================================*/

#pragma once

#include <string>
using std::string;

namespace RITA {

   static const string Data_help = "Command 'data' is useful to list all defined data "
                                   "and entities in a rita execution. It has no arguments.";
   static const string Fct_help = "name, variable, nb, definition";
   static const string Fct_Help = "function [name=<f>] [variable=<x>] [nb=<n>] [definition=<d>]\n\n"
                                  "f: Name to give to function\n"
                                  "x: Name of variable\n"
                                  "n: Number of variables. If many variables and variable name is x,"
                                  " actual variables are x1, x2, ...\n"
                                  "d: Regular expression defining function\n";
   static const string Tab_help = "name, file, min, max, ne, vector";
   static const string Tab_Help = "tabulation name=<n> file=<f> min=<m> max=<M> ne=<e> vector=<v>\n\n"
                                  "n: Name to give to tabulation\n"
                                  "f: tabulation file containing tabulated data\n"
                                  "m: Minimum value of data\n"
                                  "M: Maximum value of data\n"
                                  "e: Number of intervals\n"
                                  "v: Vector containing data to tabulate\n";
   static const string Sample_help = "function, nb, x, y, z, grid, mesh, vector, tabulation";
   static const string Sample_Help = "sample function=<f> [nb=<nx>[,<ny>,<nz>] [x=<mx>,<Mx>] [x=<my>,<My>] [z=<mz>,<Mz>] |"
                                     " grid=<g> | mesh=<ms>] vector=<v>\n\n"
                                     "f:          Name of function to sample\n"
                                     "nx, ny, nz: Number of sampling equidistant points (for 1-D, 2-D or 3-D cases)\n"
                                     "mx, Mx:     Interval of points [mx,Mx]\n"
                                     "my, My:     Interval of points [my,My] (2-D and 3-D cases)\n"
                                     "mz, Mz:     Interval of points [mz,Mz] (3-D case)\n"
                                     "g:          Name of grid used in sampling\n"
                                     "ms:         Name of mesh used in sampling\n"
                                     "v:          Name of sampling vector to create\n";
   static const string Mesh_help = "1d, rectangle, cube, point, curve, surface, volume, contour\n"
                                   "code, generate, nbdof, list, plot, save, read";
   static const string Mesh_Help = "1d        : Data to generate a 1-D mesh\n"
                                   "rectangle : Data to mesh a rectangle\n"
                                   "cube      : Data to mesh a cube (parallelepiped)\n"
                                   "point     : Define a point\n"
                                   "curve     : Define a curve\n"
                                   "surface   : Define a surface\n"
                                   "volume    : Define a volume\n"
                                   "contour   : Define a contour as a sequence of curves or surfaces\n"
//                            "subdomain : Define a subdomain\n"
                                   "code      : Set code for points, lines, surfaces, volumes\n"
                                   "generate  : Generate mesh of a polygon\n"
                                   "plot      : Plot mesh\n"
                                   "clear     : Clear mesh\n"
                                   "read      : Read mesh from file\n";
   static const string Mesh_1d_help = "domain, ne, codes, nbdof, save";
   static const string Mesh_1d_Help = "1d [domain=<m>,<M>] [ne=<n>] [codes=<c1>,<c2>] [nbdof=<d>]\n\n"
                                      "m, M: Extremal coordinates of interval ends. The default values are 0., 1.\n"
                                      "n: Number of elements in the interval. Its default value is 10.\n"
                                      "c1, c2: Codes associated to the first and the last node. These integer values\n"
                                      "        are necessary to enforce boundary conditions. A code 0 (Default value) means\n"
                                      "        no condition to prescribe.\n"
                                      "d: Number of degrees of freedom associated to any generated node. Default value is 1.";
   static const string Mesh_Rect_help = "min, max, ne, codes, nbdof";
   static const string Mesh_Rect_Help = "rectangle [min=<mx>,<my>] [max=<Mx>,<My>] [ne=<nx>,<ny>]  [codes=<c1>,<c2>,<c3>,<c4>]\n"
                                        "          [nbdof=<d>]\n\n"
                                        "mx, my: coordinates of the lower left corner of the rectangle.\n"
                                        "        The default values are 0., 0.\n"
                                        "Mx, My: Coordinates of the upper right corner of the rectangle.\n"
                                        "        The default values are 1., 1.\n"
                                        "nx, ny: Number of elements in x and y direction respectively.\n"
                                        "        Their default value is 10.\n"
                                        "c1, c2, c3, c4: Codes associated to the nodes generated on the lines y=my,\n"
                                        "                x=Mx, y=My, x=mx respectively. These integer values are necessary\n"
                                        "                to enforce boundary conditions. A code 0 (Default value) means no\n"
                                        "                condition to prescribe.\n"
                                        "d: Number of degrees of freedom associated to any generated node. Default value is 1.\n";
   static const string Mesh_cube_help = "min, max, ne, codes, nbdof";
   static const string Mesh_cube_Help = "cube [min=<mx>,<my>,<mz>] [max=<Mx>,<My>,<Mz>] [ne=<nx>,<ny>,<nz>] [codes=<cxm>,<cxM>,<cym>,<cyM>,<czm>,<czMs>]\n"
                                        "     [nbdof=<d>]\n"
                                        "mx, my, mz: Minimal coordinates in each direction.\n"
                                        "            The default values are 0., 0., 0.\n"
                                        "Mx, My, Mz: Maximal coordinates in each direction.\n"
                                        "            The default values are 1., 1.,1.\n"
                                        "nx, ny, nz: Number of elements in x, y and z direction respectively.\n"
                                        "            Their default value is 10.\n"
                                        "cxm, cxM, cym, cyM, czm, czM: Codes associated to the nodes generated on the face x=mx\n" 
                                        "x=Mx, y=my, y=My, z=mz, z=Mz respectively.\n"
                                        "These integer values are necessary to enforce boundary\n"
                                        "conditions. A code 0 (Default value) means no condition to prescribe.\n"
                                        "d: Number of degrees of freedom associated to any generated node. Default value is 1.\n";
   static const string Mesh_code_help = "value, points, curves, surfaces, volumes";
   static const string Mesh_code_Help = "code value=<v> [points=<p1>,<p2>,... | curves=<c1>,<c2>,... | surfaces=<s1>,<s2>,... | volumes=<v1>,<v2>,...]\n\n"
                                        "v: Code value (integer) to assign to an entity\n"
                                        "p1, p2, ...: Points to which the code v is assigned\n"
                                        "c1, c2, ...: Curves to which the code v is assigned\n"
                                        "s1, s2, ...: Surfaces to which the code v is assigned\n"
                                        "v1, v2, ...: Volumes to which the code v is assigned.";
   static const string Mesh_point_help = "label, n, coord, size";
   static const string Mesh_point_Help = "point label=<n> coord=<x>,<y>,<z> size=<h>\n"
                                         "n: Point's label\n"
                                         "x, y, z: Point coordinates. If y and/or z are not given, their value is set to 0?\n"
                                         "This can be the case for 1-D and 2-D.\n"
                                         "h: Mesh size around the point.";
   static const string Mesh_curve_help = "label, line, circle";
   static const string Mesh_curve_Help = "curve label=<n> [line=<n1>,<n2>] [circle=<n1>,<n2>,<n3>]\n\n"
                                         "n: Curve's label\n"
                                         "n1, n2: Labels of points defining a line\n"
                                         "n1, n2, n3: Labels of points defining a circular arc, n1 and n2 are the two \n"
                                         "            points defining the extremities of the arc and n3 is the point defining\n"
                                         "            the center of the circular arc.";
   static const string Mesh_contour_help = "label, curve";
   static const string Mesh_contour_Help = "contour label=<n> curve=<c1>,<c2>,...\n\n"
                                           "n: Contour's label\n"
                                           "c1, c2, ...: Labels of curves defining the contour";
   static const string Mesh_surf_help = "label, contours";
   static const string Mesh_surf_Help = "surface label=<n> contours=<c1>,<c2>,...\n"
                                        "n: Surface's label\n"
                                        "c1, c2, ...: Labels of contours defining the surface.";
   static const string Mesh_subdomain_help = "line, orientation";
   static const string Mesh_subdomain_Help = "line        : Label for a line in subdomain\n"
                                             "orientation : Orientation\n"
                                             "code        : Code to associate to subdomain";
   static const string Mesh_read_help = "mesh, geo, gmsh";
   static const string Mesh_read_Help = "mesh: Read mesh in OFELI mesh file\n"
                                        "geo:  Read mesh in OFELI mesh file\n"
                                        "gmsh: Read mesh in gmsh file";
   static const string Mesh_save_help = "domain, geo, mesh, gmsh, vtk, gnuplot, matlab, tecplot";
   static const string Mesh_save_Help = "domain:  Save domain file\n"
                                        "geo:     Save geo file\n"
                                        "mesh:    Save mesh in OFELI format\n"
                                        "gmsh:    Save mesh in gmsh format\n"
                                        "vtk:     Save mesh in vtk format\n"
                                        "gnuplot: Save mesh in gnuplot format\n"
                                        "matlab:  Save mesh in matlab format\n"
                                        "tecplot: Save mesh in tecplot format";
   static const string Approx_help = "lagrange, fitting, bspline, bezier, nurbs";
   static const string Approx_Help = "lagrange: Polynomial Lagrange interpolation\n"
                                     "fitting:  Least square fitting\n"
                                     "bspline:  Basis spline interpolation\n"
                                     "bezier:   Bézier approximation\n"
                                     "nurbs:    NURBS (Non-uniform rational basis spline) approximation";
   static const string Approx_HHelp = "lagrange: Polynomial Lagrange interpolation\n"
                                      "fitting:  Least square fitting\n"
                                      "bspline:  Basis spline interpolation\n"
                                      "bezier:   Bézier approximation\n"
                                      "nurbs:    NURBS (Non-uniform rational basis spline) approximation";
   static const string Approx_Lagrange_help = "file, tabulation, function";
   static const string Approx_Lagrange_Help = "lagrange tabulation=<t>|file=<f> [function=<fct>]\n\n"
                                              "t:   Tabulation defining data to interpolate\n"
                                              "f:   File where data to interpolate are stored\n"
                                              "fct: Interpolation function name\n";
   static const string Approx_Fitting_help = "file, tabulation, degree, basis, function";
   static const string Approx_Fitting_Help = "fitting tabulation=<t>|file=<f> degree=<d>|basis=<f1,f2,...> [function=<fct>]\n\n"
                                             "t:           Tabulation defining data to interpolate\n"
                                             "f:           File where data to interpolate are stored\n"
                                             "d:           Degree of polynomial fitting (must be >0,default: 1)\n"
                                             "b1, b2, ...: List of basis functions"
                                             "fct:         Fitting function name\n";
   static const string Approx_BSpline_help = "file, vertices, order, ncu, ncv, npu, npv, vector";
   static const string Approx_BSpline_Help = "bspline file=<f>|vertices=<v> [order=<o>] [nc=<ncu>,<ncv>] np=<npu>[,<npv>] [vector=<ve>]\n\n"
                                             "f:   File containing control polygon vertex data: 3 coordinates for each point\n"
                                             "v:   List of vertices of control polygon: 3 coordinates for each point\n"
                                             "o:   Order of the B-spline basis function (Default: 2)\n"
                                             "ncu: Number of polygon vertices in the u-direction. Mandatory if vector v is given\n"
                                             "ncv: Number of polygon vertices in the v-direction if a bspline surface is to be defined\n"
                                             "npu: Number of points to define the resulting curve, or number of points in the u-direction for the resulting surface\n"
                                             "npv: Number of points in the v-direction for a resulting surface\n"
                                             "ve:  Name of vector defining the resulting curve or surface: : 3 coordinates for each point\n";
   static const string Approx_Bezier_help = "file, vertices, ncu, ncv, npu, npv, vector";
   static const string Approx_Bezier_Help = "bezier file=<f>|vertices=<v> [nc=<ncu>,<ncv>] np=<npu>[,<npv>] [vector=<ve>]\n\n"
                                            "f:   File where data to approximate are stored\n"
                                            "v:   List of vertices of control polygon: 3 coordinates for each point\n"
                                            "ncu: Number of polygon vertices in the u-direction. Mandatory if vector v is given\n"
                                            "ncv: Number of polygon vertices in the v-direction if a bspline surface is to be defined\n"
                                            "npu: Number of points to define the resulting curve, or number of points in the u-direction for the resulting surface\n"
                                            "npv: Number of points in the v-direction for a resulting surface\n"
                                            "ve:  Name of vector defining the resulting curve or surface: : 3 coordinates for each point\n";
   static const string Approx_Nurbs_help = "file, tabulation, degree, function";
   static const string Approx_Nurbs_Help = "nurbs file=<f>|vertices=<v> [order=<o>] [nc=<ncu>,<ncv>] np=<npu>[,<npv>] [vector=<ve>]\n\n"
                                            "f:   File where data to approximate are stored\n"
                                            "o:   Order of the Nurbs approximation (Default: 2)\n"
                                            "v:   List of vertices of control polygon: 3 coordinates for each point\n"
                                            "ncu: Number of polygon vertices in the u-direction. Mandatory if vector v is given\n"
                                            "ncv: Number of polygon vertices in the v-direction if a bspline surface is to be defined\n"
                                            "npu: Number of points to define the resulting curve, or number of points in the u-direction for the resulting surface\n"
                                            "npv: Number of points in the v-direction for a resulting surface\n"
                                            "ve:  Name of vector defining the resulting curve or surface: : 3 coordinates for each point\n";
   static const string Int_help = "function, definition, x, variable, ne, formula, result";
   static const string Int_Help = "integration [function=<f>|definition=<exp>] [x=<xmin>,<xmax>] [result=<r>]\n"
                                  "            [variable=<x>] [ne=<nx>] [formula=<m>]\n\n"
                                  "f:          Name of already defined function to integrate\n"
                                  "exp:        Expression defining function to integrate\n"
                                  "xmin, xmax: Integration is made on the interval (min,max)\n"
                                  "x:          name of variable\n"
                                  "ne:         Number of subdivisions of the interval\n"
                                  "m:          Numerical integration formula. To choose among the values: left-rectangle,\n"
                                  "            right-rectangle, mid-point, trapezoidal, simpson, gauss-legendre, gauss-lobatto\n"
	                                "            For the Gauss formulae, the number of points can be specified by typing \n"
 	                                "            gauss-legendre, 2 for instance\n"
                                  "r:          Parameter where result is stored\n";
   static const string LS_help = "matrix, rhs, solution, solver, preconditioner, initial";
   static const string LS_Help = "ls matrix=<A> rhs=<b> [solution=<x>] [solver=<s>] [preconditioner=<p>] [initial=<i>] [name=<n>]\n\n"
                                 "A: Matrix of the linear system (already defined)\n"
                                 "b: Vector containing right-hand side of the system\n"
                                 "x: Vector that contains solution (after execution). By default, the solution is stored in the\n"
                                 "   right-hand side vector\n"
                                 "s: Linear solver, can be direct or iterative. By default the direct solver is used\n"
                                 "p: Preconditioner for the iterative solution case. By default the diagonal preconditioner is used.\n"
                                 "i: Initial guess for iterations. By default the 0-vector is used.\n\n";
   static const string AE_help = "function, definition, init, variable, nls\n";
   static const string AE_hhelp = "size, function, definition, jacobian, init, variable, nls\n";
   static const string AE_Help = "algebraic [size=<n>] [function=<f>] [definition=<d>] [jacobian=<j>] variable=<x> [init=<i>] [nls=<s>]\n\n"
                                 "n: Size of algebraic system (Number of equations)\n"
                                 "f: Name of a defined function\n"
                                 "d: Give function expression F to define the algebraic equation F(x)=0\n"
                                 "j: Define jacobian of mapping (for the Newton's algorithm)\n"
                                 "x: Variable (vector) name as unknown of the equation\n"
                                 "i: Initial guess for iterations\n"
                                 "s: Nonlinear equation iteration solver\n";
   static const string AE_HHelp = "size:       Size of algebraic system (Number of equations)\n"
                                  "function:   Name of already defined function\n"
                                  "definition: Give function expression F to define the algebraic equation F(x)=0\n"
                                  "jacobian:   Define jacobian of mapping (for the Newton's algorithm)\n"
                                  "variable:   Variable (vector) name as unknown of the equation\n"
                                  "init:       Initial guess for iterations\n"
                                  "nls:        Nonlinear equation iteration solver\n";
   static const string ODE_help = "size, function, definition, variable, initial, final-time, time-step\n"
                                  "scheme, phase, analytic, analytic-function, err";
   static const string ODE_Help = "ode [size=<n>] [function=<f>] [definition=<d>] [variable=<v>] [initial=<i>] [final-time=<ft>]\n"
                                  "   [time-step=<ts>] [scheme=<s>] [analytic=<a>] [analytic-function=<af>] [error=<e>]\n\n"
                                 "n:  Size of differential system: (Number of equations)\n"
                                 "f:  Function defining ode (already defined)\n"
                                 "d:  Expression defining ode\n"
                                 "v:  Variable name as unknown of the equation\n"
                                 "i:  Initial condition\n"
                                 "ft: Final time\n"
                                 "ts: Time step\n"
                                 "s:  Time integration scheme\n"
                                 "a:  Expression of analytical solution of the ODE as function of x and t\n"
                                 "af: Function defining the analytical solution of the ODE\n"
                                 "e:  Name of parameter containing error\n";
   static const string ODE_HHelp = "size:              Size of differential system: (Number of equations)\n"
                                   "function:          Function defining ode (already defined)\n"
                                   "definition:        Expression defining ode\n"
                                   "variable:          Variable (or vector) name as unknown of the equation\n"
                                   "initial:           Initial condition\n"
                                   "final-time:        Final time\n"
                                   "time-step:         Time step\n"
                                   "scheme:            Time integration scheme\n"
                                   "phase:             Compute phase for a given component and store it in a history vector file\n"
                                   "analytic:          Expression of analytical solution of the ODE as function of x and t\n"
                                   "analytic-function: Function defining the analytical solution of the ODE\n"
                                   "err:               Name of parameter containing error\n";
   static const string PDE_help = "variable, coef, axi, init, bc, bf, source, sf, traction, space, la, nls, final-time\n" 
                                  "time-step, scheme, name, analytic, analytic-function, err, save-every";
   static const string PDE_Help = "vector:            Vector name of an unknown of the equation\n"
                                  "coef:              PDE coefficients\n"
                                  "axi:               Choose axisymmetric geometry\n"
                                  "init:              Set initial condition or guess for pde (if transient analysis)\n"
                                  "bc:                Set boundary conditions\n"
                                  "source:            Set sources or body forces\n"
                                  "sf:                Set side (boundary) forces\n"
                                  "space:             Space discretization method\n"
                                  "ls:                Set linear system solver\n"
                                  "nls:               Set nonlinear system iteration procedure\n"
                                  "final-time:        Final time (if transient analysis)\n"
                                  "time-step:         Time step (if transient analysis)\n"
                                  "scheme:            Time integration scheme (if transient analysis)\n"
                                  "name:              Set name of the defined PDE\n"
                                  "analytic:          Define expression of analytic solution\n"
                                  "analytic-function: Define name of already defined function as analytic solution\n"
                                  "err:               Names of 2 variables containing L2 and Inifinity norms errors\n"
                                  "save-every:        Save PDE solution every n time steps\n\n"
                                  "Available pde's are so far:\n"
                                  "linear:                       Generic linear partial differential equations of first and second order in time and space\n"
                                  "laplace:                      Laplace (or Poisson) equation\n"
                                  "heat:                         The linear heat equation\n"
                                  "wave:                         The linear wave equation\n"
                                  "transport:                    The linear 1-D transport (convection) equation\n"
                                  "linear-elasticity:            Linearized elasticity in 2-D and 3-D\n"
                                  "truss:                        Plane truss equation in structural mechanics\n"
                                  "beam:                         3-D beam equation in structural mechanics\n"
                                  "incompressible-navier-stokes: Time-dependent incompressible Navier-Stokes equations\n";
   static const string PDE_Coef_help = "c00, c10, c01, c20, c02, rho, density, Cp, specific-heat, kappa, thermal-conductivity\n"
                                       "Mu, magnetic-permeability, sigma, electric-conductivity, mu, viscosity, epsilon\n"
                                       "electric-permittivity, omega, angular-frequency, beta, thermal-dilatation, v,\n"
                                       "velocity, young, poisson";
   static const string PDE_Coef_Help = "The arguments of this command must be in function of the chosen PDE.\n\n"
                                       "For a generic linear pde (of the form: a00*u + a10*du/dt + a20*d2u/dt2 + a01*nabla(u) - div(a02 nabla(u)):\n"
                                       "[c00=<a00>] [c10=<a10>] [c01=<a01>] [c20=<a20>] |c02=<a02>]\n"
                                       "Every value a00, a10, ... can be a constant or a string defining the expression of a function (with variables x,y,z,t)\n"
                                       "When a vector is to be given, its components are to be given component by component, separated by a comma\n\n"
                                       "For diffusion-convection problems:\n"
                                       "[rho|density=<r>] [Cp|specific-heat=c] [kappa|thermal-conductivity=<k>] [v|velocity=vv]\n"
                                       "r:  Density\n"
                                       "c:  Specific heat\n"
                                       "k:  Thermal conductivity\n"
                                       "vv: velocity\n\n"
                                       "For fluid dynamics problems:\n"
                                       "[rho|density=<r>] [mu|viscosity=<m>] [beta|thermal-dilatation=<b>]\n"
                                       "r: Density\n"
                                       "m: Dynamic viscosity\n"
                                       "b: Thermal dilatation\n\n"
                                       "For solid mechanics problems:\n"
                                       "[young=<y>] [poisson=<p>]\n"
                                       "y: Young modulus\n"
                                       "p: Poisson ratio\n\n"
                                       "For electromagnetic problems:\n"
                                       "[rho|density=<r>] [Mu|magnetic-permeability=m] [sigma|electric-conductivity=<s>] [epsilon|electric-permittivity=e]\n"
                                       "[omega-angular-frequency=o]\n"
                                       "m: Magnetic permeability\n"
                                       "s: Electric conductivity\n"
                                       "e: Electric permittivity\n"
                                       "o: Angular frequency\n";
   static const string Opt_hhelp = "size, function, objective, lp, gradient, hessian, low-bound, up-bound\n"
	                                 "ineq-constraint, eq-constraint, penalty, variable, init, algorithm";
   static const string Opt_Help = "optimization [size=<n>] [function=<f>] [objective=<c>] [lp] [gradient=<g>] [hessian=<h>]\n"
                                  "[low-bound=<lb>] [up-bound=<ub>] [variable=<v>] [initial=<i>] [algorithm=<a>]\n\n"
                                  "n:  Size of optimization problem (Number of optimization variables)\n"
                                  "f:  Objective (cost) function (already defined)\n"
                                  "c:  Expression of objective (cost) function\n"
                                  "lp: Linear programming is set\n"
                                  "g:  Gradient of objective function\n"
                                  "h:  Hessian of objective function\n"
                                  "lb: Define a lower bound for a given variable as constraint\n"
                                  "ub: Define an upper bound for a given variable as constraint\n"
                                  "v:  Variable name as unknown of the optimization proble\n"
                                  "i:  Initial guess for iterations\n"
                                  "a:  Optimization algorithm\n";
   static const string Opt_HHelp = "size:          Size of optimization problem (Number of optimization variables)\n"
                                   "function:      Objective (cost) function (already defined)\n"
                                   "objective:     Expression of objective (cost) function\n"
                                   "lp:            Linear programming is set\n"
                                   "gradient:      Gradient of objective function\n"
                                   "hessian:       Hessian of objective function\n"
                                   "low-bound:     Define a lower bound for a given variable as constraint\n"
                                   "up-bound:      Define an upper bound for a given variable as constraint\n"
                                   "ge-constraint: Define a (>=) inequality constraint (for Linear Programming problems only)\n"
                                   "le-constraint: Define a (<=) inequality constraint\n"
                                   "eq-constraint: Define an equality constraint\n"
                                   "penalty:       Penalty parameter (small) to enforce constraints\n"
                                   "variable:      Variable name as unknown of the optimization problem\n"
                                   "init:          Initial guess for iterations\n"
                                   "algorithm:     Set optimization algorithm\n";
   static const string Eig_help = "matrix, symmetric, method, nb, eigv, evec";
   static const string Eig_Help = "eigen matrix=<M> [symmetric] [method=<m>] [nb=<n>] [eigv] [evect=<p>]\n\n"
                                  "M:         Name of matrix for which eigenvalues are to be extracted\n"
                                  "symmetric: If matrix is symmetric\n"
                                  "m:         Method to extract eigenvalues: subspace or qr\n"
                                  "n:         Number of eigenvalues. By default, all eigenvalues are extracted\n"
                                  "eigv:      Switch to compute eigenvectors as well\n"
                                  "p:         Prefix for names of files containing eigenvectors\n";
   static const string Solve_help = "problem, analysis";
   static const string Solve_Help = "solve [problem=<p1>,<p2>,...] [analysis=<a>]\n"
                                    "p1, p2: List of coupled problems to solve\n"
                                    "a:      stationary or transient (All other types are automatically determined)\n";
   static const string Solve_HHelp = "problem:  Set equation (or system of equations) to solve\n"
                                     "analysis: Define analysis type: stationary or transient (All other types are automatically determined)\n";
   static const string Plot_help = "function, file, tabulation, vector, history, mesh, x, y, isolines, contour, use";
   static const string Plot_Help = "plot [function=<f> | file=<fi> | tabulation=<t> | vector=<v> | history=<h> | mesh=<ms>]\n"
                                   "     [x=<mx>,<Mx>] [y=<my>,<My>] [log=<lx>,<ly>] [title=<ti>] [label=<l>] [use=<u>] [component=<c>]\n\n"
                                   "f:        Name of function to plot\n"
                                   "fi:       File with contents to plot\n"
                                   "v:        Vector to plot\n"
                                   "h:        History vector to plot\n"
                                   "t:        Tabulation to plot\n"
                                   "ms:       Mesh to plot\n"
//                            "i:        Plot vector isolines\n"
//                            "ct:       Plot vector contour\n"
                                   "mx,Mx:    Minimal and maximal values of first variable\n"
                                   "my,My:    Minimal and maximal values of second variable\n"
                                   "lx,ly:    1 (log scale) or 0 (linear scale) for x and/or y variables\n"
                                   "ti:       Plot title\n"
                                   "l:        Label of curve\n"
                                   "u:        Software for plotting. Must be gnuplot (Default) or gmsh if available\n"
                                   "c:        Component to plot\n";
   static const string Plot_HHelp = "function:  Name of function to plot\n"
                                    "file:      File with contents to plot\n"
                                    "vector:    Vector to plot\n"
                                    "history:   History vector to plot\n"
                                    "tab:       Tabulation to plot\n"
                                    "mesh:      Mesh to plot\n"
                                    "isolines:  \n"
                                    "contour:   \n"
                                    "x:         Minimal and maximal values of first variable\n"
                                    "y:         Minimal and maximal values of second variable\n"
                                    "log:       1 (log scale) or 0 (linear scale) for x and/or y variables\n"
                                    "mark:      Marks on pkot points\n"
                                    "title:     Plot title\n"
                                    "use:       Software for plotting. Must be gnuplot (Default) or gmsh if available\n"
                                    "component: Component of a field to plot\n";
   static const string Config_help = "verbosity, save-results, history, log, echo";
   static const string Config_Help = "set [verbosity=<v>] [history=<h>] [log=<l>] [echo on|off]\n\n"
                                     "v:    Verbosity parameter (0, 1 or 2)\n"
                                     "h:    Name of history file\n"
                                     "l:    Name of log file\n"
                                     "echo: Choose if entered commands are echoed or not\n";
   static const string Grid_help = "name, x, y, z, ne";
   static const string Grid_Help = "grid name=n x=mx,Mx y=my,My z=mz=Mz [ne=nx,ny,nz]\n\n"
                                   "n: Name of the grid entity\n"
                                   "mx, my, mz: Minimum value of coordinate in each direction\n"
                                   "            If only one value is given the grid is 1D, if 2 values\n"
                                   "            are given the grid is 2D, \n"
                                   "Mx, My, Mz: Maximum value of coordinate in each direction\n"
                                   "nx, ny, nz: Number of sub-intervals in each direction. ";
   static const string Vect_help = "name, size, grid, mesh, file, define, set";
   static const string Vect_Help = "vector [name=nm] [size=n] [grid=g] [mesh=m] [file=f] [define=d] [set=s]\n"
                                   "nm: Name of vector\n"
                                   "n: Size of vector. Vector is then initialized to zero\n"
                                   "g: Name of grid to associate to vector\n"
                                   "m: Name of mesh to associate to vector\n"
                                   "f: File name, if vector is read from a file\n";
//                                   "d: \n"
//                                   "s: \n";
   static const string HVect_Help = "history <v> <V>\n"
                                    "v: Name of (existing) vector to add to a history vector\n"
                                    "V: Name of a hisory vector. If this one does not exist, it will be created by this command";
   static const string Matrix_help = "name, file, storage, nr, nc, define, set";
   static const string Matrix_Help = "matrix [name=<m>] [file=<f>] [storage=<s>] [nr=<n1>] [nc=<n2>] [define=<d>] [set=<st>]\n\n"
                                     "m: Name to give to matrix. If this one is read in file, matrix name can be retrieved from the file\n"
                                     "f: File (OFELI format) containing matrix data\n"
                                     "s: Storage type\n"
                                     "n1: Number of rows\n"
                                     "n2: Number of columns\n";
//                                     "d: \n"
//                                     "st: \n";
   static const string Save_help = "name, file, format, every";
   static const string Save_Help = "save name=<n> file=<fi> [format=<fo>] [every=<e>]\n"
                                   "n: Name of entity to save\n"
                                   "fi: File where to save data\n"
                                   "f: File format\n"
                                   "e: Frequency of saving (Default: 1)\n";
    static const string Calc_Help = "The application rita works also as a simple calculator that provides basic operations on\n"
                                    "scalars, vectors and matrices.\n"
                                    "- Scalars can be defined simply by their initialization or by using standard functions.\n"
                                    "  Examples:\n"
                                    "  d = 12\n"
                                    "  x = d*log(d)\n"
                                    "  z = exp(x)-1\n"
                                    "- Vectors can be initialized as in the examples:\n"
                                    "  v = {x,3,log(7)} # This is a row vector\n"
                                    "  w = vector(3)    # This is a column vector\n"
                                    "  u = v'+w\n"
                                    "- Matrices can be defined as in the following examples:\n"
                                    "  M = matrix(3,3)\n"
                                    "  M[0,1]=-1\n"
                                    "  w = M*v\n"
                                    "  Note that matrix entries are denoted by M[i,j] where the indices i,j start at 0.\n"
                                    "- All the computed entities (scalars, vectors, matrices) are stored in rita and can be\n"
                                    "  used for any other purposes.\n";
    static const string Print_Help = "Command '=' (or 'print') prints the content of a given entity defined by its name. For instance:\n"
                                     "'= PDE-1' (or print 'PDE-1') displays information about the entity named 'PDE-1'\n"
                                     "which can be a partial differential equation.\n";
    static const string Help_Topics = "help, license, set, load, end, exit, quit, calculator, vector, matrix, tabulation,\n"
                                      "mesh, grid, function, sample, plot, history, data, list, approximation, integration, ls,\n"
                                      "algebraic (or ae), ode, pde, optimization, eigen, solve, rename, remove, calc, print\n";

} /* namespace RITA */
