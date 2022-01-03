#!/bin/sh

# ==============================================================================
#
#                                    O  F  E  L  I
#
#                           Object  Finite  Element  Library
#
# ==============================================================================
#
#   Copyright (C) 1998 - 2022 Rachid Touzani
#
#   This file is part of OFELI.
#
#   OFELI is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   OFELI is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with OFELI. If not, see <http://www.gnu.org/licenses/>.
#
# ==============================================================================
#
#                 A script to test various ofeli applications
#                           (Demos and utilities)
#
# ==============================================================================

echo "Testing OFELI Demos and Utilities ..."

echo "================================================================="
echo "Testing OFELI Demos ..."
echo "================================================================="
cd demos/1D
echo
echo "-----------------------------------------------------------------"
echo "Test 1-D demo programs"
echo "-----------------------------------------------------------------"
echo "Test elliptic equation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./elliptic
fi
echo "Test heat equation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./heat 10 0.1
fi
echo "Test transport equation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./transport 10 0.1
fi

cd ../laplace
echo "-----------------------------------------------------------------"
echo "Test Laplace equation programs"
echo "-----------------------------------------------------------------"
echo "Test the 2-D Laplace equation using P1 finite elements (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./laplace_demo1 proj2.dat
fi
echo "Test the 2-D Laplace equation using P2 finite elements (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./laplace_demo2 proj2.dat
fi
echo "Test the 3-D Laplace equation using P1 finite elements (y/n) ? \c"
if test "$ans" = "y" ; then
   ./laplace_demo3 proj3.dat
fi
echo "Test the 2-D Steklov-Poincaré problem using P0 boundary elements (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./laplace_demo4 proj.dat
fi

echo "-----------------------------------------------------------------"
echo "Testing Heat Transfer Demos ..."
echo "-----------------------------------------------------------------"
cd ../thermal/stdc2
echo "Test demo for 2-D steady state thermal computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./stdc2 proj.dat
fi
cd ../ttd2
echo "Test demo for 2-D transient thermal computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ttd2 proj.dat
fi
cd ../std3
echo "Test demo for 3-D steady state thermal computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./std3 beam.dat
fi

echo "-----------------------------------------------------------------"
echo "Testing Solid Mechanics Demos ..."
echo "-----------------------------------------------------------------"

cd ../../solid/truss
echo "Test demo for 2-D truss structure (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./truss truss.dat
fi

cd ../beam
echo "Test demo for 3-D beams (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./beam beam.dat
fi

cd ../lelas2d
echo "Test demo for 2-D linear elasticity computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./lelas2d beam.dat
fi

cd ../contact
echo "Test demo for 2-D linear elasticity with contact (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./contact beam.dat
fi

cd ../lelas3d
echo "Test demo for 3-D linear elasticity computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./lelas3d beam.dat
fi

echo "-----------------------------------------------------------------"
echo "Testing Fluid Dynamics Demos ..."
echo "-----------------------------------------------------------------"
cd ../../fluid/tins2
echo "Test demo for 2-D transient incompressible fluid flow computation (Projection method) (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./tins2 cavity.dat
fi

cd ../tiff2
echo "Test demo for 2-D transient incompressible fluid flow computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./tiff2 cavity.dat
fi

cd ../lh2d
echo "Test demo for linear hyperbolic equation in 2-D (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./lh2d proj.dat
fi

cd ../euler-2d
echo "Test demo for 2-D Compressible Euler equations (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./euler2d proj.dat
fi

echo "-----------------------------------------------------------------"
echo "Testing Electromagnetic Demos ..."
echo "-----------------------------------------------------------------"
cd ../../electromagnetics/Helmholtz
echo "Test demo for Helmholtz equation in bounded media (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./helmholtz proj.dat
fi

echo "-----------------------------------------------------------------"
echo "Testing Linear System solver Demos ..."
echo "-----------------------------------------------------------------"
cd ../../solvers/LS
echo "Test demo for a direct solver (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ls_demo1 10
fi

echo "Test demo for an iterative solver (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ls_demo2 mesh.m
fi

echo "-----------------------------------------------------------------"
echo "Testing ODE solver Demos ..."
echo "-----------------------------------------------------------------"
cd ../../solvers/ODE
echo "Test demo for 1-st order ode given by expression (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ode_demo1 0.1
fi

echo "Test demo for 1-st order ode given by data (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ode_demo2 0.1
fi

echo "Test demo for a system of 1-st order ode (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ode_demo3 0.1
fi

echo "Test demo for a system of 2-nd order linear ode (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ode_demo4 0.1
fi

echo "Test demo for a parabolic pde (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
    ./ts_demo1 mesh1.m
fi

echo "Test demo for a hyperbolic pde (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
    ./ts_demo2 mesh2.m
fi

echo "-----------------------------------------------------------------"
echo "Testing Eigen solver Demos ..."
echo "-----------------------------------------------------------------"
cd ../../solvers/EIGEN
echo "Test demo for computation of eigenvalues of a matrix (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./eigen_demo1 
fi

echo "Test demo for computation of eigenvalues of the Laplace operator (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
    ./eigen_demo2 eigen_demo2.dat
fi

echo "-----------------------------------------------------------------"
echo "Testing Optimization solver Demos ..."
echo "-----------------------------------------------------------------"
cd ../../solvers/OPTIM
echo "Test demo for a one-variable optimization problem (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./opt_demo1
fi

echo "Test demo for a multi-variable optimization problem (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
    ./opt_demo2
fi

echo "Test demo for a pde based optimization problem (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
    ./opt_demo3 test.dat
fi

echo "Test demo for the Brachistochrone optimization problem (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
    ./opt_demo4 30
fi

echo "Test demo for a linear programming problem (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
    ./opt_demo5
fi

echo "Test demo for another linear programming problem (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
    ./opt_demo6
fi

echo "-----------------------------------------------------------------"
echo "Testing Nonlinear Systems solver Demos ..."
echo "-----------------------------------------------------------------"
cd ../../solvers/NLAS
echo "Test demo for a one-variable nonlinear equation with a user defined function (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./nl_demo1
fi

echo "Test demo for a one-variable nonlinear equation defined by a regular expression (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
    ./nl_demo2
fi

echo "Test demo for a multi-variable nonlinear equation with a user defined function (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
    ./nl_demo3
fi

echo "Test demo for a multi-variable nonlinear equation defined by a regular expression (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
    ./nl_demo4
fi

echo "-----------------------------------------------------------------"
echo "Testing Mesh Adaptation Demos ..."
echo "-----------------------------------------------------------------"
cd ../../adapt
echo "Test mesh adaptation demo 1 (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ad1 rect.dom
fi

echo "Test mesh adaptation demo 2 (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ad2 1
fi

cd ../../util/cmesh
echo "================================================================="
echo "Testing OFELI Utilities ..."
echo "-----------------------------------------------------------------"

echo "Test mesh conversion utility (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./cmesh --input disk.bamg --from bamg --to ofeli
fi

echo "Test field conversion utility (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   cd ../cfield
   ./cfield -m cavity.m -i cavity.s -f gmsh
fi

echo "Test mesh generation utility (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   cd ../g2m
   ./g2m -d test.dom -o test.m
fi

cd ..
echo "==========================================================================="
echo "OFELI Testing completed"
echo "==========================================================================="
