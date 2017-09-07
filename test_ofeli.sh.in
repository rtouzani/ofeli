#!/bin/sh

# test_ofeli.sh 
# A script to test OFELI

echo "Testing OFELI Examples, Demos and Utilities ..."

echo "================================================================="
echo "Testing OFELI Examples ..."
echo "-----------------------------------------------------------------"
cd tutorial/lesson1
echo "Test tutorial example 1 (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./lesson1
fi

cd ../lesson2
echo "Test tutorial example 2 (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./lesson2 test.m
fi

cd ../lesson3
echo "Test tutorial example 3 (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./lesson3 test.m
fi

cd ../lesson4
echo "Test tutorial example 4 (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./lesson4 test.dat
fi

cd ../../demos


echo "================================================================="
echo "Testing OFELI Demos ..."

echo "-----------------------------------------------------------------"
echo "Testing Heat Transfer Demos ..."

cd thermal/stdc2
echo "Test demo for 2-D steady state thermal computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./stdc2 proj.dat
   /bin/rm proj.t
fi

cd ../ttd2
echo "Test demo for 2-D transient thermal computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ttd2 proj.dat
   /bin/rm proj.pl
fi

cd ../std3
echo "Test demo for 3-D steady state thermal computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./std3 beam.dat
   /bin/rm beam.sav
fi

echo "-----------------------------------------------------------------"
echo "Testing Solid Mechanics Demos ..."

cd ../../solid/truss
echo "Test demo for 2-D truss structure (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./truss truss.dat
   /bin/rm truss.d truss-1.m
fi

cd ../beam
echo "Test demo for 3-D beams (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./beam beam.dat
  /bin/rm deformed_beam.m beam_tecplot.dat
fi

cd ../lelas2d
echo "Test demo for 2-D linear elasticity computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./lelas2d beam.dat
   /bin/rm beam.d
fi

cd ../contact
echo "Test demo for 2-D linear elasticity with contact computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./contact beam.dat
   /bin/rm beam.d
fi

cd ../lelas3d
echo "Test demo for 3-D linear elasticity computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./lelas3d beam.dat
   /bin/rm beam.d beam-1.m
fi

echo "-----------------------------------------------------------------"
echo "Testing Fluid Dynamics Demos ..."

cd ../../fluid/tins2
echo "Test demo for 2-D transient incompressible fluid flow computation (Projection method) (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./tins2 cavity.dat
   /bin/rm cavity.v cavity.p
fi

cd ../tiff2
echo "Test demo for 2-D transient incompressible fluid flow computation (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./tiff2 cavity.dat
   /bin/rm cavity.v cavity.p
fi

cd ../lh2d
echo "Test demo for linear hyperbolic equation in 2-D (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./lh2d proj.dat
   /bin/rm rectangle.t
fi

cd ../euler-2d
echo "Test demo for 2-D Compressible Euler equations (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./euler2d proj.dat
   /bin/rm rect.*
fi

echo "-----------------------------------------------------------------"
echo "Testing Electromagnetic Demos ..."

cd ../../electromagnetics/Helmholtz
echo "Test demo for Helmholtz equation in bounded media (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./helmholtz proj.dat
fi

cd ../EddyCurrent2d-1
echo "Test demo for 2-D Eddy Current computation in bounded media (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ec2d1 proj.dat
fi

echo "-----------------------------------------------------------------"
echo "Testing Linear System solver Demos ..."

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
    /bin/rm *.vtk mm-*.m
fi

echo "-----------------------------------------------------------------"
echo "Testing Eigen solver Demos ..."

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
    /bin/rm *.pos
fi

echo "-----------------------------------------------------------------"
echo "Testing Optimization solver Demos ..."

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

echo "-----------------------------------------------------------------"
echo "Testing Mesh Adaptation Demos ..."

cd ../../adapt
echo "Test mesh adaptation demo 1 (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ad1 rect.dom
    /bin/rm *.pos *.m
fi

echo "Test mesh adaptation demo 2 (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ./ad2 L.dom
fi

cd ../../
echo "================================================================="
echo "Testing OFELI Utilities ..."
echo "-----------------------------------------------------------------"

cd util/conv/examples
echo "Test mesh conversion utility (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ../src/cmesh --input disk.bamg --from bamg --to ofeli
   /bin/rm disk.m
fi

echo "Test field conversion utility (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ../src/cfield -m cavity.m -i cavity.s -f gmsh
   /bin/rm cavity.pos
fi

cd ../../g2m/examples
echo "Test mesh generation utility (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
   ../src/g2m -d test.dom -o test.m
   /bin/rm test.m
fi

cd ../..
echo "==========================================================================="
echo "OFELI Testing completed"
