/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2024 Rachid Touzani

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

                          Implementation of class 'help'

  ==============================================================================*/

#include "help.h"
#include "helps.h"

using std::string;
using std::cout;
using std::endl;
using std::vector;

namespace RITA {

help::help(rita* r,
           cmd*  command)
     : _rita(r), _cmd(command)
{ }


int help::run(int opt)
{
   _nb_args = _cmd->getNbArgs();
   if (_nb_args==0) {
      if (opt==0)
         cout << H0 << endl;
      else if (opt==1)
         cout << H2 << endl;
      cout << "Type: help <keyword> to display help about a keyword.\n\n";
      cout << "Available help topics are:\n" << Help_Topics << endl;
      return 0;
   }
   static const vector<string> kw {"help","lic$ense","set","load","end","exit","quit",
                                   "vect$or","mat$rix","tab$ulation","mesh","grid","func$tion",
                                   "sample","plot","hist$ory","data","list","desc$ription",
                                   "approx$imation","integr$ation","ls","algebraic","ae","ode",
                                   "pde","optim$ization","eigen","solve","rename",
                                   "remove","delete","calc$ulator","print"};
   _cmd->set(kw);
   CHK_MSGR(_nb_args>2,_pr,"Only two arguments are allowed for command 'help'")
   switch (_cmd->getArg(" ")) {

      case   0:
         HELP_ARG2
         cout << "Command 'help' to display help with or without argument." << endl;
         break;

      case   1:
         HELP_ARG2
         cout << "Command 'license' to display rita licensing agreement." << endl;
         break;

      case   2:
         cout << "Command 'set' to set various configuration parameters." << endl;
         cout << Config_Help << endl;
         break;

      case   3:
         HELP_ARG2
         cout << "Command 'load' enables choosing a script file for rita." << endl;
         cout << "load <file>" << endl;
         break;

      case   4:
         HELP_ARG2
         cout << "Command 'end' or '<' brings back to command level." << endl;
         break;

      case   5:
      case   6:
         HELP_ARG2
         cout << "End execution of program rita." << endl;
         break;

      case   7:
         HELP_ARG3("vector")
         Vector();
         break;

      case   8:
         HELP_ARG3("matrix")
         Matrix();
         break;

      case   9:
         HELP_ARG3("tabulation")
         cout << "Command 'tabulation'" << endl;
         cout << Tab_Help << endl;
         break;

      case  10:
         HELP_ARG3("mesh")
         mesh();
         break;

      case  11:
         HELP_ARG3("grid")
         cout << "Command 'grid' enables defining a 1-D, 2-D or 2-D uniform grid. Synopsis:\n\n";
         cout << Grid_Help << endl;
         break;

      case  12:
         HELP_ARG3("function")
         cout << "Command 'function' enables defining analytically a function. Synopsis:\n\n";
         cout << Fct_Help << endl;
         break;

      case  13:
         cout << "Command 'sample' enables sampling a function, i.e. defining its value on a given grid or mesh.\n";
         cout << "Synopsis:\n\n" << Sample_Help << endl;
         break;

      case  14:
         plot();
         break;

      case  15:
         cout << "Command 'history' enables adding a given vector (e.g. solution of ode, pde, ... at given time step)\n";
         cout << "to a given (or to define) history vector. Synopsis:\n\n" << HVect_Help << endl;
         break;

      case  16:
         cout << "Command 'data' is useful to list all defined data and entities in a rita execution. It has no arguments." << endl;
         break;

      case  17:
         cout << "Command 'list' is useful to list specified data or entities that has been defined by the user.\nSynopsis:\n\n";
         cout << "list <e>\n\n";
         cout << "e must have one of the values: param, grid, mesh, vector, hvector, function, tabulation,\n";
         cout << "matrix, ls, algebraic, ode, pde, optimization, eigen" << endl;
         break;

      case  18:
         cout << "Command 'description' enables associating a description of a given datum or entity. Synopsis:\n\n";
         cout << "description <text>" << endl;
         break;

      case  19:
         cout << "Command 'approximation'" << endl;
         cout << Approx_HHelp << endl;
         break;

      case  20:
         cout << "Command 'integration' enables computing a numerical integration to approximate the defined integral of a function.\nSynopsis:" << endl;
         cout << Int_Help << endl;
         break;

      case  21:
         ls();
         break;

      case  22:
      case  23:
         AE();
         break;

      case  24:
         ODE();
         break;

      case  25:
         PDE();
         break;

      case  26:
         Opt();
         break;

      case  27:
         Eig();
         break;

      case  28:
         solve();
         break;

      case  29:
         cout << "Command 'rename' enables renaming data or entities. Synopsis:\n\n";
         cout << "rename <name>\n";
         cout << "name is the name of the entity to rename. rita indicates if the entity does not exist." << endl;
         break;

      case  30:
      case  31:
         cout << "Command 'remove' enables deleting data or entities. Note that for efficiency reasons,\n";
         cout << "specified data are not really deleted. Internally, they are marked as inactive. This is transparent to\n";
         cout << "the user. Synopsis:\n\n";
         cout << "remove <name>\n";
         cout << "name is the name of the entity to remove. rita indicates if the entity does not exist." << endl;
         break;

      case  32:
         Calc();
         break;

      case  33:
         cout << "Command 'help' \n" << Print_Help << endl;
         break;

      default:
         MSGR(_pr,"Unrecognized help topic: "+_cmd->getToken()+"\nHelp topics are:\n"+Help_Topics)
   }
   return 0;
}


int help::Vector()
{
   cout << "Command 'vector' enables defining a vector. A vector in rita can be used as datum\n";
   cout << "for a system of equations, a solution of a PDE defined at nodes, etc.\n";
   cout << "A vector can be defined by one of the three methods:\n";
   cout << "- By calculator. For instance, the command\n";
   cout << "  v = {1.0,-1.2,2.3}\n";
   cout << "  defines a row vector with dimension 3, with given entries.\n";
   cout << "  v = vector(3)\n";
   cout << "  defines a column vector with dimension 3. Its entries can be defined as v[0], v[1], v[2]\n";
   cout << "- From a solver. For instance when defining a partial differential equation,\n";
   cout << "  declaring a variable name defines it as a vector in rita and size it.\n";
   cout << "- By using the present command.\n" << "Synopsis:\n\n" << Vect_Help << endl;
   return 0;
}


int help::Calc()
{
   cout << "rita can be used as a simple calculator to manipulate various data. Calculator instructions can be "
           "introduced like any command.\n\n";
   cout << "- Variables can be introduced like in the following examples:\n";
   cout << "  x = 10\n  y = 2*x - exp(x)\n  z = y^2-sin(pi*x)\n";
   cout << "  x, y, z are then stored as `parameters'.\n\n";
   cout << "- Vectors can also be defined as in:\n";
   cout << "    v = {1.0,-2.0,6.0}\n    w = {0.0,-1.,2.0}\n    A = v'*w\n    z = v*w'\n";
   cout << "  Here v and w are row vectors, v' is the transpose of v. Hence A is a 3x3-matrix and z is a scalar.\n";
   cout << "  Some functions are defined for vectors:\n";
   cout << "    p1 = norm1(v)\n    p2 = norm2(v)\n    pI = normI(v)\n";
   cout << "  Here p1, p2, pI are respectively the 1-norm, the 2-norm (euclidean) and the max-norm\n\n";
   cout << "- Matrices can be defined as in:\n    M=matrix(2,3)\n  Here M is a 2x3-matrix and N is a 4x4-matrix\n";
   cout << "  Entries of these matrices can be refered to or defined by:\n";
   cout << "    M[i,j] = 3\n  Here the index i runs from 0 to 1, and j from 0 to 2.\n";
   cout << "  Beside classical algebraic operations on matrices (sums, products, matrix-vector products, ...), other\n";
   cout << "  are available:\n";
   cout << "    M = eye(2)\n    N = ones(4,2)\n    A = diag(2,-1)\n    K = laplace1d(3,2.0)\n";
   cout << "  Here M is the identity 2x2-matrix, N is the 4x2-matrix with all entries equal to 1.\n";
   cout << "  A is the 2x2 diagonal matrix with diagonal entries equal to -1.\n";
   cout << "  K is the matrix of 3-stencil finite difference approximation of the operator (-u\") multiplied by 2.0, i.e.\n";
   cout << "          |  4.0  -2.0  0.0 |\n     K =  | -2.0   4.0 -2.0 |\n          |  0.0  -2.0  4.0 |\n\n";
   return 0;
}


int help::Matrix()
{
   cout << "Command 'matrix' enables defining a matrix.\nA matrix can be defined by one of the three methods:\n";
   cout << "- By calculator. For instance, the command\n";
   cout << "  M = matrix(2,3)\n";
   cout << "  defines a 2x3-matrix initialized by <span class=var>0</span>. Once declared, an entry of the matrix\n";
   cout << "  can be given for instance by\n";
   cout << "  M[0,2]=5\n";
   cout << "  Note that the first row, first column entry is M[0,0], and generally the n-th row, m-th column entry\n";
   cout << "  is M[n-1,m-1]\n";
   cout << "- By using the present command.\nSynposis:\n" << Matrix_Help << endl;
   return 0;
}


int help::mesh()
{
   if (_nb_args==1) {
      cout << "Command 'mesh' constructs a finite element mesh. More precisely, it prepares data to construct the mesh\n";
      cout << "using gmsh. The command has no argument. It creates a submenu with the following subcommands:\n\n";
      cout << Mesh_Help << endl;
      cout << "Type: 'help mesh subcommand' to get help on the chosen subcommand." << endl; 
      return 0;
   }
   static const vector<string> kw {"1d","rect$angle","cube","point","curve","surface","volume","contour","code","read"};
   switch (_cmd->getKW(kw)) {

      case   0:
         return set1D();

      case  1:
         return setRectangle();

      case  2:
         return setCube();

      case  3:
         return setPoint();

      case  4:
         return setCurve();

      case  5:
         return setSurface();

      case  6:
         return setVolume();

      case  7:
         return setContour();

      case  8:
         return setCode();

      case  9:
         return MeshRead();

      default:
         MSGR(_pr,"Unknown mesh sub-command.")
   }
   return 0;
}


int help::set1D()
{
   cout << "Subcommand '1d' constructs a 1-D mesh. Synopsis:\n\n";
   cout << Mesh_1d_Help << endl;
   return 0;
}


int help::setRectangle()
{
   cout << "Subcommand 'rectangle' constructs triangulation of a rectangle. Synopsis:\n\n";
   cout << Mesh_Rect_Help << endl;
   return 0;
}


int help::setCube()
{
   cout << "Subcommand 'cube' constructs mesh of a cube-like domain. Synopsis:\n\n";
   cout << Mesh_cube_Help << endl;
   return 0;
}


int help::setPoint()
{
   cout << "Subcommand 'point' defines a point in the domain to mesh. Synopsis:\n\n";
   cout << Mesh_point_Help << endl;
   return 0;
}


int help::setCurve()
{
   cout << "Subcommand 'curve' defines a curve in the domain to mesh. Synopsis:\n\n";
   cout << Mesh_curve_Help << endl;
   return 0;
}


int help::setSurface()
{
   cout << "Subcommand 'surface' defines a surface in the domain to mesh. Synopsis:\n\n";
   cout << Mesh_surf_Help << endl;
   return 0;
}


int help::setVolume()
{
//   cout << Mesh_volume_Help << endl;
   return 0;
}


int help::setContour()
{
   cout << "Subcommand 'contour' defines a contour in the domain to mesh. Synopsis:\n\n";
   cout << Mesh_contour_Help << endl;
   return 0;
}


int help::setCode()
{
   cout << "Subcommand 'code' enables assigning values to codes associated to specific points, curves, surfaces or volumes.\n";
   cout << "Synopsis:\n\n" << Mesh_code_Help << endl;
   return 0;
}


int help::MeshRead()
{
   cout << "Subcommand 'read' enables reading mesh data from a given file. Synopsis:\n\n";
   cout << Mesh_read_Help << endl;
   return 0;
}


int help::ls()
{
   cout << "Command 'ls' enables defining a linear system of equations.\nSynopsis:\n\n";
   cout << LS_Help << endl;
   return 0;
}


int help::AE()
{
   cout << "Command 'algebraic' (or simply 'ae') enables solving (systems of) algebraic equations. The command has a\n";
   cout << "an extended form. In its short version, the Synopsis is:\n\n";
   cout << AE_Help << endl;
   cout << "In its extended version, the command has no arguments but a subcommand menu provides the following commands:\n\n";
   cout << AE_HHelp << endl;
   return 0;
}


int help::ODE()
{
   cout << "Command 'ode' enables solving (systems of) ordinary differential equations. The command has a\n";
   cout << "an extended form. In its short version, the synopsis is:\n\n";
   cout << ODE_Help << endl;
   cout << "In its extended version, the command has no arguments but a subcommand menu provides the following commands:\n\n";
   cout << ODE_HHelp << endl;
   return 0;
}


int help::PDE()
{
   cout << "Command 'pde' enables defining a partial differential equation to solve later. The list of\n";
   cout << "available partial differential equations is in constant growing.\n";
   cout << "The command 'pde' has one argument, which is the name of the pde to define. The list of\n";
   cout << "available pde names is at the end of this text.";
   cout << "Once the command 'pde' entered with its argument.\n";
   cout << "A submenu is available with the following subcommands:" << endl;
   cout << PDE_Help << endl;
   return 0;
}


int help::Eig()
{
   cout << "Command 'eigen' enables solving an eigenvalue problem. Synopsis:\n" << endl;
   cout << Eig_Help << endl;
   return 0;
}


int help::Opt()
{
   cout << "Command 'optimization' enables defining an optimization problem. The command can be given\n";
   cout << "either in a short or an extended form.\n";
   cout << "In its short version, the command synopsis is:\n\n";
   cout << Opt_Help << endl;
   cout << "In its extended version, the command has no arguments but a subcommand menu provides the following commands:\n\n";
   cout << Opt_HHelp << endl;
   return 0;
}


int help::solve()
{
   cout << "Command 'solve' launches the solution of a selection of problems among defined once.\n";
   cout << "The command has short and extended versions. For the short version, the synopsis is:\n\n";
   cout << Solve_Help << endl;
   cout << "In its extended version, the command has no arguments but a subcommand menu provides the following commands:\n\n";
   cout << Solve_HHelp << endl;
   return 0;
}


int help::plot()
{
   cout << "Command 'plot' enables plotting specific data. Various types of data can be plotted\n";
   cout << "using various plotting tools. The command can be invoked either in its short or its\n";
   cout << " extended version.\n";
   cout << "In its short version, the synopsis is:\n\n";
   cout << Plot_Help << endl;
   cout << "In its extended version, the synopsis is:\n\n";
   cout << Plot_HHelp << endl;
   return 0;
}

} /* namespace RITA */