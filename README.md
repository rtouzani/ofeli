OFELI
An Object Finite Element Library

Copyright (c) 1998-2024 Rachid Touzani


1. Introduction
2. Documentation
3. Installation
   - Unix Installation
   - Windows Installation

1. Introduction
---------------

OFELI is an object oriented library of C++ classes for development of finite element codes.
Its main features are:
* Various storage schemes of matrices (dense, sparse, skyline).
* Direct methods of solution of linear systems of equations as well as
  various combinations of iterative solvers and preconditioners.
* Shape functions of most "popular" finite elements
* Element arrays of most popular problems (Heat Transfer, Fluid Flow, Solid
  Mechanics, Electromagnetics, ...).

The OFELI package is not only a library of classes for Finite Element developments. The package
contains in addition :

* Very simple examples for "academic" finite element codes
* More elaborated codes for various types of problems
* Utility programs to convert mesh and output files and to generate simple meshes

OFELI is free but Copyrighted software. It is distributed under the
terms of the GNU General Public License (GPL). Details of the GPL are in
the file 'COPYING' that comes with the OFELI software package.


2. Installation
---------------

The OFELI package can be installed on Linux, Unix and Windows (95/98/NT/2000/XP) systems.

Unix Installation:

The following instructions apply to UNIX (like) systems only.
The only thing you need is a C++ compiler.

A clean installation can be performed in the following steps:
* 'gunzip' and 'untar' the downloaded file.
* Edit the file OFELI_Config.h in directory 'include' and modify the line defining 'PATH_MATERIAL'
  according to your directory installation. This procedure will be automatized in the forthcoming
  releases.
* Execute the configuration script by typing:
     ./configure
  A first execution with argument --help displays a short manual of the script.
  A typical execution looks like:
     ./configure --enable-release --prefix=/home/me/ofeli --libdir=/home/me/lib --bindir=/home/me/bin CXXFLAGS=-O4


Windows Installation :

To install on Windows system, you must have the Visual C++ compiler.
The installation procedure depends on the file you have downloaded.

* If this is a zip file then you must unzip it, go the subdirectory windows where you can find
  Visual C++ workspace to install the library, utilities, examples, and demos.
  You can also install the package by using the makefile by executing
     nmake -f makefile.vc
  Here also the Visual C++ compiler is required.

* If this is an exe file, than this is a standard windows setup installator. The library and
  the examples and demo files are already compiled. You can also recompile using Visual C++
  workspace.

* Edit the file OFELI_Config.h in directory 'include' and modify the line defining 'PATH_MATERIAL'
  according to your directory installation. This procedure will be automatized in the forthcoming
  releases.

