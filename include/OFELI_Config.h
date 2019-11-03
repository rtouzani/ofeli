/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

   This file is part of OFELI.

   OFELI is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OFELI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OFELI. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================

                             OFELI's Configuration File

  ==============================================================================*/

#pragma once

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdexcept>
#include "util/constants.h"
#include "util/macros.h"
#include "io/fparser/fparser.h"
#include <complex>

/*! \file OFELI_Config.h
 *  \ingroup Util
 *  \brief File that contains some macros.
 *
 * All these macros can be modified before compiling the library.
 */

/*!
 * \def OFELI_VERSION
 * gives the current version of the library
 */
#define OFELI_VERSION                 "4.0.0"

/*!
 * \def OFELI_RELEASE_DATE
 * gives the date (month-year) of current release
 */
#define OFELI_RELEASE_DATE            "12-2019"

#define MY_RANDOM             55085111

/*! \typedef lsize_t
 *  \ingroup Util
 * \brief This type stands for type \p unsigned \p long
 */
typedef  unsigned long         lsize_t;

/*! \typedef real_t
 *  \ingroup Util
 * \brief This type stands for \p double
 */
typedef  double                real_t;

/*! \typedef complex_t
 *  \ingroup Util
 * \brief This type stands for type \p std::complex<double>
 */
typedef  std::complex<double>  complex_t;


#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define XGRAPH_                0
#undef  INIT_PETSC
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \def GRAPH_MEMORY
 *  \ingroup Mesh
 * \brief Memory necessary to store matrix graph.
 * \details This value is necessary only if nodes are to be renumbered.
 */
#define GRAPH_MEMORY           1000000

/*! \def MAX_NB_EQUATIONS
 *  \ingroup Solver
 * \brief Maximum number of equations
 * \details Useful for coupled problems
 */
#define MAX_NB_EQUATIONS                   5

/*! \def MAX_NB_INPUT_FIELDS
 *  \ingroup Solver
 * \brief Maximum number of fields for an equation
 * \details Useful for coupled problems
 */
#define MAX_NB_INPUT_FIELDS                3

/*! \def MAX_NB_MESHES
 *  \ingroup Solver
 * \brief Maximum number of meshes
 * \details Useful for coupled problems
 */
#define MAX_NB_MESHES                     10

/*! \def MAX_NB_ELEMENTS
 *  \ingroup Mesh
 * \brief Maximal Number of elements
 */
#define MAX_NB_ELEMENTS                10000

/*! \def MAX_NB_NODES
 *  \ingroup Mesh
 * \brief Maximal number of nodes
 */
#define MAX_NB_NODES                   10000

/*! \def MAX_NB_SIDES
 *  \ingroup Mesh
 * \brief Maximal number of sides in.
 */
#define MAX_NB_SIDES                   30000

/*! \def MAX_NB_EDGES
 *  \ingroup Mesh
 * \brief Maximal Number of edges.
 */
#define MAX_NB_EDGES                   30000

/*! \def MAX_NBDOF_NODE
 *  \ingroup Mesh
 * \brief Maximum number of DOF supported by each node
 */
#define MAX_NBDOF_NODE                     6

/*! \def MAX_NBDOF_SIDE
 *  \ingroup Mesh
 * \brief Maximum number of DOF supported by each side
 */
#define MAX_NBDOF_SIDE                     6

/*! \def MAX_NBDOF_EDGE
 *  \ingroup Mesh
 * \brief Maximum number of DOF supported by each edge
 */
#define MAX_NBDOF_EDGE                     2

/*! \def MAX_NB_ELEMENT_NODES
 *  \ingroup Mesh
 * \brief Maximum number of nodes by element
 */
#define MAX_NB_ELEMENT_NODES              20

/*! \def MAX_NB_ELEMENT_EDGES
 *  \ingroup Mesh
 * \brief Maximum number of edges by element
 */
#define MAX_NB_ELEMENT_EDGES              10

/*! \def MAX_NB_SIDE_NODES
 *  \ingroup Mesh
 *  \brief Maximum number of nodes by side
 */
#define MAX_NB_SIDE_NODES                  9

/*! \def MAX_NB_ELEMENT_SIDES
 *  \ingroup Mesh
 *  \brief Maximum number of sides by element
 */
#define MAX_NB_ELEMENT_SIDES               8

/*! \def MAX_NB_ELEMENT_DOF
 *  \ingroup Mesh
 *  \brief Maximum number of dof by element
 */
#define MAX_NB_ELEMENT_DOF                27

/*! \def MAX_NB_SIDE_DOF
 *  \ingroup Mesh
 *  \brief Maximum number of dof by side
 */
#define MAX_NB_SIDE_DOF                    4

/*! \def MAX_NB_INT_PTS
 *  \ingroup Mesh
 *  \brief Maximum number of integration points in element
 */
#define MAX_NB_INT_PTS                    20

/*! \def MAX_NB_MATERIALS
 *  \ingroup Mesh
 *  \brief Maximum number of materials
 */
#define MAX_NB_MATERIALS                  10

/*! \def MAX_NB_PAR
 *  \ingroup IO
 *  \brief Maximum number of parameters
 *  \details Used in class IPF
 */
#define MAX_NB_PAR                        50

/*! \def MAX_ARRAY_SIZE
 *  \ingroup IO
 *  \brief Maximum array size
 *  \details Used in class IPF
 */
#define MAX_ARRAY_SIZE                   100

/*! \def MAX_INPUT_STRING_LENGTH
 *  \ingroup IO
 *  \brief Maximum string length
 *  \details Used in class IPF
 */
#define MAX_INPUT_STRING_LENGTH          100

/*! \def FILENAME_LENGTH
 *  \ingroup IO
 *  \brief Length of a string defining a file name
 */
#define FILENAME_LENGTH                  150

/*! \def MAX_FFT_SIZE
 *  \ingroup IO
 *  \brief Maximal size for the FFT Table 
 *  This table can be used by the FFT for any number of points from 2 up to MAX_FFT_SIZE.
 *  For example, if MAX_FFT_SIZE = 14, then we can transform anywhere from 
 *  2 to 2^15 = 32,768 points, using the same sine and cosine table.
 */
#define MAX_FFT_SIZE                     15

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// Macro
#define XGRAPH_          0
#define ANY              123.45678901e05 // A macro to define any number
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \enum FieldName
 * \brief Enumerate variable that selects field
 */
enum FieldName {
   TEMPERATURE                  =   0,    ///< Temperature scalar field
   DISPLACEMENT                 =   1,    ///< Displacement vector field
   VELOCITY                     =   2,    ///< Velocity vector field
   PRESSURE                     =   3,    ///< Pressure scalar field
   ELECTRIC_FIELD               =   4,    ///< Electric vector field
   MAGNETIC_FIELD               =   5,    ///< Magnetic vector field
   POTENTIAL                    =   6,    ///< Potential vector field
   CONCENTRATION                =   7,    ///< Concentration field
   STRESS                       =   8,    ///< Stress vector field
   INTERNAL_ENERGY              =   9,    ///< Internal energy
   DENSITY                      =  10,    ///< Density
   MOMENTUM                     =  11,    ///< Momentum
};

/*! \enum FieldType
 * \brief Enumerate variable that selects field type
 */
 enum FieldType {
    NONE                = 0, ///< No support information
    NODE_FIELD          = 1, ///< DOFs are supported by nodes
    ELEMENT_FIELD       = 2, ///< DOFs are supported by elements
    SIDE_FIELD          = 3, ///< DOFs are supported by sides (faces in 3-D, edges in 2-D)
    BOUNDARY_SIDE_FIELD = 4, ///< DOFs are supported by boundary sides (faces in 3-D, edges in 2-D)
    EDGE_FIELD          = 5  ///< DOFs are supported by edges
 };

/*! \enum DOFSupport
 * \brief Choose Support of degrees of freedom
 */
enum DOFSupport {
    NODE_DOF          = 1,   ///< DOFs are supported by nodes
    ELEMENT_DOF       = 2,   ///< DOFs are supported by elements
    SIDE_DOF          = 3,   ///< DOFs are supported by sides
    BOUNDARY_SIDE_DOF = 4,   ///< DOFs are supported by sides
    EDGE_DOF          = 5    ///< DOFs are supported by edges
};


/*! \enum ElementShape
 * \brief Enumerate list for element shapes
 */
enum ElementShape {
   NO_ELEMENT    = 0,            /*!< Mesh with no element              */
   POINT         = 1,            /*!< Elements are single points        */
   LINE          = 2,            /*!< Elements are segment lines        */
   TRIANGLE      = 3,            /*!< Elements are triangles            */
   QUADRILATERAL = 4,            /*!< Elements are quadrilaterals       */
   TETRAHEDRON   = 5,            /*!< Elements are tetrahedra           */
   HEXAHEDRON    = 6,            /*!< Elements are hexahedra (bricks)   */
   PENTAHEDRON   = 7,            /*!< Elements are pentahedra (prisms)  */
   PRISM         = 8,            /*!< Elements are prisms               */
   PYRAMID       = 9             /*!< Elements are pyramids             */
};

/*! \brief ExternalFileFormat
 * \brief Enumerate variable that selects external file formats
 */
enum ExternalFileFormat {
   OFELI_FF,    ///< OFELI file format
   GMSH,        ///< Gmsh file format
   GNUPLOT,     ///< Gnuplot file format
   MATLAB,      ///< Matlab m-file
   VTK,         ///< VTK file format
   TECPLOT,     ///< Tecplot file format
   EASYMESH,    ///< Easymesh file format
   GAMBIT,      ///< Gambit file format
   BAMG,        ///< Bamg file format
   NETGEN,      ///< Netgen file format
   TETGEN,      ///< Tetgen file format
   TRIANGLE_FF  ///< Triangle file format
};

/*! \enum ElementType
 * \brief Choose finite element type
 */
enum ElementType {
   LINE2    =  0,            /*!< Line element with 2 nodes (P1) */
   TRIANG3  =  1,            /*!< Triangular element with 3 nodes (P1) */
   QUAD4    =  2,            /*!< Quadrilateral element with 4 nodes (Q1) */
   TETRA4   =  3,            /*!< Tetrahedral element with 4 nodes (P1) */
   HEXA8    =  4,            /*!< Hexahedral element with 8 nodes (Q1) */
   PENTA6   =  5             /*!< Pentahedral element with 6 nodes (P1*Q1) */
};

/*! \enum NonLinearIter
 * Selects iteration method for solving nonlinear problems
 */
enum NonLinearIter {
   BISECTION     =  0,    /*!< Bisection method                          */
   REGULA_FALSI  =  1,    /*!< Regula Falsi method                       */
   PICARD        =  2,    /*!< Picard's iteration method                 */
   SECANT        =  3,    /*!< Secant method                             */
   NEWTON        =  4,    /*!< Newton's method                           */
};

#ifdef WIN32
#define PATH_MATERIAL "c:\\Program Files\\ofeli-"OFELI_VERSION"\\material\\"
#else
#include "datadir.h"
#endif


#if defined(_MSC_VER) && !defined(__MWERKS__)
#define _MSVCPP_ _MSC_VER
#endif

#if defined(OFELI_COMPILED_WITH_DEBUGGING_)
#define _DEBUG
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define OFELI_CONJ std::conj
#define OFELI_ABS  std::abs

#define NOMINMAX

#ifdef WIN32
#define PATH_SEP "\\";
#else
#define PATH_SEP "/";
#endif
#define MATERIAL_EXT      ".md"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @namespace OFELI
 * \brief Namespace OFELI groups all %OFELI library classes, functions and global variables
 */

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/**
 *  \ingroup Global
 *  \brief Verbosity parameter
 *  \details The value of Verbosity can be modified anywhere in the calling programs.
 *  It allows outputting messages in function of the used class or function.
 *  To see how this parameter is used in any class, the OFELI user has to read
 *  corresponding documentation. 
 */
    extern int Verbosity;

/**
 *  \ingroup Global
 *  \brief Time step counter
 *  \details This counter must be initialized by the user if the
 *  macro timeLoop is not used
 *  @remark May be used in conjunction with the macro TimeLoop.
 *  In this case, it has to be initialized before. Its default value is 1
 */
    extern int theStep;

/**
 *  \ingroup Global
 *  \brief Iteration counter
 *  \details This counter must be initialized by the user
 *  @remark May be used in conjunction with the macro IterationLoop.
 *  Its default value is 1
 */
    extern int theIteration;

/**
 *  \ingroup Global
 *  \brief Number of time steps
 *  @remark May be used in conjunction with the macro TimeLoop.
 */
    extern int NbTimeSteps;

/**
 *  \ingroup Global
 *  \brief Maximal number of iterations
 *  @remark May be used in conjunction with the macro IterationLoop.
 *  Its default value is 1000
 */
    extern int MaxNbIterations;

/**
 *  \ingroup Global
 *  \brief Parameter for verbosity of message outputting
 *  \details Its default value is 1
 */
    extern int Verbosity;

/**
 *  \ingroup Global
 *  \brief Time step label
 *  @remark May be used in conjunction with the macro TimeLoop.
 *  In this case, it has to be initialized before
 */
   extern real_t theTimeStep;

/**
 *  \ingroup Global
 *  \brief Time value
 *  @remark May be used in conjunction with the macro TimeLoop.
 *  Its default value is \c 0.0
 */
   extern real_t theTime;

/**
 *  \ingroup Global
 *  \brief Final time value
 *  @remark May be used in conjunction with the macro TimeLoop.
 *  In this case, it has to be initialized before
 */
   extern real_t theFinalTime;

/**
 *  \ingroup Global
 *  \brief Tolerance value for convergence
 *  @remark May be used within an iterative procedure.
 *  Its default value is \c 1.e-8
 */
   extern real_t theTolerance;

/**
 *  \ingroup Global
 *  \brief Value of discrepancy for an iterative procedure
 *  Its default value is \c 1.0
 */
   extern real_t theDiscrepancy;

/**
 *  \ingroup Global
 *  \brief Boolean variable to say if an iterative procedure has converged
 *  \details Its default value is \c false
 */
    extern bool Converged;

/**
 *  \ingroup Global
 *  Boolean to say if PETSc use was initialized.
 *  Useful only if PETSc is used
 */
    extern bool InitPetsc;

/*! @} End of Doxygen Groups */
} /* namespace OFELI */
