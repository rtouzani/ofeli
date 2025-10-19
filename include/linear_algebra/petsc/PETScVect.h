/*==============================================================================

                                    O  F  E  L  I
                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

                 Definition of Template class PETScVect for vectors
                               using the library PETSc

  ==============================================================================*/


#ifndef __PETSC_VECT_H
#define __PETSC_VECT_H

#ifdef USE_PETSC
#include <petsc.h>
#include <petscvec.h>

#include "OFELI_Config.h"
#include "util/macros.h"
#include "mesh/Mesh.h"
#include "mesh/MeshUtil.h"
#include "mesh/Grid.h"
#include "linear_algebra/Point.h"
#include "linear_algebra/LocalVect.h"
#include "linear_algebra/LocalMatrix.h"
#include "io/Fct.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Tetra4.h"
#include "util/util.h"
#include "linear_algebra/Vect.h"


/*! \file PETScVect.h
 *  \brief Definition file for class PETScVect.
 */

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! 
 *  \class PETScVect
 *  \ingroup VectMat
 *  \brief To handle general purpose vectors using Petsc
 *
 * This template class enables considering vectors of various data types.
 * Operators \b =, \b [] and \b () are overloaded so that one can write for instance:
 *
 * \verbatim
    PETScVect<double> u(10), v(10);
    v = -1.0;
    u = v;
    u.set(3,-2.0);
   \endverbatim
 *
 * to set vector \b v entries to \b -1, copy vector \b v into vector \b u and assign
 * third entry of \b v to \b -2. Note that entries of \b v are here \b v(1), \b v(2),
 * ..., \b v(10), <i>i.e.</i> vector entries start at index \b 1.
 *
 * @remark A PETScVect instance can be 1-D, 2-D or 3-D, <tt>i.e.</tt> one can have 1, 2 or 
 * 3 indices. This is set while the vector is constructed. This can be helpful for instance 
 * in the case of a structured grid.
 *
 * @warning This class is available only when OFELI has been installed with Petsc
 * In this case, only vectors used for building and solving linear systems need to be
 * instances of PETScVect. 
 *
 * \tparam T_ Data type (double, int, complex<double>, ...)
 */

template<class T_> class Vect;
using std::to_string;


template<class T_>
class PETScVect
{

 public:

/// \brief Default Constructor.
/// Initialize a zero-length vector
    PETScVect();

/// \brief Constructor setting vector size.
/// @param [in] n Size of vector
    PETScVect(size_t n);

/** \brief Constructor of a 2-D index vector.
 *  \details This constructor can be used for instance for a 2-D grid vector
 *  @param [in] nx Size for the first index
 *  @param [in] ny Size for the second index
 *  @remark The size of resulting vector is nx*ny
 */
    PETScVect(size_t nx,
              size_t ny);

/** \brief Constructor of a 3-D index vector.
 *  \details This constructor can be used for instance for a 3-D grid vector
 *  @param [in] nx Size for the first index
 *  @param [in] ny Size for the second index
 *  @param [in] nz Size for the third index
 *  @remark The size of resulting vector is <tt>nx*ny*nz</tt>
 */
    PETScVect(size_t nx,
              size_t ny,
              size_t nz);

/** \brief Constructor of a 4-D index vector.
 *  \details This constructor can be used for instance for a 4-D grid vector
 *  @param [in] nx Size for the first index
 *  @param [in] ny Size for the second index
 *  @param [in] nz Size for the third index
 *  @param [in] nt Size for the fourth index
 *  @remark The size of resulting vector is nx*ny*nz*nt
 */
    PETScVect(size_t nx,
              size_t ny,
              size_t nz,
              size_t nt);

/** \brief Create an instance of class PETScVect as an image of a C/C++ array.
 *  @param [in] n Dimension of vector to construct
 *  @param [in] x C-array to copy
 */
    PETScVect(size_t n,
              T_*    x);

/** \brief Constructor with a Grid instance
 *  \details The constructed vector has as size the total number of grid nodes
 *  @param [in] g Grid instance
 */
    PETScVect(Grid& g);

/** \brief Constructor of a MPI vector using its global size.
 *  @param [in] comm Communicator which represents all the processs that PETSc knows about
 *  @param [in] n Global size of vector
 */
    PETScVect(MPI_Comm comm,
              size_t   n);

/** \brief Constructor with a mesh instance
 *  @param [in] m Mesh instance
 *  @param [in] nb_dof Number of degrees of freedom per node, element or side
 *  If <tt>nb_dof</tt> is set to <tt>0</tt> (default value) the constructor picks this number from
 *  the Mesh instance
 *  @param [in] dof_type Type of degrees of freedom. To be given among the enumerated
 *  values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>, <tt>SIDE_DOF</tt> or <tt>EDGE_DOF</tt>
 *  [Default: <tt>NODE_DOF</tt>]
 */
    PETScVect(Mesh& m,
              int   nb_dof=0,
              int   dof_type=NODE_DOF);

/** \brief Constructor with a mesh instance giving name and time for vector
 *  @param [in] m Mesh instance
 *  @param [in] name Name of the vector
 *  @param [in] t Time value for the vector
 *  @param [in] nb_dof Number of degrees of freedom per node, element or side
 *  If <tt>nb_dof</tt> is set to <tt>0</tt> the constructor picks this number from the Mesh instance
 *  @param [in] dof_type Type of degrees of freedom. To be given among the enumerated
 *  values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>, <tt>SIDE_DOF</tt> or <tt>EDGE_DOF</tt>
 *  [Default: <tt>NODE_DOF</tt>]
 */
    PETScVect(Mesh&  m,
              string name,
              real_t t=0.0,
              int    nb_dof=0,
              int    dof_type=NODE_DOF);

/** \brief Constructor using boundary conditions.
 *  \details Boundary condition values contained in <tt>bc</tt> are reported to vector <tt>v</tt>
 *  @param [in] v PETScVect instance to update
 *  @param [in] bc PETScVect instance containing imposed valued at desired DOF
 */
    PETScVect(const PETScVect<T_>& v,
              const PETScVect<T_>& bc);

/** \brief Constructor to select some components of a given vector.
 *  @param [in] v PETScVect instance to extract from
 *  @param [in] nb_dof Number of DOF to extract
 *  @param [in] first_dof First DOF to extract
 *  For instance, a choice <tt>first_dof=2</tt> and <tt>nb_dof=1</tt> means that
 *  the second DOF of each node is copied in the vector
 */
    PETScVect(const PETScVect<T_>& v,
              size_t               nb_dof,
              size_t               first_dof);

/// \brief Copy constructor
    PETScVect(const PETScVect<T_>& v);

/** \brief Constructor to select one component from a given 2 or 3-component vector
 *  @param [in] v PETScVect instance to extract from
 *  @param [in] n Component to extract (must be > 1 and < 4 or).
 */
    PETScVect(const PETScVect<T_>& v,
              size_t               n);

/** \brief Constructor that extracts some degrees of freedom (components) from given instance of PETScVect
 *  \details This constructor enables constructing a subvector of a given
 *  PETScVect instance. It selects a given list of degrees of freedom
 *  and put it according to a given order in the instance to construct.
 *  @param [in] d Integer number giving the list of degrees of freedom. This
 *  number is made of <tt>n</tt> digits where <tt>n</tt> is the number of degrees of freedom.
 *  Let us give an example: Assume that the instance <tt>v</tt> has 3 DOF by entity (node, element 
 *  or side). The
 *  choice <tt>d=201</tt> means that the constructed instance has 2 DOF where the first
 *  DOF is the third one of <tt>v</tt>, and the second DOF is the first one of f <tt>v</tt>.
 *  Consequently, no digit can be larger than the number of DOF the constructed instance.
 *  In this example, a choice <tt>d=103</tt> would produce an error message.
 *  @param [in] v PETScVect instance from which extraction is performed.
 *  @param [in] name Name to assign to vector instance [Default: " "].
 *  \warning Don't give zeros as first digits for the argument <tt>d</tt>. The number is
 *  in this case interpreted as octal !!
*/
    PETScVect(size_t               d,
              const PETScVect<T_>& v,
              const string&        name=" ");

/// \brief Destructor
    ~PETScVect();

/** \brief Initialize vector with a c-array
 *  @param [in] v c-array (pointer) to initialize PETScVect
 *  @param [in] n size of array
 */
    void set(const T_* v,
             size_t    n);

/** \brief Initialize a local vector using MPI
 *  @param [in] comm
 *  @param [in] n local size of vector
 *  @param [in] N global size of vector
 */
    void setMPI(MPI_Comm comm,
                size_t   n,
                size_t   N);

/** \brief Initialize vector with another PETScVect instance
 *  @param [in] v PETScVect instance to extract from
 *  @param [in] nb_dof Number of DOF per node, element or side (By default, 0: Number of degrees
 *  of freedom extracted from the Mesh instance)
 *  @param [in] first_dof First DOF to extract (Default: 1)
 *  For instance, a choice <tt>first_dof=2</tt> and <tt>nb_dof=1</tt> means that
 *  the second DOF of each node is copied in the vector
 */
    void select(const PETScVect<T_>& v,
                size_t               nb_dof=0,
                size_t               first_dof=1);

/** \brief Initialize vector with an algebraic expression
 *  @param [in] exp Regular algebraic expression that defines a function of x, y and z
 *  which are coordinates of nodes.
 *  @param [in] dof Degree of freedom to which the value is assigned [Default: <tt>1</tt>]
 */
    void set(const string& exp,
             size_t        dof=1);

/** \brief Initialize vector with an algebraic expression with providing mesh data
 *  @param [in] ms Mesh instance
 *  @param [in] exp Regular algebraic expression that defines a function of x, y and z
 *  which are coordinates of nodes.
 *  @param [in] dof Degree of freedom to which the value is assigned [Default: <tt>1</tt>]
 */
    void set(Mesh&         ms,
             const string& exp,
             size_t        dof=1);

/** \brief Initialize vector with an algebraic expression 
 *  @param [in] x PETScVect instance that contains coordinates of points
 *  @param [in] exp Regular algebraic expression that defines a function of x and i
 *  which are coordinates of nodes and indices starting from <tt>1</tt>.
 */
    void set(const PETScVect<real_t>& x,
             const string&            exp);

/** \brief Define mesh class to size vector
 *  @param [in] m Mesh instance
 *  @param [in] nb_dof Number of degrees of freedom per node, element or side
 *  If <tt>nb_dof</tt> is set to <tt>0</tt> the constructor picks this number from the Mesh instance
 *  @param [in] dof_type Parameter to precise the type of degrees of freedom. To be chosen
 *  among the enumerated values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>,
 *  <tt>SIDE_DOF</tt>, <tt>EDGE_DOF</tt> [Default: <tt>NODE_DOF</tt>]
 */
    void setMesh(Mesh& m,
                 int   nb_dof=0,
                 int   dof_type=NODE_DOF);

/** \brief Define grid class to size vector
 *  @param [in] g Grid instance
 */
    void setGrid(Grid& g);

/// \brief Return vector (global) size
    size_t size() const { return _size; }

/// \brief Return vector local size
/// \details Local size is the size on the current processor
    PetscInt getLocalSize() const;

/** \brief Set vector size (for 1-D, 2-D or 3-D cases)
 *  \details This function allocates memory for the vector but does not initialize its components
 *  @param [in] nx Number of grid points in <tt>x</tt>-direction
 *  @param [in] ny Number of grid points in <tt>y</tt>-direction [Default: <tt>1</tt>]
 *  @param [in] nz Number of grid points in <tt>z</tt>-direction [Default: <tt>1</tt>]
 */
    void setSize(size_t nx,
                 size_t ny=1,
                 size_t nz=1);

/** \brief Set vector size
 *  \details This function allocates memory for the vector but does not initialize its components
 *  @param [in] n Size of vector
 */
    void resize(size_t n);

/** \brief Set vector size and initialize to a constant value
 *  \details This function allocates memory for the vector
 *  @param [in] n Size of vector
 *  @param [in] v Value to assign to vector entries
 */
    void resize(size_t n,
                T_     v);

/** \brief Set DOF type of vector
 *  \details The DOF type combined with number of DOF per component enable
 *  determining the size of vector
 *  @param [in] dof_type Type of degrees of freedom. Value to be chosen among 
 *  the enumerated values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>, <tt>SIDE_DOF</tt> 
 *  or <tt>EDGE_DOF</tt>
 */
    void setDOFType(int dof_type) { _dof_type = dof_type; }

/// \brief Return vector number of degrees of freedom
    size_t getNbDOF() const { return _nb_dof; }

/// \brief Return vector number of entities (nodes, elements or sides)
    size_t getNb() const { return _size/_nb_dof; }

/// \brief Return Mesh instance
    Mesh &getMesh() const { return *_theMesh; }

/// \brief Return <tt>true</tt> if vector contains a Mesh pointer, <tt>false</tt> if not
/// \details A PETScVect instance can be constructed using mesh information 
    bool WithMesh() const { return _with_mesh; }

/// \brief Return <tt>true</tt> if vector contains a Grid pointer, <tt>false</tt> if not
/// \details A PETScVect instance can be constructed using grid information 
    bool WithGrid() const { return _with_grid; }

/** Return DOF type of vector
 *  @return dof_type Type of degrees of freedom. Value among the enumerated
 *  values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>, <tt>SIDE_DOF</tt> or <tt>EDGE_DOF</tt>
 */
    int getDOFType() const { return _dof_type; }

/// \brief Set time value for vector
    void setTime(real_t t) { _time = t; }

/// \brief Get time value for vector
    real_t getTime() const { return _time; }

/// \brief Set name of vector
    void setName(string name) { _name = name; }

/// \brief Get name of vector
    string getName() const { return _name; }

/// \brief Calculate 1-norm of vector
    PetscScalar getNorm1() const;

/** \brief Compute a norm of vector
 *  @param [in] t Norm type to compute: To choose among enumerate values:
 *                NORM1: 1-norm
 *                WNORM1: Weighted 1-norm (Discrete L1-norm)
 *                NORM2: 2-norm
 *                WNORM2: Weighted 2-norm (Discrete L2-norm)
 *                NORM_MAX: max norm (Infinity norm)
 *  @return Value of norm
 *  @warning This function is available for real valued vectors only
 */
    real_t Norm(NormType t) const;

/// \brief Calculate 2-norm (Euclidean norm) of vector
    PetscScalar getNorm2() const;

/// \brief Calculate weighted 1-norm of vector
/// The wighted 1-norm is the 1-Norm of the vector divided by its size
    PetscScalar getWNorm1() const { return getNorm1()/_size; }

/** \brief Calculate weighted 2-norm of vector
 *  \details The weighted 2-norm is the 2-Norm of the vector divided by the
 *  square root of its size
 */
    PetscScalar getWNorm2() const { return getNorm2()/sqrt(real_t(_size)); }

/// \brief Calculate Max-norm (Infinite norm) of vector
    PetscScalar getNormMax() const;

/// \brief Calculate Min value of vector entries
    T_ getMin() const;

/// \brief Calculate Max value of vector entries
    T_ getMax() const;

/// \brief Return number of grid points in the <tt>x</tt>-direction if grid indexing is set
    size_t getNx() const { return _nx; }

/// \brief Return number of grid points in the <tt>y</tt>-direction if grid indexing is set
    size_t getNy() const { return _ny; }

/// \brief Return number of grid points in the <tt>z</tt>-direction if grid indexing is set
    size_t getNz() const { return _nz; }

/// \brief Return number of grid points in the <tt>t</tt>-direction if grid indexing is set
    size_t getNt() const;

/** \brief Assign a given function (given by an interpretable algebraic expression) of indices 
 *  components of vector.
 *  \details This function enable assigning a value to vector entries as function of indices 
 *  @param [in] exp  Regular algebraic expression to assign. It must involve the variables <tt>i</tt>,
 *  <tt>j</tt> and/or <tt>k</tt>.
 */
    void setIJK(const string& exp);

/** \brief Assign a given function (given by an interpretable algebraic expression) of indices 
 *  components of vector.
 *  \details This function enable assigning a value to vector entries as function of indices 
 *  @param [in] exp  Regular algebraic expression to assign. It must involve the variables <tt>i</tt>,
 *  <tt>j</tt>, <tt>k</tt> and/or <tt>l</tt>.
 */
    void setIJKL(const string& exp);

/** \brief Assign a given value to components of vector with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] m Instance of mesh
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] val Value to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned [default: <tt>1</tt>]
 */
    void setNodeBC(Mesh&  m,
                   int    code,
                   T_     val,
                   size_t dof=1);

/** \brief Assign a given function (given by an interpretable algebraic expression) to 
 *  components of vector with given code.
 *  \details Vector components are assumed nodewise
 *  @param [in] m    Instance of mesh
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned [default: <tt>1</tt>]
 */
    void setNodeBC(Mesh&         m,
                   int           code,
                   const string& exp,
                   size_t        dof=1);

/** \brief Assign a given value to components of vector with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] val  Value to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned [Default: <tt>1</tt>]
 */
    void setNodeBC(int    code,
                   T_     val,
                   size_t dof=1);

/** \brief Assign a given function (given by an interpretable algebraic expression) to 
 *  components of vector with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned [Default: <tt>1</tt>]
 */
    void setNodeBC(int           code,
                   const string& exp,
                   size_t        dof=1);

/** \brief Assign a given value to components of vector corresponding to sides with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] m Instance of mesh
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] val Value to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned [default: <tt>1</tt>]
 */
    void setSideBC(Mesh&  m,
                   int    code,
                   T_     val,
                   size_t dof);
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setSideBC(Mesh& m,
                   int   code,
                   T_    val);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Remove boundary conditions.
 *  \details This member function copies to current vector a vector
 *  where only non imposed DOF are retained.
 *  @param [in] ms Mesh instance
 *  @param [in] v %Vector (PETScVect instance to copy from)
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned [Default: <tt>0</tt>] or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only one 
 *  degree of freedom
 */
    void removeBC(const Mesh&          ms,
                  const PETScVect<T_>& v,
                  int                  dof=0);

/** \brief Remove boundary conditions.
 *  \details This member function copies to current vector a vector
 *  where only non imposed DOF are retained.
 *  @param [in] ms Mesh instance
 *  @param [in] v %Vector (Vect instance to copy from)
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned [Default: <tt>0</tt>] or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only one 
 *  degree of freedom
 */
    void removeBC(const Mesh&     ms,
                  const Vect<T_>& v,
                  int             dof=0);

/** \brief Remove boundary conditions.
 *  \details This member function copies to current vector a vector
 *  where only non imposed DOF are retained.
 *  @param [in] v Vector (PETScVect instance to copy from)
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned [Default: <tt>0</tt>] or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only one
 *  degree of freedom.
 *  @remark This function is to be used only when the PETScVect instance was constructed by using the Mesh instance
 */
    void removeBC(const PETScVect<T_>& v,
                  int                  dof=0);

/** \brief Remove boundary conditions.
 *  \details This member function copies to current vector a vector
 *  where only non imposed DOF are retained.
 *  @param [in] v Vector (Vect instance to copy from)
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned [Default: <tt>0</tt>] or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only one
 *  degree of freedom.
 *  @remark This function is to be used only when the PETScVect instance was constructed by using the Mesh instance
 */
    void removeBC(const Vect<T_>& v,
                  int             dof=0);

/** \brief Transfer boundary conditions to the vector
 *  @param [in] bc PETScVect instance from which imposed degrees of freedom are copied to current instance
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned [Default: <tt>0</tt>] or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only one 
 *  degree of freedom.
 */
    void transferBC(const PETScVect<T_>& bc,
                    int                  dof=0);

/** \brief Insert boundary conditions
 *  @param [in] m Mesh instance.
 *  @param [in] v PETScVect instance from which free degrees of freedom are copied to current instance.
 *  @param [in] bc PETScVect instance from which imposed degrees of freedom are copied to current instance.
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned [Default: <tt>0</tt>] or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only 
 *  one degree of freedom by node or side
 */
    void insertBC(Mesh&                m, 
                  const PETScVect<T_>& v,
                  const PETScVect<T_>& bc,
                  int                  dof=0);

    void insertBC(Mesh&                m, 
                  const PETScVect<T_>& v,
                  const Vect<T_>&      bc,
                  int                  dof=0);

/** \brief Insert boundary conditions
 *  \details DOF with imposed boundary conditions are set to zero.
 *  @param [in] m Mesh instance.
 *  @param [in] v PETScVect instance from which free degrees of freedom are copied to current instance.
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned [Default: <tt>0</tt>] or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only 
 *  one degree of freedom by node or side
 */
    void insertBC(Mesh&                m,
                  const PETScVect<T_>& v,
                  int                  dof=0);

/** \brief Insert boundary conditions
 *  @param [in] v PETScVect instance from which free degrees of freedom are copied to current instance.
 *  @param [in] bc PETScVect instance from which imposed degrees of freedom are copied to current instance.
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned [Default: <tt>0</tt>] or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only 
 *  one degree of freedom by node or side
 */
    void insertBC(const PETScVect<T_>& v,
                  const PETScVect<T_>& bc,
                  int                  dof=0);

    void insertBC(const PETScVect<T_>& v,
                  const Vect<T_>&      bc,
                  int                  dof=0);

/** \brief Insert boundary conditions
 *  \details DOF with imposed boundary conditions are set to zero.
 *  @param [in] v PETScVect instance from which free degrees of freedom are copied to current instance.
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned [Default: <tt>0</tt>] or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only 
 *  one degree of freedom by node or side
 */
    void insertBC(const PETScVect<T_>& v, 
                  int                  dof=0);

/** \brief Assembly of element vector (as C-array) into Vect instance.
 *  @param [in] el Reference to element instance
 *  @param [in] b Local vector to assemble (C-Array)
 */
    void Assembly(const Element& el,
                  const T_*      b);

/** \brief Assembly of side vector (as C-array) into PETScVect instance.
 *  @param [in] sd Reference to side instance
 *  @param [in] b Local vector to assemble (C-Array)
 */
    void Assembly(const Side& sd,
                  T_*         b);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void DGAssembly(const Element&                          el,
                    const LocalVect<T_,MAX_NB_ELEMENT_DOF>& b);
    void DGAssembly(const Side&                          sd,
                    const LocalVect<T_,MAX_NB_SIDE_DOF>& b);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Evaluate the discrete Gradient vector of the current vector.
 *  \details The resulting gradient is stored in a PETScVect instance
 *  This function handles node vectors assuming P<sub>1</sub> approximation
 *  The gradient is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the gradient, where <tt>v(n,1)</tt>,
 *  <tt>v(n,2)</tt> and <tt>v(n,3)</tt> are respectively the <tt>x</tt> and <tt>y</tt> and 
 *  <tt>z</tt> derivatives at element <tt>n</tt>.
 */
    void getGradient(PETScVect<T_>& v);

/** \brief Evaluate the discrete Gradient vector of the current vector.
 *  \details The resulting gradient is stored in a PETScVect instance
 *  This function handles node vectors assuming P<sub>1</sub> approximation
 *  The gradient is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the gradient, where <tt>v(n,1).x</tt>,
 *  <tt>v(n,2).y</tt> and <tt>v(n,3).z</tt> are respectively the <tt>x</tt> and <tt>y</tt> and 
 *  <tt>z</tt> derivatives at element <tt>n</tt>.
 */
    void getGradient(PETScVect<Point<T_> >& v);

/** \brief Evaluate the discrete curl vector of the current vector.
 *  \details The resulting curl is stored in a PETScVect instance
 *  This function handles node vectors assuming P<sub>1</sub> approximation
 *  The curl is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the curl, where <tt>v(n,1)</tt>,
 *  <tt>v(n,2)</tt> and <tt>v(n,3)</tt> are respectively the <tt>x</tt> and <tt>y</tt> and
 *  <tt>z</tt> <tt>curl</tt> components at element <tt>n</tt>.
 */
    void getCurl(PETScVect<T_>& v);

/** \brief Evaluate the discrete curl vector of the current vector.
 *  \details The resulting curl is stored in a PETScVect instance
 *  This function handles node vectors assuming P<sub>1</sub> approximation
 *  The curl is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the curl, where <tt>v(n,1).x</tt>,
 *  <tt>v(n,2).y</tt> and <tt>v(n,3).z</tt> are respectively the <tt>x</tt> and <tt>y</tt> and 
 *  <tt>z</tt> <tt>curl</tt> components at element <tt>n</tt>.
 */
    void getCurl(PETScVect<Point<T_> >& v);

/** \brief Evaluate the discrete scalar curl in 2-D of the current vector.
 *  \details The resulting curl is stored in a PETScVect instance
 *  This function handles node vectors assuming P<sub>1</sub> approximation
 *  The curl is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the scalar curl.
 */
    void getSCurl(PETScVect<T_>& v);

/** \brief Evaluate the discrete Divergence of the current vector.
 *  \details The resulting divergence is stored in a PETScVect instance
 *  This function handles node vectors assuming P<sub>1</sub> approximation
 *  The divergence is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the divergence.
 */
    void getDivergence(PETScVect<T_>& v);

/** \brief Return average value of vector in a given element
 *  @param [in] el Element instance
 *  @param [in] type Type of element. This is to be chosen
 *  among enumerated values:
 *  <tt>LINE2</tt>, <tt>TRIANG3</tt>, <tt>QUAD4</tt>
 *  <tt>TETRA4</tt>, <tt>HEXA8</tt>
 */
    real_t getAverage(const Element& el,
                            int      type) const;

/** \brief Save vector in a file according to a given format
 *  @param [in] file Output file where to save the vector
 *  @param [in] opt Option to choose file format to save. This is to be chosen
 *  among enumerated values: <tt>GMSH</tt>, <tt>GNUPLOT</tt>, <tt>MATLAB</tt>,
 *  <tt>TECPLOT</tt> and <tt>VTK</tt>
 */
   void save(string file,
             int    opt);
   
/** \brief Multiply by a constant then add to a vector.
 *  @param [in] x PETScVect instance to add
 *  @param [in] a Constant to multiply before adding
 */
    PETScVect<T_> &MultAdd(const PETScVect<T_>& x,
                           const T_&            a);

/** \brief Add to vector the product of a vector by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] x Vect instance by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                   a,
              const PETScVect<T_>& x);

/// \brief Assign a value to an entry for a 1-D vector
/// @param [in] i Rank index in vector (starts at <tt>1</tt>)
/// @param [in] a Value to assign
    void set(size_t i,
             T_     a);

/** \brief Assign a value to an entry for a 2-D vector
 *  @param [in] i First index in vector (starts at <tt>1</tt>)
 *  @param [in] j Second index in vector (starts at <tt>1</tt>)
 *  @param [in] a Value to assign
 */
    void set(size_t i,
             size_t j,
             T_     a);

/** \brief Assign a value to an entry for a 3-D vector
 *  @param [in] i First index in vector (starts at <tt>1</tt>)
 *  @param [in] j Second index in vector (starts at <tt>1</tt>)
 *  @param [in] k Third index in vector (starts at <tt>1</tt>)
 *  @param [in] a Value to assign
 */
    void set(size_t i,
             size_t j,
             size_t k,
             T_     a);

/// \brief Add a value to an entry for a 1-index vector
/// @param [in] i Rank index in vector (starts at <tt>1</tt>)
/// @param [in] a Value to assign
    void add(size_t i,
             T_     a);

/** \brief Add a value to an entry for a 2-index vector
 *  @param [in] i First index in vector (starts at <tt>1</tt>)
 *  @param [in] j Second index in vector (starts at <tt>1</tt>)
 *  @param [in] a Value to assign
 */
    void add(size_t i,
             size_t j,
             T_     a);

/** \brief Assign a value to an entry for a 3-index vector
 *  @param [in] i First index in vector (starts at <tt>1</tt>)
 *  @param [in] j Second index in vector (starts at <tt>1</tt>)
 *  @param [in] k Third index in vector (starts at <tt>1</tt>)
 *  @param [in] a Value to assign
 */
    void add(size_t i,
             size_t j,
             size_t k,
             T_     a);

/// \brief Set all vector entries to zero
    void clear();

/// \brief Operator <tt>[]</tt>
/// @param [in] i Rank index in vector (starts at <tt>0</tt>)
    T_ operator[](size_t i) const;

/** \brief Operator <tt>()</tt>
 *  @param [in] i Rank index in vector (starts at <tt>1</tt>)
 *  <ul>
 *    <li><tt>v(i)</tt> starts at <tt>v(1)</tt> to <tt>v(size())</tt>
 *    <li><tt>v(i)</tt> is the same element as <tt>v[i-1]</tt>
 *  </ul>
 */
    T_ operator()(size_t i) const;

/** \brief Operator () with 2-D indexing (Case of a grid vector)
 *  @param [in] i first index in vector (Number of vector components in the <tt>x</tt>-grid)
 *  @param [in] j second index in vector (Number of vector components in the <tt>y</tt>-grid)
 *  <tt>v(i,j)</tt> starts at <tt>v(1,1)</tt> to <tt>v(getNx(),getNy())</tt>
 */
    T_ operator()(size_t i,
                  size_t j) const;

/** \brief Operator <tt>()</tt> with 3-D indexing (Case of a grid vector)
 *  @param [in] i first index in vector (Number of vector components in the <tt>x</tt>-grid)
 *  @param [in] j second index in vector (Number of vector components in the <tt>y</tt>-grid)
 *  @param [in] k third index in vector (Number of vector components in the <tt>z</tt>-grid)
 *  <tt>v(i,j,k)</tt> starts at <tt>v(1,1,1)</tt> to <tt>v(getNx(),getNy(),getNz())</tt>
 */
    T_ operator()(size_t i,
                  size_t j,
                  size_t k) const;

/// \brief Operator <tt>=</tt> between vectors
    PETScVect<T_> & operator=(const PETScVect<T_>& v);

/** \brief Operator <tt>=</tt>
 *  \details Assign a constant to vector entries
 *  @param [in] a Value to set
 */
    PETScVect<T_> & operator=(const T_& a);

/** \brief Operator <tt>+=</tt>
 *  \details Add vector <tt>x</tt> to current vector instance.
 *  @param [in] v PETScVect instance to add to instance
 */
    PETScVect<T_> & operator+=(const PETScVect<T_>& v);

/** \brief Operator <tt>+=</tt>
 *  \details Add a constant to current vector entries.
 *  @param [in] a Value to add to vector entries
 */
    PETScVect<T_> & operator+=(const T_& a);

/// \brief Operator <tt>-=</tt>
/// @param [in] v Vect instance to subtract from
    PETScVect<T_> & operator-=(const PETScVect<T_>& v);

/** \brief Operator <tt>-=</tt>
 *  \details Subtract constant from vector entries.
 *  @param [in] a Value to subtract from
 */
    PETScVect<T_> & operator-=(const T_& a);

/// \brief Operator <tt>*=</tt>
/// @param [in] a Value to multiply by
    PETScVect<T_> &operator*=(const T_& a);

/// \brief Operator <tt>/=</tt>
/// @param [in] a Value to divide by
    PETScVect<T_> &operator/=(const T_& a);

/// \brief Return reference to Mesh instance
    const Mesh &getMeshPtr() const { return *_theMesh; }
    
/** \brief Return Dot (scalar) product of two vectors
 *  \details A typical use of this operator is <tt>double a = (v,w)</tt>
 *  where <tt>v</tt> and <tt>w</tt> are 2 instances of <tt>PETScVect<double></tt>
 *  @param [in] v PETScVect instance by which the current instance is multiplied
 */
    T_ operator,(const PETScVect<T_>& v) const;

/** \brief Casting operator
 *  \details This member functions enables casting an instance of class PETScVect
 *  into the <tt>Petsc</tt> vector type <tt>Vec</tt>. This is useful when one wants
 *  to usr any Petsc function that is not available in the wrapper (class PETScWrapper)
 *  or PETScVect.
 */
    operator Vec() const { return _v; }

/// \brief Vector assembly.
/// \details This function assembles vector (begins and ends)
    void setAssembly() {
       VecAssemblyBegin(_v);
       VecAssemblyEnd(_v);
    }

/** \brief Insert values into certain locations of the vector.
 *  @param [in] ii Vector containing indices where to insert (Note the indices start from 
 *  0 like any C-array)
 *  @param [in] v Vector of values to insert, corresponding to indices in ii. Here 
 *  the vector has entries of type Point<T_>.
 */
    void Insert(const vector<int>&        ii,
                const vector<Point<T_> >& v);

/** \brief Add values into certain locations of the vector.
 *  @param [in] ii Vector containing indices where to add (Note the indices start from 
 *  0 like any C-array)
 *  @param [in] v Vector of values to add, corresponding to indices in ii
 */
    void Add(const vector<int>& ii,
             const vector<T_>&  v);

 private:
    vector<int>    _ix;
    vector<T_>     _w;
    const vector<string> _var {"x","y","z","t"};
    const vector<string> _var_xit {"x","i","t"};
    const vector<string> _var_ijkt {"i","j","k","t"};
    size_t         _nx, _ny, _nz, _nt, _size, _nb, _dof_type, _nb_dof;
    int            _dg_degree;
    bool           _with_mesh, _with_grid, _with_regex[10], _created;
    Mesh           *_theMesh;
    Grid           *_theGrid;
    string         _name, _regex[10];
    real_t         _time;
    PetscInt       _low, _high;
    Vec            _v;
    T_             _val;
    Fct            _theFct;
    PetscErrorCode _err;

    static Vec Create(size_t   n,
                      MPI_Comm comm=PETSC_COMM_WORLD)
    {
       Vec v;
       VecCreate(comm,&v);
       VecSetSizes(v,PETSC_DECIDE,n);
       VecSetFromOptions(v);
       return v;
    }
    void dof_select(size_t d, vector<size_t> &dof_list);
    int ijk(size_t i, size_t j)           const { return _ny*(i-1)+j-1; }
    int ijk(size_t i, size_t j, size_t k) const { return _ny*_nz*(i-1)+_nz*(j-1)+k-1; }
};


///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_>
PETScVect<T_>::PETScVect()
              : _nx(0), _ny(1), _nz(1), _nt(1), _size(0), _dof_type(NODE_DOF),
                _nb_dof(1), _with_grid(false), _with_mesh(false),
                _created(false), _theMesh(NULL), _name("u"), _time(0.),
                _low(0), _high(0)
{
}


template<class T_>
PETScVect<T_>::PETScVect(size_t n)
              : _nx(n), _ny(1), _nz(1), _size(n), _dof_type(NONE),
                _nb_dof(1), _with_grid(true), _with_mesh(false),
                _created(true), _theMesh(NULL), _name("u"), _time(0.)
{
   _v = Create(_size,PETSC_COMM_WORLD);
   _err = VecGetOwnershipRange(_v,&_low,&_high);
}


template <class T_>
PETScVect<T_>::PETScVect(MPI_Comm comm,
                         size_t   n)
              : _nx(n), _ny(1), _nz(1), _size(n), _dof_type(NONE),
                _nb_dof(1), _with_grid(true), _with_mesh(false),
                _created(true), _theMesh(NULL), _name("u"), _time(0.)
{
   VecCreateMPI(comm,PETSC_DECIDE,_size);
}


template<class T_>
PETScVect<T_>::PETScVect(size_t nx,
                         size_t ny)
              : _nx(nx), _ny(ny), _nz(1), _size(nx*ny), _dof_type(NONE),
                _nb_dof(1), _with_grid(true), _with_mesh(false),
                _created(true), _theMesh(NULL), _name("u"), _time(0.)
{
   _v = Create(_size,PETSC_COMM_WORLD);
   _err = VecGetOwnershipRange(_v,&_low,&_high);
}


template<class T_>
PETScVect<T_>::PETScVect(size_t nx,
                         size_t ny,
                         size_t nz)
              : _nx(nx), _ny(ny), _nz(nz), _size(nx*ny*nz), _dof_type(NONE),
                _nb_dof(1), _with_grid(true), _with_mesh(false),
                _created(true), _theMesh(NULL), _name("u"), _time(0.)
{
   _v = Create(_size,PETSC_COMM_WORLD);
   _err = VecGetOwnershipRange(_v,&_low,&_high);
}


template<class T_>
PETScVect<T_>::PETScVect(size_t nx,
                         size_t ny,
                         size_t nz,
                         size_t nt)
              : _nx(nx), _ny(ny), _nz(nz), _nt(nt), _size(nx*ny*nz*nt), _dof_type(NONE),
                _nb_dof(1), _with_grid(true), _with_mesh(false),
                _created(true), _theMesh(NULL), _name("u"), _time(0.)
{
   _v = Create(_size,PETSC_COMM_WORLD);
   _err = VecGetOwnershipRange(_v,&_low,&_high);
}


template<class T_>
PETScVect<T_>::PETScVect(size_t n,
                         T_*    x)
              : _nx(n), _ny(1), _nz(1), _size(n), _dof_type(NONE),
                _nb_dof(1), _with_grid(true), _with_mesh(false),
                _created(true), _theMesh(NULL), _name("u"), _time(0.)
{
   _v = Create(_size,PETSC_COMM_WORLD);
   _err = VecGetOwnershipRange(_v,&_low,&_high);
   vector<int> ix;
   vector<T_> v;
   for (size_t i=0; i<n; ++i)
      set(i+1,x[i]);
   _err = VecSetValues(_v,ix.size(),&ix[0],&v[0],INSERT_VALUES);
   _err = VecAssemblyBegin(_v);
   _err = VecAssemblyEnd(_v);
}


template<class T_>
PETScVect<T_>::PETScVect(Grid& g)
              : vector<T_>((g.getNx()+1)*(g.getNy()+1)*(g.getNz()+1)),
                _dof_type(NODE_DOF), _nb_dof(1), _dg_degree(-1),
                _with_grid(true), _with_mesh(false), _theMesh(nullptr), _name("#"), _time(0)
{
   setGrid(g);
   for (size_t i=0; i<10; ++i)
      _with_regex[i] = false;
}


template<class T_>
PETScVect<T_>::PETScVect(class Mesh& m,
                         int         nb_dof,
                         int         dof_type)
              : _with_grid(false), _name("u"), _time(0.)
{
   setMesh(m,nb_dof,dof_type);
}


template<class T_>
PETScVect<T_>::PETScVect(class Mesh& m,
                         string      name,
                         real_t      t,
                         int         nb_dof,
                         int         dof_type)
              : _with_grid(false), _name(name), _time(t)
{
   setMesh(m,nb_dof,dof_type);
}


template<class T_>
PETScVect<T_>::PETScVect(const PETScVect<T_>& v,
                         const PETScVect<T_>& bc)
{
   _theMesh = v._theMesh;
   _with_mesh = v._with_mesh;
   _ny = _nb_dof = v._nb_dof;
   _nx = _nb = v._nb;
   _nz = 1;
   _size = _nx*_ny;
   _dof_type = v._dof_type;
   _time = v._time;
   _name = v._name;
   size_t i=1, n=0;
   _created = false;
   setSize(_nx,_ny);
   MESH_ND {
      for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
         set(i,bc[n++]);
         if (The_node.getCode(k) == 0)
            set(i,v(The_node.getDOF(k)));
         i++;
      }
   }
   _with_grid = v._with_grid;
}


template<class T_>
PETScVect<T_>::PETScVect(const PETScVect<T_>& v,
                         size_t               nb_dof,
                         size_t               first_dof)
{
   _time = v._time;
   _name = v._name;
   _theMesh = v._theMesh;
   _with_mesh = v._with_mesh;
   _dof_type = v._dof_type;
   _ny = _nb_dof = nb_dof;
   _nx = _nb = v._nb;
   _nz = 1;
   _created = false;
   setSize(_nx,_ny);
   for (size_t i=1; i<=_nb; i++)
      for (size_t j=1; j<=_nb_dof; j++)
         set(i,j,v(i,j+first_dof-1));
   _with_grid = v._with_grid;
}


template<class T_>
PETScVect<T_>::PETScVect(const PETScVect<T_>& v)
              : _nx(v._nx), _ny(v._ny), _nz(v._nz)
{
   _theMesh = v._theMesh;
   _with_mesh = v._with_mesh;
   _time = v._time;
   _name = v._name;
   _nb_dof = v._nb_dof;
   _dof_type = v._dof_type;
   _nb = v._nb;
   _with_grid = v._with_grid;
   _created = false;
   setSize(_nx,_ny,_nz);
   _err = VecCopy(v._v,_v);
}


template<class T_>
PETScVect<T_>::PETScVect(const PETScVect<T_>& v,
                         size_t               n)
              : _nx(v._nx), _ny(v._ny), _nz(v._nz)
{
   _theMesh = v._theMesh;
   _with_mesh = v._with_mesh;
   _time = v._time;
   _name = v._name;
   _nb_dof = v._nb_dof;
   _with_grid = v._with_grid;
   if (n==1) {
      _ny = _nz = 1;
      _created = false;
      setSize(_nx,_ny,_nz);
      for (size_t i=1; i<=_nx; i++)
         set(i,v(i,1,1));
   }
   else if (n==2) {
      _nx = _nz = 1;
      _created = false;
      setSize(_nx,_ny,_nz);
      for (size_t j=1; j<=_ny; j++)
         set(j,v(1,j,1));
   }
   else if (n==3) {
      _nx = _ny = 1;
      _created = false;
      setSize(_nx,_ny,_nz);
      for (size_t k=1; k<=_nz; k++)
         set(k,v(1,1,k));
   }
}


template <class T_>
PETScVect<T_>::PETScVect(size_t               d,
                         const PETScVect<T_>& v,
                         const string&        name)
{
  if (d<=0)
     throw OFELIException("PETScVect::PETScVect(size_t,PETScVect<T_>,string): Illegal value of nb_dof = "+to_string(d));
   size_t nd=v.getNbDOF();
   vector<size_t> dof_list(nd);
   dof_select(d,dof_list);
   if (_nb_dof>nd)
      throw OFELIException("PETScVect::PETScVect(size_t,PETScVect<T_>,string): Illegal value of dof = "+to_string(nd));
   _time = v._time;
   _theMesh = v._theMesh;
   _with_mesh = v._with_mesh;
   _name = name;
   _nb = v._nb;
   _theMesh = &(v.getMesh());
   _dof_type = v._dof_type;
   _nx = _nb; _ny = _nb_dof; _nz = 1;
   _created = false;
   setSize(_nx,_ny,_nz);
   for (size_t i=1; i<=_nb; i++) {
      for (size_t k=0; k<nd; k++) {
         if (dof_list[k]!=0)
	    set(i,dof_list[k],v(i,k+1));
      }
   }
}


template<class T_>
PETScVect<T_>::~PETScVect()
{
   if (_created)
      _err = VecDestroy(&_v);
}


template<class T_>
void PETScVect<T_>::setMesh(class Mesh& m,
                            int         nb_dof,
                            int         dof_type)
{
   _theMesh = &m;
   _with_mesh = true;
   size_t n = _theMesh->getNbDOF();
   _nb_dof = nb_dof;
   _dof_type = dof_type;
   if (dof_type==NODE_DOF || dof_type==NODE_DOF)
      _nb = _theMesh->getNbNodes();
   else if (dof_type==SIDE_DOF || dof_type==SIDE_DOF)
      _nb = _theMesh->getNbSides();
   else if (dof_type==ELEMENT_DOF || dof_type==ELEMENT_DOF)
      _nb = _theMesh->getNbElements();
   if (nb_dof==0)
      _nb_dof = n/_nb;
   _size = _nb_dof * _nb;
   _nx = _nb, _ny = _nb_dof, _nz = 1;
   _created = false;
   setSize(_nx,_ny,_nz);
}


template<class T_>
void PETScVect<T_>::setGrid(Grid& g)
{
   _theGrid = &g;
   setSize(g.getNx()+1,g.getNy()+1,g.getNz()+1);
   clear();
}

 
template<class T_>
void PETScVect<T_>::set(const T_* v,
                        size_t    n)
{
   setSize(n);
   PetscInt i=1;
   _err = VecSetValues(_v,n,&i,v,INSERT_VALUES);
}


template<class T_>
void PETScVect<T_>::select(const PETScVect<T_>& v,
                           size_t               nb_dof,
                           size_t               first_dof)
{
   _size = nb_dof*v._nb;
   setSize(_size);
   _ix.clear();
   _w.clear();
   size_t i=0;
   for (size_t n=1; n<=v._nb; ++n) {
      for (size_t j=first_dof; j<=nb_dof+first_dof-1; j++) {
         _ix.push_back(i);
         _w.push_back(v(n,j));
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template <class T_>
void PETScVect<T_>::set(const string& exp,
                        size_t        dof)
{
   if (_theMesh==nullptr)
      throw OFELIException("In PETScVect::set(string,dof): No mesh defined");
   Fct& f = _theFct;
   f.set(exp,_var);
   int err;
   vector<real_t> xv = {0.,0.,0.,_time};
   _ix.clear();
   _w.clear();
   if (_dof_type==NODE_DOF)
      throw OFELIException("In PETScVect::set(string,size_t): This member function is for nodewise vectors only.");
   set(The_node.getNbDOF()*(node_label-1)+dof,f(The_node.getCoord(),_time));
   MESH_ND {
      xv[0] = The_node.getCoord(1);
      xv[1] = The_node.getCoord(2);
      xv[2] = The_node.getCoord(3);
      _ix.push_back(The_node.getNbDOF()*(node_label-1)+dof-1);
      _w.push_back(f(xv));
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template <class T_>
void PETScVect<T_>::set(Mesh&         ms,
                        const string& exp,
                        size_t        dof)
{
   setMesh(ms);
   if (_theMesh==nullptr)
      throw OFELIException("In PETScVect::set(ms,string,dof): No mesh is defined");
   if (_dof_type!=NODE_DOF)
      throw OFELIException("In PETScVect::set(string,size_t): This member function is for nodewise vectors only.");
   vector<real_t> xv = {0.,0.,0.,_time};
   _ix.clear();
   _w.clear();
   Fct &f = _theFct[0];
   f.set(exp,_var_xit);
   MESH_ND {
      xv[0] = The_node.getCoord(1);
      xv[1] = The_node.getCoord(2);
      xv[2] = The_node.getCoord(3);
      _ix.push_back(The_node.getNbDOF()*(node_label-1)+dof-1);
      _w.push_back(f(xv));
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template <class T_>
void PETScVect<T_>::set(const PETScVect<real_t>& x,
                        const string&            exp)
{
   setSize(x._nx,x._ny,x._nz,x._nt);
   _ix.clear();
   _w.clear();
   Fct &f = _theFct[0];
   f.set(exp,_var_xit);
   vector<real_t> xv = {0.,0.,_time};
   for (size_t i=0; i<_size; i++) {
      _ix.push_back(i);
      xv[0] = x[i], xv[1] = i+1;
      _w.push_back(f(xv));
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
PetscInt PETScVect<T_>::getLocalSize() const
{
   PetscInt n;
   VecGetLocalSize(_v,&n);
   return n;
}


template<class T_>
void PETScVect<T_>::setSize(size_t nx,
                            size_t ny,
                            size_t nz)
{
   _nx = nx, _ny = ny, _nz = nz;
   _size = _nx*_ny*_nz;
   if (_created)
      _err = VecDestroy(&_v);
   _v = Create(_size,PETSC_COMM_WORLD);
   _err = VecGetOwnershipRange(_v,&_low,&_high);
   _created = true;
}


template<class T_>
PetscScalar PETScVect<T_>::getNorm1() const
{
   PetscScalar s;
   VecNormBegin(_v,NORM_1,&s);
   VecNormEnd(_v,NORM_1,&s);
   return s;
}


template<>
inline PetscScalar PETScVect<real_t>::getNorm2() const
{
   PetscScalar s;
   VecNormBegin(_v,NORM_2,&s);
   VecNormEnd(_v,NORM_2,&s);
   return s;
}


template<>
inline PetscScalar PETScVect<complex_t>::getNorm2() const
{
   PetscScalar s;
   VecNormBegin(_v,NORM_2,&s);
   VecNormEnd(_v,NORM_2,&s);
   return s;
}


template<class T_>
PetscScalar PETScVect<T_>::getNormMax() const
{
   PetscScalar s;
   _err = VecNormBegin(_v,NORM_INFINITY,&s);
   _err = VecNormEnd(_v,NORM_INFINITY,&s);
   return s;
}


template<>
inline real_t PETScVect<real_t>::Norm(NormType t) const
{
   if (t==NORM1)
      return getNorm1();
   else if (t==WNORM1)
      return getWNorm1();
   else if (t==NORM2)
      return getNorm2();
   else if (t==WNORM2)
      return getWNorm2();
   else if (t==NormType(NORM_MAX))
      return getNormMax();
   else
      return 0.;
}


template<class T_>
T_ PETScVect<T_>::getMin() const
{
   T_ s;
   int i;
   _err = VecMin(_v,&i,&s);
   return s;
}


template<class T_>
T_ PETScVect<T_>::getMax() const
{
   T_ s;
   int i;
   _err = VecMax(_v,&i,&s);
   return s;
}


template<class T_>
void PETScVect<T_>::removeBC(const class Mesh&    ms,
                             const PETScVect<T_>& v,
                             int                  dof)
{
   _ix.clear();
   _w.clear();
   if (dof==0) {
      size_t n=0;
      node_loop(&ms) {
         for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
            _ix.push_back(The_node.getDOF(k)-1);
            if (The_node.getCode(k) == 0)
               _w.push_back(v[n]);
            n++;
         }
      }
   }
   else {
      node_loop(&ms) {
         _ix.push_back(The_node.getDOF(dof)-1);
         if (The_node.getCode(dof) == 0)
            _w.push_back(v(node_label));
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::removeBC(const class Mesh& ms,
                             const Vect<T_>&   v,
                             int               dof)
{
   _ix.clear();
   _w.clear();
   if (dof==0) {
      size_t n=0;
      node_loop(&ms) {
         for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
            _ix.push_back(The_node.getDOF(k)-1);
            if (The_node.getCode(k) == 0)
               _w.push_back(v[n]);
            n++;
         }
      }
   }
   else {
      node_loop(&ms) {
         _ix.push_back(The_node.getDOF(dof)-1);
         if (The_node.getCode(dof) == 0)
            _w.push_back(v(node_label));
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::removeBC(const PETScVect<T_>& v,
                             int                  dof)
{
   _ix.clear();
   _w.clear();
   if (dof==0) {
      size_t n=0;
      MESH_ND {
         for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
           _ix.push_back(The_node.getDOF(k)-1);
           if (The_node.getCode(k)==0)
              _w.push_back(v[n]);
            n++;
         }
      }
   }
   else {
      MESH_ND {
         _ix.push_back(The_node.getDOF(dof)-1);
         if (The_node.getCode(dof)==0)
            _w.push_back(v(node_label));
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::removeBC(const Vect<T_>& v,
                             int             dof)
{
   _ix.clear();
   _w.clear();
   if (dof==0) {
      size_t n=0;
      MESH_ND {
         for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
            _ix.push_back(The_node.getDOF(k)-1);
            if (The_node.getCode(k)==0)
               _w.push_back(v[n]);
            n++;
         }
      }
   }
   else {
      MESH_ND {
         _ix.push_back(The_node.getDOF(dof)-1);
         if (The_node.getCode(dof)==0)
            _w.push_back(v(node_label));
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::transferBC(const PETScVect<T_>& bc,
                               int                  dof)
{
   _ix.clear();
   _w.clear();
   size_t i=0;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         MESH_ND {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_node.getCode(k)>0)
                  _w.push_back(bc[i]);
               i++;
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         MESH_ND {
            for (size_t k=1; k<=TheSide.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (TheSide.getCode(k)>0)
                  _w.push_back(bc[i]);
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (_theMesh->NodesAreDOF()) {
         size_t k=0;
         MESH_ND {
            _ix.push_back(i);
            if (The_node.getCode(dof)>0)
               _w.push_back(bc[k]);
            i++;
            k += The_node.getNbDOF();
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         size_t k=0;
         MESH_ND {
            _ix.push_back(i);
            if (TheSide.getCode(dof)>0)
               _w.push_back(bc[k]);
            i++, k += The_node.getNbDOF();
         }
      }
      else
         ;
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::insertBC(class Mesh&          m,
                             const PETScVect<T_>& v,
                             const PETScVect<T_>& bc,
                             int                  dof)
{
   _ix.clear();
   _w.clear();
   size_t i=0;
   if (dof==0) {
      if (m.NodesAreDOF()) {
         node_loop(&m) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_node.getCode(k)==0)
                  _w.push_back(v(The_node.getDOF(k)));
               else
                  _w.push_back(bc[i]);
               i++;
            }
         }
      }
      else if (m.SidesAreDOF()) {
         side_loop(&m) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_side.getCode(k)>=0)
                  _w.push_back(v(The_side.getDOF(k)));
               else
                  _w.push_back(bc[i]);
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (m.NodesAreDOF()) {
         size_t k=dof;
         node_loop(&m) {
            _ix.push_back(i);
            if (The_node.getCode(dof)==0)
               _w.push_back(v(The_node.getDOF(dof)));
            else
               _w.push_back(bc(k));
            i++, k += The_node.getNbDOF();
         }
      }
      else if (m.SidesAreDOF()) {
         size_t k=dof;
         side_loop(&m) {
            _ix.push_back(i);
            if (The_side.getCode(dof) >= 0)
               _w.push_back(v(The_side.getDOF(dof)));
            else
               _w.push_back(bc(k));
            i++, k+=The_side.getNbDOF();
         }
      }
      else
         ;
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::insertBC(class Mesh&          m,
                             const PETScVect<T_>& v,
                             const Vect<T_>&      bc,
                             int                  dof)
{
   _ix.clear();
   _w.clear();
   size_t i=0;
   if (dof==0) {
      if (m.NodesAreDOF()) {
         node_loop(&m) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_node.getCode(k)==0)
                  _w.push_back(v(The_node.getDOF(k)));
               else
                  _w.push_back(bc[i]);
               i++;
            }
         }
      }
      else if (m.SidesAreDOF()) {
         side_loop(&m) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_side.getCode(k) >= 0)
                  _w.push_back(v(The_side.getDOF(k)));
               else
                  _w.push_back(bc[i]);
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (m.NodesAreDOF()) {
         size_t k=dof;
         node_loop(&m) {
            _ix.push_back(i);
            if (The_node.getCode(dof)==0)
               _w.push_back(v(The_node.getDOF(dof)));
            else
               _w.push_back(bc(k));
            i++, k += The_node.getNbDOF();
         }
      }
      else if (m.SidesAreDOF()) {
         size_t k=dof;
         side_loop(&m) {
            _ix.push_back(i);
            if (The_side.getCode(dof)>=0)
               _w.push_back(v(The_side.getDOF(dof)));
            else
               _w.push_back(bc(k));
            i++, k+=The_side.getNbDOF();
         }
      }
      else
         ;
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::insertBC(class Mesh&          m,
                             const PETScVect<T_>& v,
                             int                  dof)
{
   _ix.clear();
   _w.clear();
   size_t i=0;
   if (dof==0) {
      if (m.NodesAreDOF()) {
         node_loop(&m) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_node.getCode(k)==0)
                  _w.push_back(v(The_node.getDOF(k)));
               else
                  _w.push_back(0);
               i++;
            }
         }
      }
      else if (m.SidesAreDOF()) {
         side_loop(&m) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_side.getCode(k)==0)
                  _w.push_back(v(The_side.getDOF(k)));
               else
                  _w.push_back(0);
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (m.NodesAreDOF()) {
         node_loop(&m) {
            _ix.push_back(i);
            if (The_node.getCode(dof)==0)
               _w.push_back(v(The_node.getDOF(dof)));
            i++;
         }
      }
      else if (m.SidesAreDOF()) {
         side_loop(&m) {
            _ix.push_back(i);
            set(i,0);
            if (The_side.getCode(dof)==0)
               _w.push_back(v(The_side.getDOF(dof)));
            else
               _w.push_back(0);
            i++;
         }
      }
      else
         ;
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::insertBC(const PETScVect<T_>& v,
                             const PETScVect<T_>& bc,
                              int                 dof)
{
   _ix.clear();
   _w.clear();
   size_t i=0;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         MESH_ND {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_node.getCode(k)==0)
                  _w.push_back(v(The_node.getDOF(k)));
               else
                  _w.push_back(bc[i]);
               i++;
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         MESH_SD {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_side.getCode(k) >= 0)
                  _w.push_back(v(The_side.getDOF(k)));
               else
                  _w.push_back(bc[i]);
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (_theMesh->NodesAreDOF()) {
         size_t k=dof;
         MESH_ND {
            _ix.push_back(i);
            if (The_node.getCode(dof) == 0)
               _w.push_back(v(The_node.getDOF(dof)));
            else
               _w.push_back(bc(k));
            i++, k+=The_node.getNbDOF();
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         size_t k=dof;
         MESH_SD {
            _ix.push_back(i);
            if (The_side.getCode(dof) >= 0)
               _w.push_back(v(The_side.getDOF(dof)));
            else
               _w.push_back(bc(The_side.getDOF(dof)));
            i++, k+=The_side.getNbDOF();
         }
      }
      else
         ;
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::insertBC(const PETScVect<T_>& v,
                             const Vect<T_>&      bc,
                             int                  dof)
{
   _ix.clear();
   _w.clear();
   size_t i=0;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         MESH_ND {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_node.getCode(k) == 0)
                  _w.push_back(v(The_node.getDOF(k)));
               else
                  _w.push_back(bc[i]);
               i++;
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         MESH_SD {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_side.getCode(k) >= 0)
                  _w.push_back(v(The_side.getDOF(k)));
               else
                  _w.push_back(bc[i]);
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (_theMesh->NodesAreDOF()) {
         size_t k=dof;
         MESH_SD {
            _ix.push_back(i);
            if (The_node.getCode(dof) == 0)
               _w.push_back(v(The_node.getDOF(k)));
            else
               _w.push_back(bc(k));
            i++, k+=The_node.getNbDOF();
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         size_t k=dof;
         MESH_SD {
            _ix.push_back(i);
            if (The_side.getCode(dof) >= 0)
               _w.push_back(v(The_side.getDOF(dof)));
            else
               _w.push_back(bc(The_side.getDOF(dof)));
            i++, k+=The_side.getNbDOF();
         }
      }
      else
         ;
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::setNodeBC(class Mesh&   m,
                              int           code,
                              const string& exp,
                              size_t        dof)
{
   _ix.clear();
   _w.clear();
   int err;
   vector<real_t> xv = {0.,0.,0.,_time};
   Fct &f = _theFct[0];
   f.set(exp,_var_xit);
   size_t k=0;
   node_loop(&m) {
      xv[0] = The_node.getCoord(1);
      xv[1] = The_node.getCoord(2);
      xv[2] = The_node.getCoord(3);
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
         if (The_node.getCode(dof)==code) {
            _ix.push_back(k);
            _w.push_back(f(xv));
         }
         k++;
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::setNodeBC(class Mesh& m,
                              int         code,
                              T_          val,
                              size_t      dof)
{
   _ix.clear();
   _w.clear();
   size_t k=0;
   if (m.getDOFSupport()==NODE_DOF) {
      node_loop(&m) {
         for (size_t i=1; i<=the_node->getNbDOF(); i++) {
            if (The_node.getCode(dof)==code) {
               _ix.push_back(k);
               _w.push_back(val);
            }
            k++;
         }
      }
   }
   if (m.getDOFSupport()==SIDE_DOF) {
      boundary_side_loop(&m) {
         for (size_t i=1; i<=theSide->getNbDOF(); i++) {
            if (theSide->getCode(dof)==code) {
               _ix.push_back(k);
               _w.push_back(val);
            }
            k++;
         }
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::setNodeBC(int           code,
                              const string& exp,
                              size_t        dof)
{
   _ix.clear();
   _w.clear();
   vector<real_t> xv = {0.,0.,0.,_time};
   Fct &f = _theFct[0];
   f.set(exp,_var_xit);
   size_t k=0;
   MESH_ND {
      xv[0] = The_node.getCoord(1);
      xv[1] = The_node.getCoord(2);
      xv[2] = The_node.getCoord(3);
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
         if (The_node.getCode(dof)==code) {
            _ix.push_back(k);
            _w.push_back(f(xv));
         }
         k++;
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::setNodeBC(int    code,
                              T_     val,
                              size_t dof)
{
   _ix.clear();
   _w.clear();
   size_t k=0;
   if (_theMesh->getDOFSupport()==NODE_DOF) {
      MESH_ND {
         for (size_t i=1; i<=The_node.getNbDOF(); i++) {
            if (The_node.getCode(dof)==code) {
               _ix.push_back(k);
               _w.push_back(val);
            }
            k++;
         }
      }
   }
   if (_theMesh->getDOFSupport()==SIDE_DOF) {
      MESH_BD_SD {
         for (size_t i=1; i<=TheSide.getNbDOF(); i++) {
            if (TheSide.getCode(dof)==code) {
               _ix.push_back(k);
               _w.push_back(val);
            }
            k++;
         }
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::insertBC(const PETScVect<T_>& v,
                             int                  dof)
{
   _ix.clear();
   _w.clear();
   size_t i=0;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         MESH_ND {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_node.getCode(k)==0)
                  _w.push_back(v(The_node.getDOF(k)));
               else
                  _w.push_back(0);
               i++;
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         MESH_SD {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               _ix.push_back(i);
               if (The_side.getCode(k) == 0)
                  _w.push_back(v(The_node.getDOF(k)));
               else
                  _w.push_back(0);
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (_theMesh->NodesAreDOF()) {
         MESH_ND {
            _ix.push_back(i);
            if (The_node.getCode(dof) == 0)
               _w.push_back(v(The_node.getDOF(dof)));
            else
               _w.push_back(0);
            i++;
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         MESH_SD {
            _ix.push_back(i);
            if (The_side.getCode(dof) == 0)
               _w.push_back(v(The_side.getDOF(dof)));
            else
               _w.push_back(0);
            i++;
         }
      }
      else
         ;
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template<class T_>
void PETScVect<T_>::DGAssembly(const Element&                          el,
                               const LocalVect<T_,MAX_NB_ELEMENT_DOF>& b)
{
   _ix.clear();
   _w.clear();
   for (size_t i=1; i<=el.getNbDOF(); ++i) {
      if (el.getDOF(i)!=0) {
         _ix.push_back(el.getDOF(i)-1);
         _w.push_back(b(i));
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],ADD_VALUES);
}


template<class T_>
void PETScVect<T_>::DGAssembly(const Side&                          sd,
                               const LocalVect<T_,MAX_NB_SIDE_DOF>& b)
{
   _ix.clear();
   _w.clear();
   for (size_t i=1; i<=sd.getNbDOF(); ++i) {
      if (sd.getDOF(i)!=0) {
         _ix.push_back(sd.getDOF(i)-1);
         _w.push_back(b(i));
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],ADD_VALUES);
}


template<class T_>
void PETScVect<T_>::Assembly(const Element& el,
                             const T_*      b)
{
   _ix.clear();
   _w.clear();
   size_t i=0;
   for (size_t n=1; n<=el.getNbNodes(); ++n) {
      Node *nd=el(n);
      for (size_t k=1; k<=nd->getNbDOF(); ++k) {
         _ix.push_back(nd->getDOF(k)-1);
         _w.push_back(b[i++]);
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],ADD_VALUES);
}


template<class T_>
void PETScVect<T_>::Assembly(const Side& sd,
                             T_*         b)
{
   _ix.clear();
   _w.clear();
   size_t i=0;
   for (size_t n=1; n<=sd.getNbNodes(); ++n) {
      Node *nd=sd(n);
      for (size_t k=1; k<=nd->getNbDOF(); ++k) {
         _ix.push_back(nd->getDOF(k)-1);
         _w.push_back(b[i++]);
      }
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],ADD_VALUES);
}


template <class T_>
void PETScVect<T_>::getGradient(class PETScVect<T_>& v)
{
   T_ a;
   real_t b;
   Point<T_> aa;
   v.setMesh(*_theMesh,_theMesh->getDim(),ELEMENT_DOF);
   v.setTime(_time);
   _ix.clear();
   _w.clear();
   size_t k=0;
   MESH_EL {
      if (The_element.getShape()==LINE) {
         a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
         b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
         _ix.push_back(element_label-1);
         _w.push_back(a/b);
      }
      else if (The_element.getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         aa = (*this)(The_element(1)->n())*dsh[0] + 
              (*this)(The_element(2)->n())*dsh[1] + 
              (*this)(The_element(3)->n())*dsh[2];
         _ix.push_back(k++); _ix.push_back(k++);
         _w.push_back(aa.x); _w.push_back(aa.y);
      }
      else if (The_element.getShape()==TETRAHEDRON) {
         Tetra4 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         aa = (*this)(The_element(1)->n())*dsh[0] + 
              (*this)(The_element(2)->n())*dsh[1] + 
              (*this)(The_element(3)->n())*dsh[2] + 
              (*this)(The_element(4)->n())*dsh[3];
         v.set(element_label,1,aa.x);
         v.set(element_label,2,aa.y);
         v.set(element_label,3,aa.z);
         _ix.push_back(k++); _ix.push_back(k++); _ix.push_back(k++);
         _w.push_back(aa.x); _w.push_back(aa.y); _w.push_back(aa.z);
      }
      else
         throw OFELIException("PETScVect::getGradient(): This function doesn't work for this element.");
   }
   VecSetValues(_v,_ix.size(),&_ix[0],&_w[0],INSERT_VALUES);
   setAssembly();
}


template <class T_>
void PETScVect<T_>::getGradient(PETScVect<Point<T_>>& v)
{
   T_ a;
   real_t b;
   _ix.clear();
   vector<Point<T_> > w;
   Point<T_> aa;
   v.setMesh(*_theMesh,_theMesh->getDim());
   v.setTime(_time);
   MESH_EL {
      if (The_element.getShape()==LINE) {
         a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
         b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
         _ix.push_back(element_label-1);
         w.push_back(Point<T_>(a/b));
      }
      else if (the_element->getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t>> dsh = t.DSh();
         aa = (*this)(The_element(1)->n())*dsh[0] + 
              (*this)(The_element(2)->n())*dsh[1] + 
              (*this)(The_element(3)->n())*dsh[2];
         _ix.push_back(element_label-1);
         w.push_back(aa);
      }
      else if (_theMesh->getShape()==TETRAHEDRON) {
         Tetra4 t(the_element);
         vector<Point<real_t> > dsh = t.DSh();
         aa = (*this)(The_element(1)->n())*dsh[0] +
              (*this)(The_element(2)->n())*dsh[1] +
              (*this)(The_element(3)->n())*dsh[2] +
              (*this)(The_element(4)->n())*dsh[3];
         _ix.push_back(element_label-1);
         w.push_back(aa);
      }
      else
         throw OFELIException("PETScVect::getGradient(): This function doesn't work for this element.");
   }
   Insert(_ix,w);
}


template <class T_>
void PETScVect<T_>::getCurl(PETScVect<T_>& v)
{
   T_ a;
   real_t b;
   Point<T_> du, dv, dw;
   v.setMesh(*_theMesh,_theMesh->getDim(),ELEMENT_DOF);
   v.setTime(_time);
   MESH_EL {
      if (The_element.getShape()==LINE) {
         a = (*this)(The_element(2)) - (*this)(The_element(1));
         b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
         v.set(element_label,a/b);
      }
      else if (The_element.getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t>> dsh = t.DSh();
         du = (*this)(The_element(1)->n())*dsh[0] + 
              (*this)(The_element(2)->n())*dsh[1] + 
              (*this)(The_element(3)->n())*dsh[2];
         v.set(element_label,1, du.y);
         v.set(element_label,2,-du.x);
      }
      else if (The_element.getShape()==TETRAHEDRON) {
         Tetra4 t(the_element);
         vector<Point<real_t>> dsh = t.DSh();
         du = (*this)(The_element(1)->n(),1)*dsh[0] + 
              (*this)(The_element(2)->n(),1)*dsh[1] + 
              (*this)(The_element(3)->n(),1)*dsh[2] + 
              (*this)(The_element(4)->n(),1)*dsh[3];
         dv = (*this)(The_element(1)->n(),2)*dsh[0] + 
              (*this)(The_element(2)->n(),2)*dsh[1] + 
              (*this)(The_element(3)->n(),2)*dsh[2] + 
              (*this)(The_element(4)->n(),2)*dsh[3];
         dw = (*this)(The_element(1)->n(),3)*dsh[0] + 
              (*this)(The_element(2)->n(),3)*dsh[1] + 
              (*this)(The_element(3)->n(),3)*dsh[2] + 
              (*this)(The_element(4)->n(),3)*dsh[3];
         v.set(element_label,1,dw.y - dv.z);
         v.set(element_label,2,du.z - dw.x);
         v.set(element_label,3,dv.x - du.y);
      }
      else
         throw OFELIException("PETScVect::getCurl(): This function doesn't work for this element.");
   }
}


template <class T_>
void PETScVect<T_>::getCurl(PETScVect<Point<T_> >& v)
{
   T_ a;
   real_t b;
   Point<T_> du, dv, dw;
   v.setMesh(*_theMesh,_theMesh->getDim(),ELEMENT_DOF);
   v.setTime(_time);
   MESH_EL {
      if (The_element.getShape()==LINE) {
         a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
         b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
         v.set(element_label,a/b);
      }
      else if (The_element.getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t>> dsh = t.DSh();
         du = (*this)(The_element(1)->n())*dsh[0] + 
              (*this)(The_element(2)->n())*dsh[1] + 
              (*this)(The_element(3)->n())*dsh[2];
         v.set(element_label,1, du.y);
         v.set(element_label,2,-du.x);
      }
      else if (The_element.getShape()==TETRAHEDRON) {
         Tetra4 t(the_element);
         vector<Point<real_t>> dsh = t.DSh();
         du = (*this)(The_element(1)->n(),1)*dsh[0] + 
              (*this)(The_element(2)->n(),1)*dsh[1] + 
              (*this)(The_element(3)->n(),1)*dsh[2] + 
              (*this)(The_element(4)->n(),1)*dsh[3];
         dv = (*this)(The_element(1)->n(),2)*dsh[0] + 
              (*this)(The_element(2)->n(),2)*dsh[1] + 
              (*this)(The_element(3)->n(),2)*dsh[2] + 
              (*this)(The_element(4)->n(),2)*dsh[3];
         dw = (*this)(The_element(1)->n(),3)*dsh[0] + 
              (*this)(The_element(2)->n(),3)*dsh[1] + 
              (*this)(The_element(3)->n(),3)*dsh[2] + 
              (*this)(The_element(4)->n(),3)*dsh[3];
         v.set(element_label,Point<T_>(dw.y-dv.z,du.z-dw.x,dv.x-du.y));
      }
      else
         throw OFELIException("PETScVect::getCurl(): This function doesn't work for this element.");
   }
}


template <class T_>
void PETScVect<T_>::getSCurl(PETScVect<T_>& v)
{
   if (_theMesh->getDim()==1 || _theMesh->getDim()==3)
      throw OFELIException("PETScVect::getSCurl(): This function is valid for 2-D only.");
   Point<T_> du, dv;
   v.setMesh(*_theMesh,1);
   v.setTime(_time);
   MESH_EL {
      if (The_element.getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t>> dsh = t.DSh();
         du = (*this)(The_element(1)->n(),1)*dsh[0] + 
              (*this)(The_element(2)->n(),1)*dsh[1] + 
              (*this)(The_element(3)->n(),1)*dsh[2];
         dv = (*this)(The_element(1)->n(),2)*dsh[0] + 
              (*this)(The_element(2)->n(),2)*dsh[1] + 
              (*this)(The_element(3)->n(),2)*dsh[2];
         v.set(element_label,dv.x - du.y);
      }
      else
         throw OFELIException("PETScVect::getSCurl(): This function doesn't work for this element");
   }
}


template <class T_>
void PETScVect<T_>::getDivergence(PETScVect<T_>& v)
{
   T_ a;
   real_t b;
   Point<T_> du, dv, dw;
   v.setMesh(*_theMesh,1);
   v.setTime(_time);
   MESH_EL {
      if (The_element.getShape()==LINE) {
         a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
         b = The_element(2)->getX() - The_element(1)->getX();
         v.set(element_label,a/b);
      }
      else if (The_element.getShape()==TRIANGLE) {
         Triang3 t(the_element);
         vector<Point<real_t>> dsh = t.DSh();
         du = (*this)(The_element(1)->n(),1)*dsh[0] + 
              (*this)(The_element(2)->n(),1)*dsh[1]+ 
              (*this)(The_element(3)->n(),1)*dsh[2];
         dv = (*this)(The_element(1)->n(),2)*dsh[0] + 
              (*this)(The_element(2)->n(),2)*dsh[1] + 
              (*this)(The_element(3)->n(),2)*dsh[2];
         v.set(element_label,du.x+dv.y);
      }
      else if (The_element.getShape()==TETRAHEDRON) {
         Tetra4 t(the_element);
         vector<Point<real_t>> dsh = t.DSh();
         du = (*this)(The_element(1)->n(),1)*dsh[0] + 
              (*this)(The_element(2)->n(),1)*dsh[1] + 
              (*this)(The_element(3)->n(),1)*dsh[2] + 
              (*this)(The_element(4)->n(),1)*dsh[3];
         dv = (*this)(The_element(1)->n(),2)*dsh[0] + 
              (*this)(The_element(2)->n(),2)*dsh[1] + 
              (*this)(The_element(3)->n(),2)*dsh[2] + 
              (*this)(The_element(4)->n(),2)*dsh[3];
         dw = (*this)(The_element(1)->n(),3)*dsh[0] + 
              (*this)(The_element(2)->n(),3)*dsh[1] + 
              (*this)(The_element(3)->n(),3)*dsh[2] + 
              (*this)(The_element(4)->n(),3)*dsh[3];
         v.set(element_label,du.x+dv.y+dw.z);
      }
      else
         throw OFELIException("PETScVect::getDivergence(): This function doesn't work for this element.");
   }
}


template<class T_>
real_t PETScVect<T_>::getAverage(const Element& el,
                                 int            type) const
{
   switch (type) {

      case LINE2:
         return 0.5*((*this)(el(1)->n())+(*this)(el(2)->n()));

      case TRIANG3: 
         return OFELI_THIRD*((*this)(el(1)->n()) +
                             (*this)(el(2)->n()) +
                             (*this)(el(3)->n()));

      case QUAD4:
         return 0.25*((*this)(el(1)->n()) +
                      (*this)(el(2)->n()) +
                      (*this)(el(3)->n()) +
                      (*this)(el(4)->n()));

      case TETRA4:
         return 0.25*((*this)(el(1)->n()) +
                      (*this)(el(2)->n()) +
                      (*this)(el(3)->n()) +
                      (*this)(el(4)->n()));

      case HEXA8:
         return 0.125*((*this)(el(1)->n()) +
                       (*this)(el(2)->n()) +
                       (*this)(el(3)->n()) +
                       (*this)(el(4)->n()) +
                       (*this)(el(5)->n()) +
                       (*this)(el(6)->n()) +
                       (*this)(el(7)->n()) +
                       (*this)(el(8)->n()));
   }
   return 0.;
}


template<class T_>
void PETScVect<T_>::save(string file,
                         int    opt)
{
   void saveField(PETScVect<real_t>& v, string output_file, int opt);
   saveField(*this,file,opt);
}


template<class T_>
PETScVect<T_> &PETScVect<T_>::MultAdd(const PETScVect<T_>& x,
                                      const T_&            a)
{
   for (size_t i=1; i<=_size; ++i)
      add(i,a*x(i));
   return *this;
}


template<class T_>
void PETScVect<T_>::Axpy(T_                   a,
                         const PETScVect<T_>& x)
{
   for (size_t i=1; i<=_size; i++)
      add(i,a*x(i));
}


template<class T_>
void PETScVect<T_>::set(size_t i,
                        T_     a)
{
   _err = VecSetValue(_v,i-1,a,INSERT_VALUES);
}


template<class T_>
void PETScVect<T_>::set(size_t i,
                        size_t j,
                        T_     a)
{
   _err = VecSetValue(_v,ijk(i,j),a,INSERT_VALUES);
}


template<class T_>
void PETScVect<T_>::set(size_t i,
                        size_t j,
                        size_t k,
                        T_     a)
{
   _err = VecSetValue(_v,ijk(i,j,k),a,INSERT_VALUES);
}


template<class T_>
void PETScVect<T_>::add(size_t i,
                        T_     a)
{
   _err = VecSetValue(_v,i-1,a,ADD_VALUES);
}


template<class T_>
void PETScVect<T_>::add(size_t i,
                        size_t j,
                        T_     a)
{
   _err = VecSetValue(_v,ijk(i,j),a,ADD_VALUES);
}


template<class T_>
void PETScVect<T_>::add(size_t i,
                        size_t j,
                        size_t k,
                        T_     a)
{
   _err = VecSetValue(_v,ijk(i,j,k),a,ADD_VALUES);
}


template<class T_>
T_ PETScVect<T_>::operator[](size_t i) const
{
   static T_ a;
   int ii=i;
   VecGetValues(_v,1,&ii,&a);
   return a;
}


template<class T_>
T_ PETScVect<T_>::operator()(size_t i) const
{
   static T_ a;
   int ii=i-1;
   VecGetValues(_v,1,&ii,&a);
   return a;
}


template<class T_>
T_ PETScVect<T_>::operator()(size_t i,
                             size_t j) const
{
   return (*this)[ijk(i,j)];
}


template<class T_>
T_ PETScVect<T_>::operator()(size_t i,
                             size_t j,
                             size_t k) const
{
   return (*this)[ijk(i,j,k)];
}


template<class T_>
PETScVect<T_> &PETScVect<T_>::operator=(const PETScVect<T_>& v)
{
   _theMesh = v._theMesh;
   _with_mesh = v._with_mesh;
   _time = v._time;
   _name = v._name;
   _nb_dof = v._nb_dof;
   _dof_type = v._dof_type;
   _nb = v._nb;
   _with_grid = v._with_grid;
   _created = false;
   setSize(v._nx,v._ny,v._nz);
   _err = VecCopy(v._v,_v);
   return *this;
}


template<class T_>
void PETScVect<T_>::clear()
{
   _err = VecZeroEntries(_v);
}


template<class T_>
PETScVect<T_> &PETScVect<T_>::operator=(const T_& a)
{
   _err = VecSet(_v,a);
   return *this;
}


template<class T_>
PETScVect<T_> &PETScVect<T_>::operator+=(const PETScVect<T_>& v)
{
   _err = VecAXPY(_v,1.0,v._v);
   return *this;
}


template<class T_>
PETScVect<T_> &PETScVect<T_>::operator+=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      add(i,a);
   return *this;
}


template<class T_>
PETScVect<T_> &PETScVect<T_>::operator-=(const PETScVect<T_>& v)
{
   _err = VecAXPY(_v,-1.0,v._v);
   return *this;
}


template<class T_>
PETScVect<T_> &PETScVect<T_>::operator-=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      add(i,-a);
   return *this;
}


template<class T_>
PETScVect<T_> &PETScVect<T_>::operator*=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      set(i,(*this)(i)*a);
   return *this;
}


template<class T_>
PETScVect<T_> &PETScVect<T_>::operator/=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      set(i,(*this)(i)/a);
   return *this;
}


template<class T_>
T_ PETScVect<T_>::operator,(const PETScVect<T_>& v) const
{
   T_ p = 0;
   for (size_t i=0; i<_size; i++)
      p += (*this)[i] * v[i];
   return p;
}


template <class T_>
void PETScVect<T_>::dof_select(size_t          d,
                               vector<size_t>& dof_list)
{
   size_t j=d, nd=dof_list.size();
   for (size_t k=0; k<nd; k++) {
      size_t kk=size_t(pow(10.,real_t(nd-k-1)));
      size_t m=j/kk;
      dof_list[k] = m;
      j -= m*kk;
   }
   _nb_dof = 0;
   for (size_t k=0; k<nd; k++)
      if (dof_list[k]!=0)
         _nb_dof++;
}


template <class T_>
void PETScVect<T_>::Insert(const vector<int>&         ii,
                           const vector<Point<T_> >&  v)
{
   _err = VecSetValues(_v,ii.size(),&ii[0],&v[0],INSERT_VALUES);
   _err = VecAssemblyBegin(_v);
   _err = VecAssemblyEnd(_v);
}


template <class T_>
void PETScVect<T_>::Add(const vector<int>& ii,
                        const vector<T_>&  v)
{
   _err = VecSetValues(_v,ii.size(),&ii[0],&v[0],ADD_VALUES);
   _err = VecAssemblyBegin(_v);
   _err = VecAssemblyEnd(_v);
}


///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////

/** \fn PETScVect<T_> operator+(const PETScVect<T_> &x, const PETScVect<T_> &y)
 *  \brief Operator + (Addition of two instances of class PETScVect)
 *  \ingroup VectMat
 *  \return  <tt>x + y</tt>
 */
template<class T_>
PETScVect<T_> operator+(const PETScVect<T_>& x,
                        const PETScVect<T_>& y)
{
   PETScVect<T_> v(x);
   for (size_t i=0; i<x.size(); ++i)
      v.add(i+1,y[i]);
   return v;
}


/** \fn PETScVect<T_> operator-(const PETScVect<T_>& x, const PETScVect<T_>& y)
 *  \brief Operator - (Difference between two instances of class PETScVect)
 *  \ingroup VectMat
 *  \return <tt>x - y</tt>
 */
template<class T_>
PETScVect<T_> operator-(const PETScVect<T_>& x,
                        const PETScVect<T_>& y)
{
   PETScVect<T_> v(x.size());
   for (size_t i=0; i<x.size(); ++i)
      v.set(i+1,x[i]-y[i]);
   return v;
}


/** \fn PETScVect<T_> operator*(const T_& a, const PETScVect<T_>& x)
 *  \brief Operator * (Premultiplication of vector by constant)
 *  \ingroup VectMat
 *  \return <tt>a*x</tt>
 */
template<class T_>
PETScVect<T_> operator*(const T_&            a,
                        const PETScVect<T_>& x)
{
   PETScVect<T_> v(x.size());
   for (size_t i=1; i<=x.size(); ++i)
      v.set(i,a*x(i));
   return v;
}


/** \fn PETScVect<T_> operator*(const PETScVect<T_>& x, const T_& a)
 *  \brief Operator * (Postmultiplication of vector by constant)
 *  \ingroup VectMat
 *  \return <tt>x*a</tt>
 */
template<class T_>
PETScVect<T_> operator*(const PETScVect<T_>& x,
                        const T_&            a)
{
   PETScVect<T_> v(x.size());
   for (size_t i=1; i<=x.size(); ++i)
      v.set(i,a*x(i));
   return v;
}


/** \fn PETScVect<T_> operator/(const T_ &a, const PETScVect<T_> &x)
 *  \brief Operator / (Divide vector entries by constant)
 *  \ingroup VectMat
 *  \return <tt>x/a</tt>
 */
template<class T_>
PETScVect<T_> operator/(const PETScVect<T_>& x,
                        const T_&            a)
{
   PETScVect<T_> v(x);
   for (size_t i=1; i<=x.size(); ++i)
      v.set(i,x(i)/a);
   v.setAssembly();
   return v;
}


/** \fn double Dot(const PETScVect<double> &x, const PETScVect<double> &y)
 *  \brief Calculate dot product of two vectors
 *  \ingroup VectMat
 *  \return Dot (inner or scalar) product
 *  Calculate dot (scalar) product of two vectors
 */
template<class T_>
PetscScalar Dot(const PETScVect<T_>& x,
                const PETScVect<T_>& y)
{
   T_ d;
   VecDotBegin(x,y,&d);
   VecDotEnd(x,y,&d);
   return d;
}


/** \fn real_t Discrepancy(PETScVect<T_> &x, const PETScVect<T_> &y, int n, int type)
 *  \ingroup Util
 *  \brief Return discrepancy between 2 vectors <tt>x</tt> and <tt>y</tt>
 *  @param [in,out] x First vector (Instance of class Vect). On output, <tt>x</tt> is 
 *  assigned the vector <tt>y</tt>
 *  @param [in] y Second vector (Instance of class Vect)
 *  @param [in] n Type of norm
 *     - 1: Weighted 1-Norm
 *     - 2: Weighted 2-Norm
 *     - 0: Max-Norm
 *  @param [in] type Discrepancy type (0: Absolute, 1: Relative [Default])
 *  @return Computed discrepancy value
 */
template<class T_>
inline PetscScalar Discrepancy(      PETScVect<T_>& x,
                               const PETScVect<T_>& y,
                                     int            n,
                                     int            type=1)
{
   size_t s=x.size();
   real_t old=0., d=0.;
   if (n==0) {
      old = x.getNormMax();
      for (size_t i=0; i<s; i++) 
         if (d<std::abs(x[i]-y[i]))
            d = std::abs(x[i]-y[i]);
   }
   else if (n==1) {
      old = x.getWNorm1();
      for (size_t i=0; i<s; i++)
         d += std::abs(x[i]-y[i]);
      d /= s;
   }
   else if (n==2) {
      old = x.getWNorm2();
      for (size_t i=0; i<s; i++)
         d += Abs2(x[i]-y[i]);
      d = sqrt(d/s);
   }
   else
      ;
   x = y;
   if (type==1 && old>0.)
      d /= old;
   return d;
}


/** \fn ostream &operator<<(ostream &s, PETScVect<T_> &v)
 *  \brief Output vector in output stream
 *  \ingroup VectMat
 */
template<class T_>
ostream &operator<<(ostream&       s,
                    PETScVect<T_>& v)
{
   s << "Vector Size: " << v.size() << std::endl;
   if (v.size() <= 0)
      return s;
   v.setAssembly();
   VecView(v,PETSC_VIEWER_STDOUT_WORLD);
   s << endl;
   return s;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif

#endif
