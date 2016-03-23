/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

                 Definition of Template class Vect for vectors

  ==============================================================================*/


#ifndef __VECT_H
#define __VECT_H

#if !defined (USE_EIGEN)
#include <vector>
using std::vector;
#endif

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
using std::ostream;
using std::istream;
using std::endl;

#include <iomanip>
using std::setw;
using std::ios;
using std::setprecision;

#include "OFELI_Config.h"
#include "util/macros.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Point.h"
#include "io/fparser/fparser.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Tetra4.h"

#if defined (USE_EIGEN)
#include <Eigen/Dense>
using Eigen::Matrix;
#endif

extern FunctionParser theParser;

/** \defgroup VectMat Vector and Matrix
 *  \brief Gathers vector and matrix related classes
 */

/*! \file Vect.h
 *  \brief Definition file for class Vect.
 */

/*!
 *  \addtogroup OFELI
 *  @{
 */

//! A namespace to group all library classes, functions, ...
namespace OFELI {

/*! 
 *  \class Vect
 *  \ingroup VectMat
 *  \brief To handle general purpose vectors.
 *
 * This template class enables defining and manipulating vectors of various data types.
 * It inherits from the class std::vector
 * An instance of class Vect can be:
 * <ul>
 *   <li> A simple vector of given size
 *   <li> A vector with up to three indices, <i>i.e.</i>, an entry of the vector
 *        can be <tt>a(i)</tt>, <tt>a(i,j)</tt> or <tt>a(i,j,k)</tt>. This feature
 *        is useful, for instance, in the case of a structured grid
 *   <li> A vector associate to a finite element mesh. In this case, a constructor
 *        uses a reference to the Mesh instance. The size of the vector is by
 *        default equal to the number of nodes <tt>x</tt> the number of degrees of freedom
 *        by node. If the degrees of freedom are supported by elements or sides,
 *        then the vector is sized accordingly
 * </ul>
 * Operators \b =, \b [] and \b () are overloaded so that one can write for instance:
 *
 * \verbatim
    Vect<real_t> u(10), v(10);
    v = -1.0;
    u = v;
    u(3) = -2.0;
   \endverbatim
 *
 * to set vector \b v entries to \b -1, copy vector \b v into vector \b u and assign
 * third entry of \b v to \b -2. Note that entries of \b v are here \b v(1), \b v(2),
 * ..., \b v(10), \e i.e. vector entries start at index \b 1.
 *
 * \tparam T_ Data type (real_t, float, complex<real_t>, ...)
 */

string itos(int i);
string dtos(real_t d);
class Mesh;
class Element;
class Side;

#if defined (USE_PETSC)
template<class T_> class PETScVect;
#endif

#if defined (USE_EIGEN)
template<class T_>
class Vect
#else
template<class T_>
class Vect
          : public vector<T_>
#endif
/// @endcond
{

 public:

#if defined (USE_EIGEN)
/*! \typedef VectorX
 *  \ingroup VectMat
 *  \brief This type is the vector type in the %Eigen library
 *  @remark: This type is available only if the %Eigen library was installed in conjunction
 *  with OFELI
 */
    typedef Eigen::Matrix<T_,Eigen::Dynamic,1> VectorX;
#endif

#if !defined (USE_EIGEN)
    using vector<T_>::size;
#endif

/// \brief Default Constructor.
/// Initialize a zero-length vector
    Vect();

/// \brief Constructor setting vector size.
/// @param [in] n Size of vector
    Vect(size_t n);

/** \brief Constructor of a 2-D index vector.
 *  \details This constructor can be used for instance for a 2-D grid vector
 *  @param [in] nx Size for the first index
 *  @param [in] ny Size for the second index
 *  @remark The size of resulting vector is nx*ny
 */
    Vect(size_t nx,
         size_t ny);

/** \brief Constructor of a 3-D index vector.
 *  \details This constructor can be used for instance for a 3-D grid vector
 *  @param [in] nx Size for the first index
 *  @param [in] ny Size for the second index
 *  @param [in] nz Size for the third index
 *  @remark The size of resulting vector is nx*ny*nz
 */
    Vect(size_t nx,
         size_t ny,
         size_t nz);

/** \brief Create an instance of class Vect as an image of a C/C++ array.
 *  @param [in] n Dimension of vector to construct
 *  @param [in] x C-array to copy
 */
    Vect(size_t n,
         T_*    x);

/** \brief Constructor with a mesh instance
 *  @param [in] m Mesh instance
 *  @param [in] nb_dof Number of degrees of freedom per node, element or side
 *  If <tt>nb_dof</tt> is set to <tt>0</tt> (default value) the constructor picks this number from
 *  the Mesh instance
 *  @param [in] dof_type Type of degrees of freedom. To be given among the enumerated
 *  values: <tt>NODE_FIELD</tt>, <tt>ELEMENT_FIELD</tt>, <tt>SIDE_FIELD</tt> or <tt>EDGE_FIELD</tt>
 *  (Default: <tt>NODE_FIELD</tt>)
 */
    Vect(Mesh& m,
         int   nb_dof=0,
         int   dof_type=NODE_FIELD);

/** \brief Constructor with a mesh instance giving name and time for vector
 *  @param [in] m Mesh instance
 *  @param [in] name Name of the vector
 *  @param [in] t Time value for the vector
 *  @param [in] nb_dof Number of degrees of freedom per node, element or side
 *  If <tt>nb_dof</tt> is set to <tt>0</tt> the constructor picks this number from the Mesh instance
 *  @param [in] dof_type Type of degrees of freedom. To be given among the enumerated
 *  values: <tt>NODE_FIELD</tt>, <tt>ELEMENT_FIELD</tt>, <tt>SIDE_FIELD</tt> or <tt>EDGE_FIELD</tt>
 *  (Default: <tt>NODE_FIELD</tt>)
 */
    Vect(Mesh&  m,
         string name,
         real_t t=0.0,
         int    nb_dof=0,
         int    dof_type=NODE_FIELD);

/** \brief Constructor of an element vector.
 *  \details The constructed vector has local numbering of nodes
 *  @param [in] el Pointer to Element to localize
 *  @param [in] v Global vector to localize
 */
    Vect(const Element*  el,
         const Vect<T_>& v);

/** \brief Constructor of a side vector
 *  \details The constructed vector has local numbering of nodes
 *  @param [in] sd Pointer to Side to localize
 *  @param [in] v Global vector to localize
 */
    Vect(const Side*     sd,
         const Vect<T_>& v);

/** \brief Constructor using boundary conditions
 *  \details Boundary condition values contained in <tt>bc</tt> are reported to vector <tt>v</tt>
 *  @param [in] v Vect instance to update
 *  @param [in] bc Vect instance containing imposed valued at desired DOF
 */
    Vect(const Vect<T_>& v,
         const Vect<T_>& bc);

/** \brief Constructor to select some components of a given vector.
 *  @param [in] v Vect instance to extract from
 *  @param [in] nb_dof Number of DOF to extract
 *  @param [in] first_dof First DOF to extract
 *  For instance, a choice <tt>first_dof=2</tt> and <tt>nb_dof=1</tt> means that
 *  the second DOF of each node is copied in the vector
 */
    Vect(const Vect<T_>& v,
         size_t          nb_dof,
         size_t          first_dof);

/// \brief Copy constructor
    Vect(const Vect<T_>& v);

/** \brief Constructor to select one component from a given 2 or 3-component vector
 *  @param [in] v Vect instance to extract from
 *  @param [in] n Component to extract (must be > 1 and < 4 or).
 */
    Vect(const Vect<T_>& v,
         size_t          n);

/** \brief Constructor that extracts some degrees of freedom (components) from given instance of Vect.
 *  \details This constructor enables constructing a subvector of a given
 *  Vect instance. It selects a given list of degrees of freedom
 *  and put it according to a given order in the instance to construct.
 *  @param [in] d Integer number giving the list of degrees of freedom. This
 *  number is made of <tt>n</tt> digits where <tt>n</tt> is the number of degrees of freedom.
 *  Let us give an example: Assume that the instance <tt>v</tt> has 3 DOF by entity (node, element 
 *  or side). The
 *  choice <tt>d=201</tt> means that the constructed instance has 2 DOF where the first
 *  DOF is the third one of <tt>v</tt>, and the second DOF is the first one of f <tt>v</tt>.
 *  Consequently, no digit can be larger than the number of DOF the constructed instance.
 *  In this example, a choice <tt>d=103</tt> would produce an error message.
 *  @param [in] v Vect instance from which extraction is performed.
 *  @param [in] name Name to assign to vector instance (Default value is " ").
 *  \warning Don't give zeros as first digits for the argument <tt>d</tt>. The number is
 *  in this case interpreted as octal !!
*/
    Vect(size_t          d,
         const Vect<T_>& v,
         const string&   name=" ");

#if defined (USE_EIGEN)
/** \brief Constructor that copies the vector from a Eigen Vector instance
 *  @param [in] v VectorX instance from which extraction is performed
 *  @warning This constructor is available only if the library <tt>eigen</tt> is used in 
 *  conjunction with OFELI
 *  @remark: This constructor is available only if the %Eigen library was installed in conjunction
 *  with OFELI
 */
    Vect(const VectorX& v);
#endif

/// \brief Destructor
    ~Vect();

/** \brief Initialize vector with a c-array
 *  @param [in] v c-array (pointer) to initialize Vect
 *  @param [in] n size of array
 */
    void set(const T_* v,
             size_t    n);

/** \brief Initialize vector with another Vect instance
 *  @param [in] v Vect instance to extract from
 *  @param [in] nb_dof Number of DOF per node, element or side (By default, 0: Number of degrees
 *  of freedom extracted from the Mesh instance)
 *  @param [in] first_dof First DOF to extract (Default: 1)
 *  For instance, a choice <tt>first_dof=2</tt> and <tt>nb_dof=1</tt> means that
 *  the second DOF of each node is copied in the vector
 */
    void select(const Vect<T_>& v,
                size_t          nb_dof=0,
                size_t          first_dof=1);

/** \brief Initialize vector with an algebraic expression
 *  \details This function is to be used is a Mesh instance is associated to the vector
 *  @param [in] exp Regular algebraic expression that defines a function of <tt>x</tt>,
 *  <tt>y</tt>, <tt>z</tt>  which are coordinates of nodes and <tt>t</tt> which is the time value.
 *  @param [in] dof Degree of freedom to which the value is assigned [Default: 1]
 *  @warning If the time variable <tt>t</tt> is involved in the expression, the time value
 *  associated to the vector instance must be defined (Default value is 0) either by using
 *  the appropriate constructor or by the member function setTime.
 */
    void set(const string& exp,
             size_t        dof=1);

/** \brief Initialize vector with an algebraic expression
 *  \details This function can be used for instance in 1-D
 *  @param [in] exp Regular algebraic expression that defines a function of <tt>x</tt>
 *  which are coordinates of nodes
 *  @param [in] x Vector
 */
    void set(const string&       exp,
             const Vect<real_t>& x);
   
/** \brief Initialize vector with an algebraic expression with providing mesh data
 *  @param [in] ms Mesh instance
 *  @param [in] exp Regular algebraic expression that defines a function of x, y and z
 *  which are coordinates of nodes.
 *  @param [in] dof Degree of freedom to which the value is assigned [Default: 1]
 */
    void set(Mesh&         ms,
             const string& exp,
             size_t        dof=1);

/** \brief Initialize vector with an algebraic expression 
 *  @param [in] x Vect instance that contains coordinates of points
 *  @param [in] exp Regular algebraic expression that defines a function of x and i
 *  which are coordinates of nodes and indices starting from <tt>1</tt>.
 */
    void set(const Vect<real_t>& x,
             const string&       exp);

/** \brief Define mesh class to size vector
 *  @param [in] m Mesh instance
 *  @param [in] nb_dof Number of degrees of freedom per node, element or side
 *  If <tt>nb_dof</tt> is set to <tt>0</tt> the constructor picks this number from the Mesh instance
 *  @param [in] dof_type Parameter to precise the type of degrees of freedom. To be chosen
 *  among the enumerated values: <tt>NODE_FIELD</tt>, <tt>ELEMENT_FIELD</tt>,
 *  <tt>SIDE_FIELD</tt>, <tt>EDGE_FIELD</tt> [Default: <tt>NODE_FIELD</tt>]
 */
    void setMesh(Mesh&  m,
                 size_t nb_dof=0,
                 size_t dof_type=NODE_FIELD);

#if defined (USE_EIGEN)
/// \brief Return vector (global) size
/// @warning This constructor is available only if the library <tt>eigen</tt> is used in 
/// conjunction with OFELI
    size_t size() const { return _size; }
#endif

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
    void resize(size_t n) { setSize(n); }

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
 *  the enumerated values: <tt>NODE_FIELD</tt>, <tt>ELEMENT_FIELD</tt>, <tt>SIDE_FIELD</tt> 
 *  or <tt>EDGE_FIELD</tt>
 */
    void setDOFType(int dof_type) { _dof_type=dof_type; }

/** \brief Set Discontinuous Galerkin type vector
 *  \details When the vector is associated to a mesh, this one is sized differently if the
 *  DG method is used.
 *  @param [in] degree Polynomial degree of the DG method [Default: <tt>1</tt>]
 */
    void setDG(int degree=1);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    bool getGrid() const { return _grid; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return vector number of degrees of freedom
    size_t getNbDOF() const { return _nb_dof; }

/// \brief Return vector number of entities (nodes, elements or sides)
    size_t getNb() const { return size()/_nb_dof; }

/// \brief Return Mesh instance
    Mesh &getMesh() const { return *_theMesh; }

/// \brief Return <tt>true</tt> if vector contains a Mesh pointer, <tt>false</tt> if not
/// \details A Vect instance can be constructed using mesh information 
    bool WithMesh() const { return _with_mesh; }

/** Return DOF type of vector
 *  @return dof_type Type of degrees of freedom. Value among the enumerated
 *  values: <tt>NODE_FIELD</tt>, <tt>ELEMENT_FIELD</tt>, <tt>SIDE_FIELD</tt> or <tt>EDGE_FIELD</tt>
 */
    int getDOFType() const { return int(_dof_type); }

/// \brief Set time value for vector
    void setTime(real_t t) { _time = t; }

/// \brief Get time value for vector
    real_t getTime() const { return _time; }

/// \brief Set name of vector
    void setName(string name) { _name = name; }

/// \brief Get name of vector
    string getName() const { return _name; }

/// \brief Calculate 1-norm of vector
/// @remark This function is available only if the template parameter is <tt>double</tt> or <tt>complex&lt;double&gt;</tt>
    real_t getNorm1() const;

/// \brief Calculate 2-norm (Euclidean norm) of vector
/// @remark This function is available only if the template parameter is <tt>double</tt> or <tt>complex&lt;double&gt;</tt>
    real_t getNorm2() const;

/// \brief Calculate Max-norm (Infinite norm) of vector
/// @remark This function is available only if the template parameter is <tt>double</tt> or <tt>complex&lt;double&gt;</tt>
    real_t getNormMax() const;

/// \brief Calculate weighted 1-norm of vector
/// The wighted 1-norm is the 1-Norm of the vector divided by its size
    real_t getWNorm1() const { return getNorm1()/size(); }

/** \brief Calculate weighted 2-norm of vector
 *  \details The weighted 2-norm is the 2-Norm of the vector divided by the
 *  square root of its size
 */
    real_t getWNorm2() const { return getNorm2()/sqrt(real_t(size())); }

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

/** \brief Assign a given function (given by an interpretable algebraic expression) to 
 *  components of vector corresponding to sides with given code.
 *  \details Vector components are assumed nodewise
 *  @param [in] m    Instance of mesh
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned [default: <tt>1</tt>]
 */
    void setSideBC(Mesh&         m,
                   int           code,
                   const string& exp,
                   size_t        dof=1);

/** \brief Assign a given value to components of vector with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] val  Value to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned [default: <tt>1</tt>]
 */
    void setNodeBC(int    code,
                   T_     val,
                   size_t dof=1) { setNodeBC(*_theMesh,code,val,dof); }

/** \brief Assign a given function (given by an interpretable algebraic expression) to 
 *  components of vector with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned [default: <tt>1</tt>]
 *  @warning This member function is to be used in the case where a constructor with a Mesh 
 *  has been used
 */
    void setNodeBC(int           code,
                   const string& exp,
                   size_t        dof=1) { setNodeBC(*_theMesh,code,exp,dof); }

/** \brief Assign a given function (given by an interpretable algebraic expression) to 
 *  components of vector with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned [default: <tt>1</tt>]
 *  @warning This member function is to be used in the case where a constructor with a Mesh 
 *  has been used
 */
    void setSideBC(int           code,
                   const string& exp,
                   size_t        dof=1) { setSideBC(*_theMesh,code,exp,dof); }

/** \brief Assign a given value to components of vector with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] val  Value to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned [default: <tt>1</tt>]
 *  @warning This member function is to be used in the case where a constructor with a Mesh 
 *  has been used
 */
    void setSideBC(int     code,
                   T_      val,
                   size_t  dof=1) { setSideBC(*_theMesh,code,val,dof); }

/** \brief Remove boundary conditions.
 *  \details This member function copies to current vector a vector
 *  where only non imposed DOF are retained.
 *  @param [in] ms Mesh instance
 *  @param [in] v %Vector (Vect instance to copy from)
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned (<tt>=0</tt>, Default) or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only one 
 *  degree of freedom
 */
    void removeBC(const Mesh&     ms,
                  const Vect<T_>& v,
                  int             dof=0);

/** \brief Remove boundary conditions.
 *  \details This member function copies to current vector a vector
 *  where only non imposed DOF are retained.
 *  @param [in] v Vector (Vect instance to copy from)
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned [Default: <tt>0</tt>] or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only one
 *  degree of freedom.
 *  @warning This member function is to be used in the case where a constructor with a Mesh 
 *  has been used
 */
    void removeBC(const Vect<T_>& v,
                  int             dof=0);

/** \brief Transfer boundary conditions to the vector
 *  @param [in] bc Vect instance from which imposed degrees of freedom are copied to current instance
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned (=0, Default) or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only one 
 *  degree of freedom.
 */
    void transferBC(const Vect<T_>& bc,
                    int             dof=0);

/** \brief Insert boundary conditions.
 *  @param [in] m Mesh instance.
 *  @param [in] v Vect instance from which free degrees of freedom are copied to current instance.
 *  @param [in] bc Vect instance from which imposed degrees of freedom are copied to current instance.
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned (=0, Default) or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only 
 *  one degree of freedom by node or side
 */
    void insertBC(Mesh&           m,
                  const Vect<T_>& v,
                  const Vect<T_>& bc,
                  int             dof=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if defined(USE_PETSC)
    void insertBC(Mesh&                m, 
                  const PETScVect<T_>& v,
                  const Vect<T_>&      bc,
                  int                  dof=0);
    void insertBC(Mesh&                m,
                  const PETScVect<T_>& v,
                  int                  dof=0);
    void insertBC(const PETScVect<T_>& v,
                  const Vect<T_>&      bc,
                  int                  dof=0);
    void insertBC(const PETScVect<T_>& v, 
                  int                  dof=0);
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Insert boundary conditions.
 *  \details DOF with imposed boundary conditions are set to zero.
 *  @param [in] m Mesh instance.
 *  @param [in] v Vect instance from which free degrees of freedom are copied to current instance.
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned (<tt>=0</tt>, Default) or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only 
 *  one degree of freedom by node or side
 */
    void insertBC(Mesh&           m,
                  const Vect<T_>& v,
                  int             dof=0);

/** \brief Insert boundary conditions.
 *  @param [in] v Vect instance from which free degrees of freedom are copied to current instance.
 *  @param [in] bc Vect instance from which imposed degrees of freedom are copied to current instance.
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned (=0, Default) or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only 
 *  one degree of freedom by node or side
 */
    void insertBC(const Vect<T_>& v,
                  const Vect<T_>& bc,
                  int             dof=0);

/** \brief Insert boundary conditions.
 *  \details DOF with imposed boundary conditions are set to zero.
 *  @param [in] v Vect instance from which free degrees of freedom are copied to current instance.
 *  @param [in] dof Parameter to say if all degrees of freedom are concerned (=0, Default) or
 *  if only one degree of freedom (<tt>dof</tt>) is inserted into vector <tt>v</tt> which has only 
 *  one degree of freedom by node or side
 *  @warning This member function is to be used in the case where a constructor with a Mesh 
 *  has been used
 */
    void insertBC(const Vect<T_>& v, 
                  int             dof=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void insertBC1(Mesh&           m, 
                   const Vect<T_>& v,
                   const Vect<T_>& bc);
    void insertBC1(const Vect<T_>& v,
                   const Vect<T_>& bc);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Assembly of element vector into current instance.
 *  @param [in] el Reference to Element instance
 *  @param [in] b Local vector to assemble (Instance of class Vect)
 */
    void Assembly(const Element&  el,
                  const Vect<T_>& b);

/** \brief Assembly of element vector (as C-array) into Vect instance.
 *  @param [in] el Reference to Element instance
 *  @param [in] b Local vector to assemble (C-Array)
 */
    void Assembly(const Element& el,
                  const T_*      b);

/** \brief Assembly of side vector into Vect instance.
 *  @param [in] sd Reference to Side instance
 *  @param [in] b Local vector to assemble (Instance of class Vect)
 */
    void Assembly(const Side&     sd,
                  const Vect<T_>& b);

/** \brief Assembly of side vector (as C-array) into Vect instance.
 *  @param [in] sd Reference to Side instance
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
 *  \details The resulting gradient is stored in a Vect instance.
 *  This function handles node vectors assuming P<sub>1</sub> approximation.
 *  The gradient is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the gradient, where <tt>v(n,1)</tt>,
 *  <tt>v(n,2)</tt> and <tt>v(n,3)</tt> are respectively the <tt>x</tt> and <tt>y</tt> and 
 *  <tt>z</tt> derivatives at element <tt>n</tt>.
 */
    void getGradient(class Vect<T_>& v);

/** \brief Evaluate the discrete Gradient vector of the current vector.
 *  \details The resulting gradient is stored in an Vect instance.
 *  This function handles node vectors assuming P<sub>1</sub> approximation.
 *  The gradient is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the gradient, where <tt>v(n,1).x</tt>,
 *  <tt>v(n,2).y</tt> and <tt>v(n,3).z</tt> are respectively the <tt>x</tt> and <tt>y</tt> and 
 *  <tt>z</tt> derivatives at element <tt>n</tt>.
 */
    void getGradient(Vect<Point<T_> >& v);

/** \brief Evaluate the discrete curl vector of the current vector.
 *  \details The resulting curl is stored in a Vect instance.
 *  This function handles node vectors assuming P<sub>1</sub> approximation.
 *  The curl is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the curl, where <tt>v(n,1)</tt>,
 *  <tt>v(n,2)</tt> and <tt>v(n,3)</tt> are respectively the <tt>x</tt> and <tt>y</tt> and
 *  <tt>z</tt> <tt>curl</tt> components at element <tt>n</tt>.
 */
    void getCurl(Vect<T_>& v);

/** \brief Evaluate the discrete curl vector of the current vector.
 *  \details The resulting curl is stored in a Vect instance.
 *  This function handles node vectors assuming P<sub>1</sub> approximation.
 *  The curl is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the curl, where <tt>v(n,1).x</tt>,
 *  <tt>v(n,2).y</tt> and <tt>v(n,3).z</tt> are respectively the <tt>x</tt> and <tt>y</tt> and 
 *  <tt>z</tt> <tt>curl</tt> components at element <tt>n</tt>.
 */
    void getCurl(Vect<Point<T_> >& v);

/** \brief Evaluate the discrete scalar curl in 2-D of the current vector.
 *  \details The resulting curl is stored in a Vect instance.
 *  This function handles node vectors assuming P<sub>1</sub> approximation.
 *  The curl is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the scalar curl.
 */
    void getSCurl(Vect<T_>& v);

/** \brief Evaluate the discrete Divergence of the current vector.
 *  \details The resulting divergence is stored in a Vect instance.
 *  This function handles node vectors assuming P<sub>1</sub> approximation.
 *  The divergence is then a constant vector for each element.
 *  @param [in] v Vect instance that contains the divergence.
 */
    void getDivergence(Vect<T_>& v);

/** \brief Return average value of vector in a given element
 *  @param [in] el Element instance
 *  @param [in] type Type of element. This is to be chosen
 *  among enumerated values:
 *  <tt>LINE2</tt>, <tt>TRIANG3</tt>, <tt>QUAD4</tt>
 *  <tt>TETRA4</tt>, <tt>HEXA8</tt>, <tt>PENTA6</tt>
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
 *  @param [in] x Vect instance to add
 *  @param [in] a Constant to multiply before adding
 */
    Vect<T_> &MultAdd(const Vect<T_>& x,
                      const T_&       a);

/** \brief Add to vector the product of a vector by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] x Vect instance by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(      T_        a,
              const Vect<T_>& x);

/// \brief Assign a value to an entry for a 1-D vector
/// @param [in] i Rank index in vector (starts at <tt>1</tt>)
/// @param [in] val Value to assign
    void set(size_t i,
             T_     val);

/** \brief Assign a value to an entry for a 2-D vector
 *  @param [in] i First index in vector (starts at <tt>1</tt>)
 *  @param [in] j Second index in vector (starts at <tt>1</tt>)
 *  @param [in] val Value to assign
 */
    void set(size_t i,
             size_t j,
             T_     val);

/** \brief Assign a value to an entry for a 3-D vector
 *  @param [in] i First index in vector (starts at <tt>1</tt>)
 *  @param [in] j Second index in vector (starts at <tt>1</tt>)
 *  @param [in] k Third index in vector (starts at <tt>1</tt>)
 *  @param [in] val Value to assign
 */
    void set(size_t i,
             size_t j,
             size_t k,
             T_     val);


/// \brief Add a value to an entry for a 1-index vector
/// @param [in] i Rank index in vector (starts at <tt>1</tt>)
/// @param [in] val Value to assign
    void add(size_t i,
             T_     val);

/** \brief Add a value to an entry for a 2-index vector
 *  @param [in] i First index in vector (starts at <tt>1</tt>)
 *  @param [in] j Second index in vector (starts at <tt>1</tt>)
 *  @param [in] val Value to assign
 */
    void add(size_t i,
             size_t j,
             T_     val);

/** \brief Assign a value to an entry for a 3-index vector
 *  @param [in] i First index in vector (starts at <tt>1</tt>)
 *  @param [in] j Second index in vector (starts at <tt>1</tt>)
 *  @param [in] k Third index in vector (starts at <tt>1</tt>)
 * @param [in] val Value to assign
 */
    void add(size_t i,
             size_t j,
             size_t k,
             T_     val);

/// \brief Clear vector: Set all its elements to zero
    void clear();

#if defined (USE_EIGEN)
/** \brief Operator <tt>[]</tt> (Non constant version)
 *  @param [in] i Rank index in vector (starts at <tt>0</tt>)
 */
    T_ &operator[](size_t i);
#endif

#if defined (USE_EIGEN)
/** \brief Operator <tt>[]</tt> (Constant version)
 *  @param [in] i Rank index in vector (starts at <tt>0</tt>)
 */
    T_ operator[](size_t i) const;
#endif

/** \brief Operator <tt>()</tt> (Non constant version)
 *  @param [in] i Rank index in vector (starts at <tt>1</tt>)
 *  <ul>
 *    <li><tt>v(i)</tt> starts at <tt>v(1)</tt> to <tt>v(size())</tt>
 *    <li><tt>v(i)</tt> is the same element as <tt>v[i-1]</tt>
 *  </ul>
 */
    T_ &operator()(size_t i);

/** \brief Operator <tt>()</tt> (Constant version)
 *  @param [in] i Rank index in vector (starts at <tt>1</tt>)
 *  <ul>
 *    <li><tt>v(i)</tt> starts at <tt>v(1)</tt> to <tt>v(size())</tt>
 *    <li><tt>v(i)</tt> is the same element as <tt>v[i-1]</tt>
 *  </ul>
 */
    T_ operator()(size_t i) const;

/** \brief Operator <tt>()</tt> with 2-D indexing (Non constant version, case of a grid vector).
 *  @param [in] i first index in vector (Number of vector components in the <tt>x</tt>-grid)
 *  @param [in] j second index in vector (Number of vector components in the <tt>y</tt>-grid)
 *  <tt>v(i,j)</tt> starts at <tt>v(1,1)</tt> to <tt>v(getNx(),getNy())</tt>
 */
    T_ &operator()(size_t i,
                   size_t j);

/** \brief Operator <tt>()</tt> with 2-D indexing (Constant version).
 *  @param [in] i first index in vector (Number of vector components in the <tt>x</tt>-grid)
 *  @param [in] j second index in vector (Number of vector components in the <tt>y</tt>-grid)
 *  <tt>v(i,j)</tt> starts at <tt>v(1,1)</tt> to <tt>v(getNx(),getNy())</tt>
 */
    T_ operator()(size_t i,
                  size_t j) const;

/** \brief Operator <tt>()</tt> with 3-D indexing (Non constant version).
 *  @param [in] i first index in vector (Number of vector components in the <tt>x</tt>-grid)
 *  @param [in] j second index in vector (Number of vector components in the <tt>y</tt>-grid)
 *  @param [in] k third index in vector (Number of vector components in the <tt>z</tt>-grid)
 *  <tt>v(i,j,k)</tt> starts at <tt>v(1,1,1)</tt> to <tt>v(getNx(),getNy(),getNz())</tt>
 */
    T_ &operator()(size_t i,
                   size_t j,
                   size_t k);

/** \brief Operator <tt>()</tt> with 3-D indexing (Constant version).
 *  @param [in] i first index in vector (Number of vector components in the <tt>x</tt>-grid)
 *  @param [in] j second index in vector (Number of vector components in the <tt>y</tt>-grid)
 *  @param [in] k third index in vector (Number of vector components in the <tt>z</tt>-grid)
 *  <tt>v(i,j,k)</tt> starts at <tt>v(1,1,1)</tt> to <tt>v(getNx(),getNy(),getNz())</tt>
 */
    T_ operator()(size_t i,
                  size_t j,
                  size_t k) const;

/// \brief Operator <tt>=</tt> between vectors
    Vect<T_> & operator=(const Vect<T_>& v);


#if defined (USE_EIGEN)
/** \brief Operator <tt>=</tt> for an instance of <tt>VectorX</tt>
 *  @param [in] v Instance of vector class in library <tt>Eigen</tt>
 *  @remark The Vect instance must have been sized before
 *  @remark This operator is available only if the %Eigen library was installed in conjunction
 *  with OFELI
 */
    Vect<T_> &operator=(const VectorX& v);
#endif

/** \brief Operator <tt>=</tt>
 * \details Assign an algebraic expression to vector entries. This operator
 * has the same effect as the member function set(s)
 * @param [in] s String defining the algebraic expression as a function 
 * of coordinates and time
 * @warning A Mesh instance must has been introduced before (<i>e.g.</i>
 * by a constructor)
 */
    void operator=(string s);

/** \brief Initialize vector entries by setting extremal values and interval
 *  @param [in] vmin Minimal value to assign to the first entry
 *  @param [in] delta Interval
 *  @param [in] vmax Maximal value to assign to the lase entry
 *  @remark Vector's size is deduced from the arguments. The vector does not need
 *  to be sized before using this function
 */
    void setUniform(T_ vmin,
                    T_ delta,
                    T_ vmax);

/** \brief Operator <tt>=</tt>
 *  \details Assign a constant to vector entries
 *  @param [in] a Value to set
 */
    Vect<T_> & operator=(const T_& a);

/** \brief Operator <tt>+=</tt>
 *  \details Add vector <tt>x</tt> to current vector instance.
 *  @param [in] v Vect instance to add to instance
 */
    Vect<T_> & operator+=(const Vect<T_>& v);

/** \brief Operator <tt>+=</tt>
 *  \details Add a constant to current vector entries.
 *  @param [in] a Value to add to vector entries
 */
    Vect<T_> & operator+=(const T_& a);

/// \brief Operator <tt>-=</tt>
/// @param [in] v Vect instance to subtract from
    Vect<T_> & operator-=(const Vect<T_>& v);

/** \brief Operator <tt>-=</tt>
 *  \details Subtract constant from vector entries.
 *  @param [in] a Value to subtract from
 */
    Vect<T_> & operator-=(const T_& a);

/// \brief Operator <tt>*=</tt>
/// @param [in] a Value to multiply by
    Vect<T_> &operator*=(const T_& a);

/// \brief Operator <tt>/=</tt>
/// @param [in] a Value to divide by
    Vect<T_> &operator/=(const T_& a);

/** \brief Add an entry to the vector
 *  \details This function is an overload of the member function push_back of
 *  the parent class vector. It adjusts in addition some vector parameters
 *  @param [in] v Entry value to add
 */
    void push_back(const T_& v);

/// \brief Return reference to Mesh instance
    const Mesh &getMeshPtr() const { return *_theMesh; }
    
/** \brief Return Dot (scalar) product of two vectors
 *  \details A typical use of this operator is <tt>double a = (v,w)</tt>
 *  where <tt>v</tt> and <tt>w</tt> are 2 instances of <tt>Vect<double></tt>
 *  @param [in] v Vect instance by which the current instance is multiplied
 */
    T_ operator,(const Vect<T_>& v) const;

#if defined (USE_EIGEN)
/** \brief Casting operator
 *  @warning This constructor is available only if the library <tt>eigen</tt> is used in 
 *  conjunction with OFELI
 */
    operator VectorX() const { return _v; }
#endif

 private:
    size_t      _size, _nx, _ny, _nz, _nb, _dof_type, _nb_dof;
    int         _dg_degree;
    bool        _grid, _with_mesh;
    Mesh        *_theMesh;
    string      _name;
    real_t      _time;
    void dof_select(size_t d, vector<size_t> &dof_list);
    size_t ijk(size_t i, size_t j)           const { return _ny*(i-1)+j-1; }
    size_t ijk(size_t i, size_t j, size_t k) const { return _ny*_nz*(i-1)+_nz*(j-1)+k-1; }
#if defined (USE_EIGEN)
    VectorX     _v;
#endif

};


///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_>
Vect<T_>::Vect() :
#if !defined (USE_EIGEN)
          vector<T_>(),
#endif
          _size(0), _nx(0), _ny(1), _nz(1), _dof_type(NONE), _nb_dof(1),
          _dg_degree(-1), _grid(true), _with_mesh(false),
          _theMesh(NULL), _name("#"), _time(0)
{ }


template<class T_>
Vect<T_>::Vect(size_t n) :
#if !defined (USE_EIGEN)
           vector<T_>(n),
#endif
           _size(n), _nx(n), _ny(1), _nz(1),
           _dof_type(NONE), _nb_dof(1), _dg_degree(-1),
           _grid(true), _with_mesh(false), _theMesh(NULL), _name("#"), _time(0)
{
#if defined (USE_EIGEN)
   _v.conservativeResize(_size);
   clear();
#else
   for (size_t i=0; i<_size; i++)
      (*this)[i] = 0;
#endif
}


template<class T_>
Vect<T_>::Vect(size_t nx, 
               size_t ny) :
#if !defined (USE_EIGEN)
           vector<T_>(nx*ny),
#endif
           _size(nx*ny), _nx(nx), _ny(ny), _nz(1),
           _dof_type(NONE), _nb_dof(1), _dg_degree(-1),
           _grid(true), _with_mesh(false), _theMesh(NULL), _name("#"), _time(0)
{
#if defined (USE_EIGEN)
   _v.conservativeResize(_size);
   clear();
#else
   for (size_t i=0; i<_size; i++)
      (*this)[i] = 0;
#endif
}


template<class T_>
Vect<T_>::Vect(size_t nx,
               size_t ny,
               size_t nz) :
#if !defined (USE_EIGEN)
           vector<T_>(nx*ny*nz),
#endif
           _size(nx*ny*nz), _nx(nx), _ny(ny), _nz(nz),
           _dof_type(NONE), _nb_dof(1), _dg_degree(-1),
           _grid(true), _with_mesh(false), _theMesh(NULL), _name("#"), _time(0)
{
#if defined (USE_EIGEN)
   _v.conservativeResize(_size);
   clear();
#else
   for (size_t i=0; i<_size; i++)
      (*this)[i] = 0;
#endif
}


template<class T_>
Vect<T_>::Vect(size_t n,
               T_*    x) :
#if !defined (USE_EIGEN)
           vector<T_>(n),
#endif
           _size(n), _nx(n), _ny(1), _nz(1),
           _dof_type(NONE), _nb_dof(1), _dg_degree(-1),
           _grid(true), _with_mesh(false), _theMesh(NULL), _name("#"), _time(0)
{
#if defined (USE_EIGEN)
   _v.conservativeResize(_size);
#endif
   for (size_t i=1; i<=n; ++i)
      set(i,x[i-1]);
}


template<class T_>
Vect<T_>::Vect(Mesh& m,
               int   nb_dof,
               int   dof_type)
         : _dg_degree(-1), _grid(false), _with_mesh(true), _name("#"), _time(0)
{
   setMesh(m,nb_dof,dof_type);
}


template<class T_>
Vect<T_>::Vect(Mesh&  m,
               string name,
               real_t t,
               int    nb_dof,
               int    dof_type) :
#if !defined (USE_EIGEN)
          vector<T_>(),
#endif
          _dg_degree(-1), _grid(false), _with_mesh(true), _name("#"), _time(t)
{
   setMesh(m,nb_dof,dof_type);
}


template<class T_>
Vect<T_>::Vect(const Element*  el,
               const Vect<T_>& v) :
#if !defined (USE_EIGEN)
          vector<T_>(),
#endif
           _nx(el->getNbNodes()), _ny(v._ny), _nz(1),
           _nb(el->getNbNodes()), _nb_dof(v._nb_dof), _dg_degree(-1),
           _grid(false), _with_mesh(false), _name(v._name), _time(v._time)
{
   setSize(_nx,_ny);
   for (size_t n=1; n<=el->getNbNodes(); ++n) {
      Node *nd=(*el)(n);
      for (size_t j=1; j<=nd->getNbDOF(); ++j)
         set(n,j,v(nd->n(),j));
   }
}


template<class T_>
Vect<T_>::Vect(const Side*     sd,
               const Vect<T_>& v) :
#if !defined (USE_EIGEN)
           vector<T_>(),
#endif
           _nx(sd->getNbNodes()), _ny(v._nb_dof), _nz(1),
           _nb(sd->getNbNodes()), _nb_dof(v._nb_dof), _dg_degree(-1),
           _grid(false), _with_mesh(false), _name(v._name), _time(v._time)
{
   setSize(_nx,_ny);
   size_t i=0;
   for (size_t n=1; n<=sd->getNbNodes(); ++n) {
      Node *nd=(*sd)(n);
      size_t k=nd->getFirstDOF()-1;
      for (size_t j=1; j<=nd->getNbDOF(); ++j)
         set(++i,v[k++]);
   }
}


template<class T_>
Vect<T_>::Vect(const Vect<T_>& v,
               const Vect<T_>& bc) :
#if !defined (USE_EIGEN)
           vector<T_>(bc.size()),
#endif
           _size(v._nb*v._nb), _nx(v._nb), _ny(v._nb_dof), _nz(1),
           _nb(v._nb), _dof_type(v._dof_type), _nb_dof(v._nb_dof),
           _dg_degree(v._dg_degree), _grid(v._grid), _with_mesh(v._with_mesh),
           _theMesh(v._theMesh), _name(v._name), _time(v._time)
{
   size_t i=1, n=0;
#if defined (USE_EIGEN)
   resize(bc.size());
#endif
   mesh_nodes(*_theMesh) {
      for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
         set(i,bc[n++]);
         if (The_node.getCode(k) == 0)
            set(i,v[The_node.getDOF(k)-1]);
         i++;
      }
   }
}


template<class T_>
Vect<T_>::Vect(const Vect<T_>& v,
                     size_t    nb_dof,
                     size_t    first_dof)
         : _size(v._size), _nx(v._nx), _ny(v._ny), _nz(v._nz),
           _nb(v._nb), _dof_type(v._dof_type), _nb_dof(v._nb_dof),
           _dg_degree(v._dg_degree), _grid(v._grid), _with_mesh(v._with_mesh),
           _theMesh(v._theMesh), _name(v._name), _time(v._time)
{
   setSize(_nb,_nb_dof,1);
   for (size_t i=1; i<=_nb; i++)
      for (size_t j=1; j<=_nb_dof; j++)
         set(i,j,v(i,j+first_dof-1));
}


template<class T_>
Vect<T_>::Vect(const Vect<T_>& v)
         : _nb(v._nb), _dof_type(v._dof_type), _nb_dof(v._nb_dof),
           _dg_degree(v._dg_degree), _grid(v._grid), _with_mesh(v._with_mesh),
           _theMesh(v._theMesh), _name(v._name), _time(v._time)
{
   setSize(v._nx,v._ny,v._nz);
   for (size_t i=1; i<=_size; i++)
      set(i,v[i-1]);
}


template<class T_>
Vect<T_>::Vect(const Vect<T_>& v,
               size_t          n)
         : _nb_dof(v._nb_dof), _grid(v._grid), _with_mesh(v._with_mesh),
           _theMesh(v._theMesh), _name(v._name), _time(v._time)
{
   if (n==1) {
      setSize(v._nx,1,1);
      for (size_t i=1; i<=_nx; i++)
         set(i,v(i,1,1));
   }
   else if (n==2) {
      setSize(1,v._ny,1);
      for (size_t j=1; j<=_ny; j++)
         set(j,v(1,j,1));
   }
   else if (n==3) {
      setSize(1,1,v._nz);
      for (size_t k=1; k<=_nz; k++)
         set(k,v(1,1,k));
   }
}


template <class T_>
Vect<T_>::Vect(size_t          d,
               const Vect<T_>& v,
               const string&   name)
         : _nb(v._nb), _dof_type(v._dof_type), _with_mesh(v._with_mesh),
           _theMesh(v._theMesh), _name(name), _time(v._time)
{
   try {
     if (d<=0)
        THROW_RT("Vect(size_t,Vect<T_>,string): Illegal value of nb_dof = "+itos(d));
   }
   CATCH("Vect");
   size_t nd=v.getNbDOF();
   vector<size_t> dof_list(nd);
   dof_select(d,dof_list);
   try {
      if (_nb_dof>nd)
         THROW_RT("Vect(size_t,Vect<T_>,string): Illegal value of dof = "+itos(nd));
   }
   CATCH("Vect");
   _theMesh = &(v.getMesh());
   setSize(_nb,_nb_dof,1);
   for (size_t i=1; i<=_nb; i++) {
      for (size_t k=0; k<nd; k++) {
         if (dof_list[k]!=0)
            set(i,dof_list[k],v(i,k+1));
      }
   }
}


#if defined (USE_EIGEN)
template<class T_>
Vect<T_>::Vect(const VectorX& v)
         : _size(v.size()), _nx(_size), _ny(1), _nz(1), _dof_type(NONE), _nb_dof(1),
           _grid(true), _with_mesh(false), _theMesh(NULL), _name("#"), _time(0)
{ }
#endif


template<class T_>
Vect<T_>::~Vect()
{ }


template<class T_>
void Vect<T_>::setMesh(Mesh&  m,
                       size_t nb_dof,
                       size_t dof_type)
{
   _theMesh = &m;
   _with_mesh = true;
   size_t n=_theMesh->getNbDOF();
   _nb_dof = nb_dof;
   _dof_type = dof_type;
   if (dof_type==NODE_FIELD || dof_type==NODE_DOF)
      _nb = _theMesh->getNbNodes();
   else if (dof_type==SIDE_FIELD || dof_type==SIDE_DOF)
      _nb = _theMesh->getNbSides();
   else if (dof_type==ELEMENT_FIELD || dof_type==ELEMENT_DOF)
      _nb = _theMesh->getNbElements();
   if (nb_dof==0)
      _nb_dof = n/_nb;
   _size = _nb_dof * _nb;
   setSize(_nb,_nb_dof,1);
   clear();
}


template<class T_>
void Vect<T_>::set(const T_* v,
                   size_t    n)
{
   setSize(n);
   for (size_t i=1; i<=n; ++i)
      set(i,v[i-1]);
}


template<class T_>
void Vect<T_>::setDG(int degree)
{
   try {
      if (_with_mesh==false)
         THROW_RT("setDG(int): To be used only if a Mesh instance is associated to the vector");
   }
   CATCH("Vect");
   _dg_degree = degree;
   if (degree<0)
      return;
   _nb_dof = 0;
   _dof_type = ELEMENT_FIELD;
   switch (_theMesh->getDim()) {

      case 1: _nb_dof = _dg_degree+1;
              break;

      case 2: if (_dg_degree<10)
                 _nb_dof = (_dg_degree+1)*(_dg_degree+2)/2;
              else
                 _nb_dof = (_dg_degree+1)*(_dg_degree+1);
              break;

      case 3: if (_dg_degree<10)
                 _nb_dof = (_dg_degree+1)*(_dg_degree+2)/2;
              else if (_dg_degree<20)
                 _nb_dof = (_dg_degree+1)*(_dg_degree+1)*(_dg_degree+1);
              else
                 _nb_dof = (_dg_degree+1)*(_dg_degree+1)*(_dg_degree+2)/2;
              break;
   }
   setSize(_theMesh->getNbElements(),_nb_dof,1);
}


template<class T_>
void Vect<T_>::select(const Vect<T_>& v,
                      size_t          nb_dof,
                      size_t          first_dof)
{
   _size = nb_dof*v._nb;
   setSize(_size);
   size_t i=1;
   for (size_t n=1; n<=v._nb; ++n)
      for (size_t j=first_dof; j<=nb_dof+first_dof-1; j++)
         set(i++,v(n,j));
}


template <class T_>
void Vect<T_>::set(const string&       exp,
                   const Vect<real_t>& x)
{
   int err;
   real_t d[4];
   PARSE(exp.c_str(),"x");
   try {
      for (size_t i=0; i<x.size(); i++) {
         d[0] = x[i];
         set(i+1,EVAL(d));
         try {
            if ((err=theParser.EvalError()))
               THROW_RT("set(string,dof): Illegal regular expression. Error code: "+itos(err));
         }
         CATCH("Vect");
      }
   }
   CATCH("Vect");
}


template <class T_>
void Vect<T_>::set(const string& exp,
                   size_t        dof)
{
   try {
      if (_theMesh==NULL)
         THROW_RT("set(string,dof): No mesh is defined");
   }
   CATCH_EXIT("Vect");
   int err;
   real_t d[4];
   PARSE(exp.c_str(),"x,y,z,t");
   try {
      if (_dof_type==NODE_FIELD) {
         mesh_nodes(*_theMesh) {
            d[0] = The_node.getCoord(1);
            d[1] = The_node.getCoord(2);
            d[2] = The_node.getCoord(3);
            d[3] = _time;
            set(The_node.getNbDOF()*(node_label-1)+dof,EVAL(d));
            try {
               if ((err=theParser.EvalError()))
                  THROW_RT("set(string,dof): Illegal regular expression. Error code: " + itos(err));
            }
            CATCH("Vect");
         }
      }
      else if (_dof_type==ELEMENT_FIELD) {
      }
      else
         THROW_RT("set(string,size_t): This member function is for nodewise vectors only.");
   }
   CATCH("Vect");
}


template <class T_>
void Vect<T_>::set(Mesh&  ms,
                   const  string& exp,
                   size_t dof)
{
   setMesh(ms);
   try {
      if (_theMesh==NULL)
         THROW_RT("set(ms,string,dof): No mesh is defined");
   }
   CATCH_EXIT("Vect");
   int err;
   real_t d[4];
   PARSE(exp.c_str(),"x,y,z,t");
   try {
      if (_dof_type==NODE_FIELD) {
         mesh_nodes(*_theMesh) {
            d[0] = The_node.getCoord(1);
            d[1] = The_node.getCoord(2);
            d[2] = The_node.getCoord(3);
            d[3] = _time;
            set(The_node.getNbDOF()*(node_label-1)+dof,EVAL(d));
            try {
               if ((err=theParser.EvalError()))
                  THROW_RT("set(string,dof): Illegal regular expression. Error code: " + itos(err));
            }
            CATCH("Vect");
         }
      }
      else
         THROW_RT("set(string,size_t): This member function is for nodewise vectors only.");
   }
   CATCH("Vect");
}


template <class T_>
void Vect<T_>::set(const Vect<real_t>& x,
                   const string&       exp)
{
   int err;
   setSize(x._nx,x._ny,x._nz);
   real_t d[2];
   theParser.Parse(exp.c_str(),"x,i");
   for (size_t i=0; i<_size; i++) {
      d[0] = x[i];
      d[1] = i + 1;
      set(i+1,theParser.Eval(d));
      try {
         if ((err=theParser.EvalError()))
            THROW_RT("set(Vect,string): Illegal regular expression. Error code: " + itos(err));
      }
      CATCH("Vect");
   }
}


template<class T_>
void Vect<T_>::resize(size_t n,
                      T_     v)
{
   _nx = n, _ny = 1, _nz = 1;
   _size = _nx*_ny*_nz;
#if defined (USE_EIGEN)
   _v.conservativeResize(_size,1);
   for (size_t i=1; i<=_size; i++)
      set(i,v);
#else
   vector<T_>::resize(_size,v);
#endif
}


template<class T_>
void Vect<T_>::setSize(size_t nx,
                       size_t ny,
                       size_t nz)
{
   _nx = nx, _ny = ny, _nz = nz;
   _size = _nx*_ny*_nz;
#if defined (USE_EIGEN)
   _v.conservativeResize(_size,1);
#else
   vector<T_>::resize(_size);
#endif
}


template<>
inline void Vect<real_t>::operator=(string s)
{
   try {
      if (_theMesh==NULL)
         THROW_RT("operator=(string): No mesh is defined");
   }
   CATCH_EXIT("Vect");
   set(s);
}


template<class T_>
void Vect<T_>::setUniform(T_ vmin,
                          T_ delta,
                          T_ vmax)
{
   setSize((vmax-vmin)/delta+delta+1);
   (*this)[0] = vmin;
   for (size_t i=1; i<_nx; i++)
      (*this)[i] = (*this)[i-1] + delta;
}


template<>
inline void Vect<real_t>::setSize(size_t nx,
                                  size_t ny,
                                  size_t nz)
{
   _nx = nx, _ny = ny, _nz = nz;
   _size = _nx*_ny*_nz;
#if defined (USE_EIGEN)
   _v.conservativeResize(_size,1);
#else
   vector<real_t>::resize(_size);
#endif
   clear();
}


template<>
inline void Vect<complex_t>::setSize(size_t nx,
                                     size_t ny,
                                     size_t nz)
{
   _nx = nx, _ny = ny, _nz = nz;
   _size = _nx*_ny*_nz;
#if defined (USE_EIGEN)
   _v.conservativeResize(_size,1);
#else
   vector<complex_t>::resize(_size);
#endif
   clear();
}


template<>
inline real_t Vect<real_t>::getNorm1() const
{
   real_t s=0.;
   for (size_t i=0; i<size(); ++i)
      s += std::abs((*this)[i]);
   return s;
}


template<>
inline real_t Vect<complex_t>::getNorm1() const
{
   real_t s=0.;
   for (size_t i=0; i<size(); ++i)
      s += std::abs((*this)[i]);
   return s;
}


template<>
inline real_t Vect<real_t>::getNorm2() const
{
   real_t s=0;
   for (size_t i=0; i<size(); ++i)
      s += (*this)[i]*(*this)[i];
   s = std::sqrt(s);
   return s;
}


template<>
inline real_t Vect<complex_t>::getNorm2() const
{
   real_t s=0;
   for (size_t i=0; i<size(); ++i) {
      complex_t z = (*this)[i];
      s += z.real()*z.real() + z.imag()*z.imag();
   }
   return std::sqrt(s);
}


template<>
inline real_t Vect<real_t>::getNormMax() const
{
   real_t s=0;
   for (size_t i=0; i<size(); ++i) {
      real_t z = std::abs((*this)[i]);
      s = (z > s ? z : s);
   }
   return s;
}


template<>
inline real_t Vect<complex_t>::getNormMax() const
{
   real_t s=0;
   for (size_t i=0; i<size(); ++i) {
      real_t z = std::abs((*this)[i]);
      s = (z > s ? z : s);
   }
   return s;
}


template<class T_>
T_ Vect<T_>::getMin() const
{
   T_ s = (*this)[0];
   for (size_t i=1; i<size(); ++i)
      s = (*this)[i] < s ? (*this)[i] : s;
   return s;
}


template<class T_>
T_ Vect<T_>::getMax() const
{
  T_ s = (*this)[0];
   for (size_t i=1; i<size(); ++i)
      s = (*this)[i] > s ? (*this)[i] : s;
   return s;
}


template<class T_>
void Vect<T_>::removeBC(const Mesh&     ms,
                        const Vect<T_>& v,
                              int       dof)
{
   if (dof==0) {
      size_t n = 1;
      mesh_nodes(ms) {
         for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
            if (The_node.getCode(k) == 0)
               set(The_node.getDOF(k),v(n));
            n++;
         }
      }
   }
   else {
      mesh_nodes(ms) {
         if (The_node.getCode(dof) == 0)
            set(The_node.getDOF(dof),v(node_label));
      }
   }
}


template<class T_>
void Vect<T_>::removeBC(const Vect<T_>& v,
                        int             dof)
{
   if (dof==0) {
      size_t n = 1;
      mesh_nodes(*_theMesh) {
         for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
            if (The_node.getCode(k) == 0)
               set(The_node.getDOF(k),v(n));
            n++;
         }
      }
   }
   else {
      mesh_nodes(*_theMesh) {
         if (The_node.getCode(dof) == 0)
            set(The_node.getDOF(dof),v(node_label));
      }
   }
}


template<class T_>
void Vect<T_>::transferBC(const Vect<T_>& bc,
                          int             dof)
{
   size_t i=1, k=1;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         mesh_nodes(*_theMesh) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               if (The_node.getCode(k)>0)
                  set(i,bc(i));
               i++;
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         mesh_nodes(*_theMesh) {
            for (size_t k=1; k<=TheSide.getNbDOF(); ++k) {
               if (TheSide.getCode(k)>0)
                  set(i,bc(i));
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (_theMesh->NodesAreDOF()) {
         mesh_nodes(*_theMesh) {
            if (The_node.getCode(dof)>0)
               set(i,bc(k));
            i++;
            k += The_node.getNbDOF();
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         mesh_nodes(*_theMesh) {
            if (TheSide.getCode(dof)>0)
               set(i,bc(k));
            i++, k += The_node.getNbDOF();
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::insertBC(Mesh&           m,
                        const Vect<T_>& v,
                        const Vect<T_>& bc,
                        int             dof)
{
   size_t i=1;
   if (dof==0) {
      if (m.NodesAreDOF()) {
         mesh_nodes(m) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               else
                  set(i,bc(i));
               i++;
            }
         }
      }
      else if (m.SidesAreDOF()) {
         mesh_sides(m) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               if (The_side.getCode(k)>=0)
                  set(i,v(The_side.getDOF(k)));
               else
                  set(i,bc(i));
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
         mesh_nodes(m) {
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_node.getNbDOF();
         }
      }
      else if (m.SidesAreDOF()) {
         size_t k=dof;
         mesh_sides(m) {
            if (The_side.getCode(dof)>=0)
               set(i,v(The_side.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_side.getNbDOF();
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::insertBC(Mesh&           m,
                        const Vect<T_>& v,
                        int             dof)
{
   size_t i=1;
   if (dof==0) {
      if (m.NodesAreDOF()) {
         mesh_nodes(m) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               set(i,0);
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               i++;
            }
         }
      }
      else if (m.SidesAreDOF()) {
         mesh_sides(m) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               set(i,0);
               if (The_side.getCode(k)==0)
                  set(i,v(The_side.getDOF(k)));
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (m.NodesAreDOF()) {
         mesh_nodes(m) {
            set(i,0);
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            i++;
         }
      }
      else if (m.SidesAreDOF()) {
         mesh_sides(m) {
            set(i,0);
            if (The_side.getCode(dof)==0)
               set(i,v(The_side.getDOF(dof)));
            i++;
         }
      }
      else
         ;
   }
}


#if defined(USE_PETSC)
template<class T_>
void Vect<T_>::insertBC(Mesh&                m,
                        const PETScVect<T_>& v,
                        const Vect<T_>&      bc,
                        int                  dof)
{
   size_t i=1;
   if (dof==0) {
      if (m.NodesAreDOF()) {
         mesh_nodes(m) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               else
                  set(i,bc(i));
               i++;
            }
         }
      }
      else if (m.SidesAreDOF()) {
         mesh_sides(m) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               if (The_side.getCode(k)>=0)
                  set(i,v(The_side.getDOF(k)));
               else
                  set(i,bc(i));
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
         mesh_nodes(m) {
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_node.getNbDOF();
         }
      }
      else if (m.SidesAreDOF()) {
         size_t k=dof;
         mesh_sides(m) {
            if (The_side.getCode(dof)>=0)
               set(i,v(The_side.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_side.getNbDOF();
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::insertBC(Mesh&                m,
                        const PETScVect<T_>& v,
                        int                  dof)
{
   size_t i=1;
   if (dof==0) {
      if (m.NodesAreDOF()) {
         mesh_nodes(m) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
	       set(i,0);
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               i++;
            }
         }
      }
      else if (m.SidesAreDOF()) {
         mesh_sides(m) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               set(i,0);
               if (The_side.getCode(k)==0)
                  set(i,v(The_side.getDOF(k)));
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (m.NodesAreDOF()) {
         mesh_nodes(m) {
            set(i,0);
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            i++;
         }
      }
      else if (m.SidesAreDOF()) {
         mesh_sides(m) {
            set(i,0);
            if (The_side.getCode(dof)==0)
               set(i,v(The_side.getDOF(dof)));
            i++;
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::insertBC(const PETScVect<T_>& v,
                        const Vect<T_>&      bc,
                        int                  dof)
{
   size_t i = 1;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         mesh_nodes(*_theMesh) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               else
                  set(i,bc(i));
               i++;
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         mesh_sides(*_theMesh) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               if (The_side.getCode(k)>=0)
                  set(i,v(The_side.getDOF(k)));
               else
                  set(i,bc(i));
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
         mesh_nodes(*_theMesh) {
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_node.getNbDOF();
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         size_t k=dof;
         mesh_sides(*_theMesh) {
            if (The_side.getCode(dof)>=0)
               set(i,v(The_side.getDOF(dof)));
            else
               set(i,bc(The_side.getDOF(dof)));
            i++, k += The_side.getNbDOF();
         }
      }
      else
         ;
   }
}
#endif


template<class T_>
void Vect<T_>::insertBC1(Mesh&           m,
                         const Vect<T_>& v,
                         const Vect<T_>& bc)
{
   size_t i=1, n=0, l=0;
   mesh_nodes(m) {
      for (size_t k=1; k<=the_node->getNbDOF(); ++k) {
         set(i,bc[n++]);
         if (The_node.getCode(k)==0)
            set(i,v[l++]);
         i++;
      }
   }
}


template<class T_>
void Vect<T_>::insertBC(const Vect<T_>& v,
                        const Vect<T_>& bc,
                        int             dof)
{
   size_t i = 1;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         mesh_nodes(*_theMesh) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               else
                  set(i,bc(i));
               i++;
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         mesh_sides(*_theMesh) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               if (The_side.getCode(k)>=0)
                  set(i,v(The_side.getDOF(k)));
               else
                  set(i,bc(i));
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
         mesh_nodes(*_theMesh) {
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            else
               set(i,bc(k));
            i++, k += The_node.getNbDOF();
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         size_t k=dof;
         mesh_sides(*_theMesh) {
            if (The_side.getCode(dof)>=0)
               set(i,v(The_side.getDOF(dof)));
            else
               set(i,bc(The_side.getDOF(dof)));
            i++, k += The_side.getNbDOF();
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::setNodeBC(Mesh&         m,
                         int           code,
                         const string& exp,
                         size_t        dof)
{
   int err;
   real_t d[4];
   theParser.Parse(exp.c_str(),"x,y,z,t");
   d[3] = _time;
   mesh_nodes(m) {
      d[0] = The_node.getCoord(1);
      d[1] = The_node.getCoord(2);
      d[2] = The_node.getCoord(3);
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
         if (The_node.getCode(dof)==code) {
            set(node_label,dof,theParser.Eval(d));
            try {
               if ((err=theParser.EvalError()))
                  THROW_RT("setNodeBC(Mesh,int,string,size_t): Illegal regular expression. Error code: " + itos(err));
            }
            CATCH("Vect");
         }
      }
   }
}


template<class T_>
void Vect<T_>::setSideBC(Mesh&   m,
                         int     code,
                         T_      val,
                         size_t  dof)
{
   int err;
   mesh_sides(m) {
     if (The_side.getCode(dof)==code) {
         for (size_t n=1; n<=The_side.getNbNodes(); n++) {
            the_node = The_side(n);
            for (size_t i=1; i<=The_side.getNbDOF(); i++) {
               set(node_label,dof,val);
               try {
                  if ((err=theParser.EvalError()))
                     THROW_RT("setSideBC(Mesh,int,T_,size_t): Illegal regular expression. Error code: " + itos(err));
               }
               CATCH("Vect");
            }
         }
      }
   }
}


template<class T_>
void Vect<T_>::setSideBC(Mesh&         m,
                         int           code,
                         const string& exp,
                         size_t        dof)
{
   int err;
   real_t d[4];
   theParser.Parse(exp.c_str(),"x,y,z,t");
   d[3] = _time;
   mesh_sides(m) {
     if (The_side.getCode(dof)==code) {
         for (size_t n=1; n<=The_side.getNbNodes(); n++) {
            the_node = The_side(n);
            d[0] = The_node.getCoord(1);
            d[1] = The_node.getCoord(2);
            d[2] = The_node.getCoord(3);
            for (size_t i=1; i<=The_side.getNbDOF(); i++) {
               set(node_label,dof,theParser.Eval(d));
               try {
                  if ((err=theParser.EvalError()))
                     THROW_RT("setSideBC(Mesh,int,string,size_t): Illegal regular expression. Error code: " + itos(err));
               }
               CATCH("Vect");
            }
         }
      }
   }
}


template<class T_>
void Vect<T_>::setNodeBC(Mesh&  m,
                         int    code,
                         T_     val,
                         size_t dof)
{
   if (m.getDOFSupport()==NODE_FIELD) {
      mesh_nodes(m) {
         for (size_t i=1; i<=the_node->getNbDOF(); i++) {
            if (The_node.getCode(dof)==code)
               set(node_label,dof,val);
         }
      }
   }
   if (m.getDOFSupport()==SIDE_FIELD) {
      MeshBoundarySides(m) {
         for (size_t i=1; i<=theSide->getNbDOF(); i++) {
            if (theSide->getCode(dof)==code)
               set(side_label,dof,val);
         }
      }
   }
}


template<class T_>
void Vect<T_>::insertBC(const Vect<T_>& v,
                        int             dof)
{
   size_t i=1;
   if (dof==0) {
      if (_theMesh->NodesAreDOF()) {
         mesh_nodes(*_theMesh) {
            for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
               set(i,0);
               if (The_node.getCode(k)==0)
                  set(i,v(The_node.getDOF(k)));
               i++;
            }
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         mesh_sides(*_theMesh) {
            for (size_t k=1; k<=The_side.getNbDOF(); ++k) {
               set(i,0);
               if (The_side.getCode(k)==0)
                  set(i,v(The_side.getDOF(k)));
               i++;
            }
         }
      }
      else
         ;
   }
   else {
      if (_theMesh->NodesAreDOF()) {
         mesh_nodes(*_theMesh) {
            set(i,0);
            if (The_node.getCode(dof)==0)
               set(i,v(The_node.getDOF(dof)));
            i++;
         }
      }
      else if (_theMesh->SidesAreDOF()) {
         mesh_sides(*_theMesh) {
            set(i,0);
            if (The_side.getCode(dof)==0)
               set(i,v(The_side.getDOF(dof)));
            i++;
         }
      }
      else
         ;
   }
}


template<class T_>
void Vect<T_>::insertBC1(const Vect<T_>& v,
                         const Vect<T_>& bc)
{
   size_t i=1, n=0, l=0;
   mesh_nodes(*_theMesh) {
      for (size_t k=1; k<=The_node.getNbDOF(); ++k) {
         set(i,bc[n++]);
         if (The_node.getCode(k)==0)
            set(i,v[l++]);
         i++;
      }
   }
}


template<class T_>
void Vect<T_>::DGAssembly(const Element&                          el,
                          const LocalVect<T_,MAX_NB_ELEMENT_DOF>& b)
{
   for (size_t i=1; i<=el.getNbDOF(); ++i) {
      if (el.getDOF(i)!=0)
         add(el.getDOF(i),b(i));
   }
}


template<class T_>
void Vect<T_>::DGAssembly(const Side&                          sd,
                          const LocalVect<T_,MAX_NB_SIDE_DOF>& b)
{
   for (size_t i=1; i<=sd.getNbDOF(); ++i) {
      if (sd.getDOF(i)!=0)
         add(sd.getDOF(i),b(i));
   }
}


template<class T_>
void Vect<T_>::Assembly(const Element&  el,
                        const Vect<T_>& b)
{
   size_t i=1;
   for (size_t n=1; n<=el.getNbNodes(); ++n) {
      Node *nd=el(n);
      for (size_t k=1; k<=nd->getNbDOF(); ++k) {
         if (nd->getDOF(k))
            add(nd->getDOF(k),b(i));
         i++;
      }
   }
}


template<class T_>
void Vect<T_>::Assembly(const Element& el,
                        const T_*      b)
{
   size_t i=0;
   for (size_t n=1; n<=el.getNbNodes(); ++n) {
      Node *nd = el(n);
      for (size_t k=1; k<=nd->getNbDOF(); ++k) {
         if (nd->getDOF(k))
            add(nd->getDOF(k),b[i]);
         i++;
      }
   }
}


template<class T_>
void Vect<T_>::Assembly(const Side&     sd,
                        const Vect<T_>& b)
{
   size_t i=1;
   for (size_t n=1; n<=sd.getNbNodes(); ++n) {
      Node *nd = sd(n);
      for (size_t k=1; k<=nd->getNbDOF(); ++k) {
         if (nd->getDOF(k))
            add(nd->getDOF(k),b(i));
         i++;
      }
   }
}


template<class T_>
void Vect<T_>::Assembly(const Side& sd,
                              T_*   b)
{
   size_t i=0;
   for (size_t n=1; n<=sd.getNbNodes(); ++n) {
      Node *nd = sd(n);
      for (size_t k=1; k<=nd->getNbDOF(); ++k) {
         if (nd->getDOF(k))
            add(nd->getDOF(k),b[i]);
         i++;
      }
   }
}


template <class T_>
void Vect<T_>::getGradient(Vect<T_>& v)
{
   T_ a;
   real_t b;
   Point<T_> aa;
   v.setMesh(*_theMesh,_theMesh->getDim(),ELEMENT_FIELD);
   v.setTime(_time);
   mesh_elements(*_theMesh) {
      try {
         if (The_element.getShape()==LINE) {
            a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
            b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
            v.set(element_label,a/b);
         }
         else if (The_element.getShape()==TRIANGLE) {
            Triang3 t(the_element);
            aa = (*this)(The_element(1)->n())*t.DSh(1) + 
                 (*this)(The_element(2)->n())*t.DSh(2) + 
                 (*this)(The_element(3)->n())*t.DSh(3);
            v.set(element_label,1,aa.x);
            v.set(element_label,2,aa.y);
         }
         else if (The_element.getShape()==TETRAHEDRON) {
	    Tetra4 t(the_element);
            aa = (*this)(The_element(1)->n())*t.DSh(1) + 
                 (*this)(The_element(2)->n())*t.DSh(2) + 
	         (*this)(The_element(3)->n())*t.DSh(3) + 
	         (*this)(The_element(4)->n())*t.DSh(4);
            v.set(element_label,1,aa.x);
            v.set(element_label,2,aa.y);
            v.set(element_label,3,aa.z);
         }
         else
            THROW_RT("getGradient(): This function doesn't work for this element.");
      }
      CATCH("Vect");
   }
}


template <class T_>
void Vect<T_>::getGradient(Vect<Point<T_> >& v)
{
   T_ a;
   real_t b;
   Point<T_> aa;
   v.setMesh(*_theMesh,_theMesh->getDim());
   v.setTime(_time);
   mesh_elements(*_theMesh) {
      try {
      if (The_element.getShape()==LINE) {
         a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
         b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
         v.set(element_label,a/b);
      }
      else if (the_element->getShape()==TRIANGLE) {
         Triang3 t(the_element);
         aa = (*this)(The_element(1)->n())*t.DSh(1) + 
              (*this)(The_element(2)->n())*t.DSh(2) + 
              (*this)(The_element(3)->n())*t.DSh(3);
         v.set(element_label,aa);
      }
      else if (_theMesh->getShape()==TETRAHEDRON) {
         Tetra4 t(the_element);
         aa = (*this)(The_element(1)->n())*t.DSh(1) +
              (*this)(The_element(2)->n())*t.DSh(2) +
              (*this)(The_element(3)->n())*t.DSh(3) +
              (*this)(The_element(4)->n())*t.DSh(4);
         v.set(element_label,aa);
      }
      else
         THROW_RT("getGradient(): This function doesn't work for this element.");
      }
      CATCH_EXIT("Vect");
   }
}


template <class T_>
void Vect<T_>::getCurl(Vect<T_>& v)
{
   T_ a;
   real_t b;
   Point<T_> du, dv, dw;
   v.setMesh(*_theMesh,_theMesh->getDim(),ELEMENT_FIELD);
   v.setTime(_time);
   mesh_elements(*_theMesh) {
      try {
         if (The_element.getShape()==LINE) {
            a = (*this)(The_element(2)) - (*this)(The_element(1));
            b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
            v.set(element_label,a/b);
         }
         else if (The_element.getShape()==TRIANGLE) {
            Triang3 t(the_element);
            du = (*this)(The_element(1)->n())*t.DSh(1) + 
                 (*this)(The_element(2)->n())*t.DSh(2) + 
                 (*this)(The_element(3)->n())*t.DSh(3);
            v.set(element_label,1, du.y);
	    v.set(element_label,2,-du.x);
         }
         else if (The_element.getShape()==TETRAHEDRON) {
            Tetra4 t(the_element);
            du = (*this)(The_element(1)->n(),1)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),1)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),1)*t.DSh(3) + 
                 (*this)(The_element(4)->n(),1)*t.DSh(4);
            dv = (*this)(The_element(1)->n(),2)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),2)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),2)*t.DSh(3) + 
                 (*this)(The_element(4)->n(),2)*t.DSh(4);
            dw = (*this)(The_element(1)->n(),3)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),3)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),3)*t.DSh(3) + 
                 (*this)(The_element(4)->n(),3)*t.DSh(4);
            v.set(element_label,1,dw.y - dv.z);
            v.set(element_label,2,du.z - dw.x);
            v.set(element_label,3,dv.x - du.y);
         }
         else
            THROW_RT("getCurl(): This function doesn't work for this element.");
      }
      CATCH("Vect");
   }
}


template <class T_>
void Vect<T_>::getCurl(Vect<Point<T_> >& v)
{
   T_ a;
   real_t b;
   Point<T_> du, dv, dw;
   v.setMesh(*_theMesh,_theMesh->getDim(),ELEMENT_FIELD);
   v.setTime(_time);
   mesh_elements(*_theMesh) {
      try {
         if (The_element.getShape()==LINE) {
            a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
            b = The_element(2)->getCoord(1) - The_element(1)->getCoord(1);
            v.set(element_label,a/b);
         }
         else if (The_element.getShape()==TRIANGLE) {
            Triang3 t(the_element);
            du = (*this)(The_element(1)->n())*t.DSh(1) + 
                 (*this)(The_element(2)->n())*t.DSh(2) + 
                 (*this)(The_element(3)->n())*t.DSh(3);
            v.set(element_label,1, du.y);
            v.set(element_label,2,-du.x);
         }
         else if (The_element.getShape()==TETRAHEDRON) {
            Tetra4 t(the_element);
            du = (*this)(The_element(1)->n(),1)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),1)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),1)*t.DSh(3) + 
                 (*this)(The_element(4)->n(),1)*t.DSh(4);
            dv = (*this)(The_element(1)->n(),2)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),2)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),2)*t.DSh(3) + 
                 (*this)(The_element(4)->n(),2)*t.DSh(4);
            dw = (*this)(The_element(1)->n(),3)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),3)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),3)*t.DSh(3) + 
                 (*this)(The_element(4)->n(),3)*t.DSh(4);
            v.set(element_label,Point<T_>(dw.y-dv.z,du.z-dw.x,dv.x-du.y));
         }
         else
            THROW_RT("getCurl(): This function doesn't work for this element.");
      }
      CATCH("Vect");
   }
}


template <class T_>
void Vect<T_>::getSCurl(Vect<T_>& v)
{
   try {
      if (_theMesh->getDim()==1 || _theMesh->getDim()==3)
         THROW_RT("getSCurl(): This function is valid for 2-D only.");
   }
   CATCH("Vect");
   Point<T_> du, dv;
   v.setMesh(*_theMesh,1);
   v.setTime(_time);
   mesh_elements(*_theMesh) {
      try {
         if (The_element.getShape()==TRIANGLE) {
            Triang3 t(the_element);
            du = (*this)(The_element(1)->n(),1)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),1)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),1)*t.DSh(3);
            dv = (*this)(The_element(1)->n(),2)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),2)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),2)*t.DSh(3);
            v.set(element_label,dv.x - du.y);
         }
         else
            THROW_RT("getSCurl(): This function doesn't work for this element");
      }
      CATCH("Vect");
   }
}


template <class T_>
void Vect<T_>::getDivergence(Vect<T_>& v)
{
   T_ a;
   real_t b;
   Point<T_> du, dv, dw;
   v.setMesh(*_theMesh,1);
   v.setTime(_time);
   mesh_elements(*_theMesh) {
      try {
         if (The_element.getShape()==LINE) {
            a = (*this)(The_element(2)->n()) - (*this)(The_element(1)->n());
            b = The_element(2)->getX() - The_element(1)->getX();
            v.set(element_label,a/b);
         }
         else if (The_element.getShape()==TRIANGLE) {
            Triang3 t(the_element);
            du = (*this)(The_element(1)->n(),1)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),1)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),1)*t.DSh(3);
            dv = (*this)(The_element(1)->n(),2)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),2)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),2)*t.DSh(3);
            v.set(element_label,du.x + dv.y);
         }
         else if (The_element.getShape()==TETRAHEDRON) {
            Tetra4 t(the_element);
            du = (*this)(The_element(1)->n(),1)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),1)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),1)*t.DSh(3) + 
                 (*this)(The_element(4)->n(),1)*t.DSh(4);
            dv = (*this)(The_element(1)->n(),2)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),2)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),2)*t.DSh(3) + 
                 (*this)(The_element(4)->n(),2)*t.DSh(4);
            dw = (*this)(The_element(1)->n(),3)*t.DSh(1) + 
                 (*this)(The_element(2)->n(),3)*t.DSh(2) + 
                 (*this)(The_element(3)->n(),3)*t.DSh(3) + 
                 (*this)(The_element(4)->n(),3)*t.DSh(4);
            v.set(element_label,du.x + dv.y + dw.z);
         }
         else
            THROW_RT("getDivergence(): This function doesn't work for this element.");
      }
      CATCH("Vect");
   }
}


template<class T_>
real_t Vect<T_>::getAverage(const Element& el,
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

      case PENTA6:
         return OFELI_SIXTH*((*this)(el(1)->n()) +
                             (*this)(el(2)->n()) +
                             (*this)(el(3)->n()) +
                             (*this)(el(4)->n()) +
                             (*this)(el(5)->n()) +
                             (*this)(el(6)->n()));

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
void Vect<T_>::save(string file,
                    int    opt)
{
   saveField(*this,file,opt);
}


template<class T_>
Vect<T_> &Vect<T_>::MultAdd(const Vect<T_>& x,
                            const T_&       a)
{
   for (size_t i=1; i<=size(); ++i)
      add(i,a*x(i));
   return *this;
}


template<class T_>
void Vect<T_>::Axpy(T_              a,
                    const Vect<T_>& x)
{
   for (size_t i=1; i<=size(); i++)
      add(i,a*x(i));
}


template<class T_>
void Vect<T_>::set(size_t i,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
#if defined (USE_EIGEN)
   _v(i-1) = val;
#else
   (*this)[i-1] = val;
#endif
}


template<class T_>
void Vect<T_>::set(size_t i,
                   size_t j,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
#endif
#if defined (USE_EIGEN)
   _v(ijk(i,j)) = val;
#else
   (*this)[ijk(i,j)] = val;
#endif
}


template<class T_>
void Vect<T_>::set(size_t i,
                   size_t j,
                   size_t k,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
   assert(k>0 && k<=_nz);
#endif
#if defined (USE_EIGEN)
   _v(ijk(i,j,k)) = val;
#else
   (*this)[ijk(i,j,k)] = val;
#endif
}


template<class T_>
void Vect<T_>::clear()
{
   for (size_t i=0; i<_size; i++)
#if defined (USE_EIGEN)
      _v[i] = static_cast<T_>(0);
#else
      (*this)[i] = static_cast<T_>(0);
#endif
}


template<class T_>
void Vect<T_>::add(size_t i,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
#if defined (USE_EIGEN)
   _v(i-1) += val;
#else
   (*this)(i) += val;
#endif
}


template<class T_>
void Vect<T_>::add(size_t i,
                   size_t j,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
#endif
#if defined (USE_EIGEN)
   _v(ijk(i,j)) += val;
#else
   (*this)(i,j) += val;
#endif
}


template<class T_>
void Vect<T_>::add(size_t i,
                   size_t j,
                   size_t k,
                   T_     val)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
   assert(k>0 && k<=_nz);
#endif
#if defined (USE_EIGEN)
   _v(ijk(i,j)) = val;
#else
   (*this)[ijk(i,j,k)] += val;
#endif
}


#if defined (USE_EIGEN)
template<class T_>
T_ &Vect<T_>::operator[](size_t i)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
   return _v[i];
}
#endif


#if defined (USE_EIGEN)
template<class T_>
T_ Vect<T_>::operator[](size_t i) const
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
   return _v[i];
}
#endif


template<class T_>
T_ &Vect<T_>::operator()(size_t i)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
#if defined (USE_EIGEN)
   return _v(i-1);
#else
   return (*this)[i-1];
#endif
}


template<class T_>
T_ Vect<T_>::operator()(size_t i) const
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=size());
#endif
#if defined (USE_EIGEN)
   return _v(i-1);
#else
   return (*this)[i-1];
#endif
}


template<class T_>
T_ &Vect<T_>::operator()(size_t i,
                         size_t j)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
#endif
#if defined (USE_EIGEN)
   return _v(ijk(i,j));
#else
   return (*this)[ijk(i,j)];
#endif
}


template<class T_>
T_ Vect<T_>::operator()(size_t i,
                        size_t j) const
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
#endif
#if defined (USE_EIGEN)
   return _v(ijk(i,j));
#else
   return (*this)[ijk(i,j)];
#endif
}


template<class T_>
T_ &Vect<T_>::operator()(size_t i,
                         size_t j,
                         size_t k)
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
   assert(k>0 && k<=_nz);
#endif
#if defined (USE_EIGEN)
   return _v(ijk(i,j,k));
#else
   return (*this)[ijk(i,j,k)];
#endif
}


template<class T_>
T_ Vect<T_>::operator()(size_t i,
                        size_t j,
                        size_t k) const
{
#ifdef _OFELI_RANGE_CHECK
   assert(i>0 && i<=_nx);
   assert(j>0 && j<=_ny);
   assert(k>0 && k<=_nz);
#endif
#if defined (USE_EIGEN)
   return _v(ijk(i,j,k));
#else
   return (*this)[ijk(i,j,k)];
#endif
}


#if defined (USE_EIGEN)
template<class T_>
Vect<T_> &Vect<T_>::operator=(const VectorX& v)
{
   for (size_t i=0; i<_size; i++)
      (*this)[i] = v[i];
   return *this;
}
#endif


template<class T_>
Vect<T_> &Vect<T_>::operator=(const Vect<T_>& v)
{
   _theMesh = v._theMesh;
   _time = v._time;
   _name = v._name;
   _nb_dof = v._nb_dof;
   _nb = v._nb;
   _grid = v._grid;
   _nx = v._nx; _ny = v._ny; _nz = v._nz;
   _size = v._size;
#if defined (USE_EIGEN)
   setSize(_nx,_ny,_nz);
   for (size_t i=0; i<_size; i++)
      (*this)[i] = v[i];
#else
   for (size_t i=0; i<_size; i++)
      (*this)[i] = v[i];
#endif
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator=(const T_& a)
{
   for (size_t i=1; i<=_size; ++i)
      set(i,a);
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator+=(const Vect<T_>& v)
{
#ifdef _OFELI_RANGE_CHECK
   assert(v.size() == _size);
#endif
   for (size_t i=1; i<=_size; ++i)
      add(i,v(i));
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator+=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      add(i,a);
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator-=(const Vect<T_>& v)
{
#ifdef _OFELI_RANGE_CHECK
   assert(v.size() == _size);
#endif
   for (size_t i=1; i<=_size; ++i)
      add(i,-v(i));
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator-=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      add(i,-a);
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator*=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      set(i,(*this)(i)*a);
   return *this;
}


template<class T_>
Vect<T_> &Vect<T_>::operator/=(const T_& a)
{
   for (size_t i=1; i<=_size; i++)
      set(i,(*this)(i)/a);
   return *this;
}


template<class T_>
T_ Vect<T_>::operator,(const Vect<T_>& v) const
{
   T_ p = 0;
   for (size_t i=0; i<_size; i++)
      p += (*this)[i] * v[i];
   return p;
}


template<class T_>
void Vect<T_>::push_back(const T_& v)
{
#if defined (USE_EIGEN)
   (*this)[_nx] = v;
#else
   vector<T_>::push_back(v);
#endif
   _nx++; _size++;
}


template <class T_>
void Vect<T_>::dof_select(size_t          d,
                          vector<size_t>& dof_list)
{
   size_t k, kk, m, j=d, nd=dof_list.size();
   for (k=0; k<nd; k++) {
      kk = size_t(pow(10.,real_t(nd-k-1)));
      m = j/kk;
      dof_list[k] = m;
      j -= m*kk;
   }
   _nb_dof = 0;
   for (k=0; k<nd; k++)
      if (dof_list[k]!=0)
         _nb_dof++;
}


///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////

/** \fn Vect<T_> operator+(const Vect<T_> &x, const Vect<T_> &y)
 *  \brief Operator + (Addition of two instances of class Vect)
 *  \ingroup VectMat
 *  \return  <tt>x + y</tt>
 */
template<class T_>
Vect<T_> operator+(const Vect<T_>& x,
                   const Vect<T_>& y)
{
#ifdef _OFELI_RANGE_CHECK
   assert(x.size() == y.size());
#endif
#if defined (USE_EIGEN)
   return Vect<T_>(Eigen::Matrix<T_,Eigen::Dynamic,1>(x)+Eigen::Matrix<T_,Eigen::Dynamic,1>(y));
#else
   Vect<T_> v(x);
   for (size_t i=0; i<x.size(); ++i)
      v.add(i+1,y[i]);
   return v;
#endif
}
/*
template<class T_>
struct VectAdd {
   VectAdd(Vect<T_> const& u, Vect<T_> const& v) : _u(u), _v(v) { }
   operator Vect<T_>() const {
      Vect<T_> w(_u.size());
      for (size_t i=0; i<_u.size(); i++)
         w[i] = _u[i] + _v[i];
      return w;
   }

   private:
     Vect<T_> const& _u, _v;
};

template<class T_>
VectAdd<T_> operator+(Vect<T_> const& u, Vect<T_> const& v) { return VectAdd<T_>(u,v); }
*/

/** \fn Vect<T_> operator-(const Vect<T_>& x, const Vect<T_>& y)
 *  \brief Operator - (Difference between two vectors of class Vect)
 *  \ingroup VectMat
 *  \return <tt>x - y</tt>
 */
template<class T_>
Vect<T_> operator-(const Vect<T_>& x,
                   const Vect<T_>& y)
{
#ifdef _OFELI_RANGE_CHECK
   assert(x.size()==y.size());
#endif
#if defined (USE_EIGEN)
   return Vect<T_>(Eigen::Matrix<T_,Eigen::Dynamic,1>(x)-Eigen::Matrix<T_,Eigen::Dynamic,1>(y));
#else
   Vect<T_> v(x);
   for (size_t i=0; i<x.size(); ++i)
      v.add(i+1,-y[i]);
   return v;
#endif
}


/** \fn Vect<T_> operator*(const T_ &a, const Vect<T_> &x)
 *  \brief Operator * (Premultiplication of vector by constant)
 *  \ingroup VectMat
 *  \return <tt>a*x</tt>
 */
template<class T_>
Vect<T_> operator*(const T_&       a,
                   const Vect<T_>& x)
{
#if defined (USE_EIGEN)
   return a*x;
#else
   Vect<T_> v(x.size());
   for (size_t i=1; i<=x.size(); ++i)
      v.set(i,a*x(i));
   return v;
#endif
}


/** \fn Vect<T_> operator*(const Vect<T_>& x, const T_& a)
 *  \brief Operator * (Postmultiplication of vector by constant)
 *  \ingroup VectMat
 *  \return <tt>x*a</tt>
 */
template<class T_>
Vect<T_> operator*(const Vect<T_>& x,
                   const T_&       a)
{
#if defined (USE_EIGEN)
   return Vect<T_>(VectorX(x)*a);
#else
   Vect<T_> v(x.size());
   for (size_t i=1; i<=x.size(); ++i)
      v.set(i,a*x(i));
   return v;
#endif
}


/** \fn Vect<T_> operator/(const T_ &a, const Vect<T_> &x)
 *  \brief Operator / (Divide vector entries by constant)
 *  \ingroup VectMat
 *  \return <tt>x/a</tt>
 */
template<class T_>
Vect<T_> operator/(const Vect<T_>& x,
                   const T_&       a)
{
#if defined (USE_EIGEN)
   return Vect<T_>(Matrix<T_,Eigen::Dynamic,1>(x)/a);
#else
   Vect<T_> v(x);
   for (size_t i=1; i<=x.size(); ++i)
      v.set(i,x(i)/a);
   return v;
#endif
}


/** \fn T_ Dot(const Vect<T_> &x, const Vect<T_> &y)
 *  \brief Calculate dot product of two vectors
 *  \ingroup VectMat
 *  \return Dot (inner or scalar) product
 *  Calculate dot (scalar) product of two vectors
 */
template<class T_>
T_ Dot(const Vect<T_>& x,
       const Vect<T_>& y)
{
#ifdef _OFELI_RANGE_CHECK
   assert(x.size() == y.size());
#endif
#if defined (USE_EIGEN)
   return Matrix<T_,Eigen::Dynamic,1>(x).dot(Matrix<T_,Eigen::Dynamic,1>(y));
#else
   T_ s=0;
   for (size_t i=0; i<x.size(); ++i)
      s += x[i]*y[i];
   return s;
#endif
}


/** \fn real_t Discrepancy(Vect<real_t> &x, const Vect<real_t> &y, int n, int type)
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
inline real_t Discrepancy(Vect<real_t>&       x,
                          const Vect<real_t>& y,
                          int                 n,
                          int                 type=1)
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
      for (size_t i=0; i<s; i++) {
         real_t z = std::abs(x[i]-y[i]);
         d += z*z;
      }
      d = std::sqrt(d/s);
   }
   else
      ;
   x = y;
   if (type==1 && old>0.)
      d /= old;
   return d;
}


/** \fn real_t Discrepancy(Vect<complex_t> &x, const Vect<complex_t> &y, int n, int type)
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
inline real_t Discrepancy(Vect<complex_t>&       x,
                          const Vect<complex_t>& y,
                          int                    n,
                          int                    type=1)
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
      for (size_t i=0; i<s; i++) {
         real_t z=std::abs(x[i]-y[i]);
         d += z*z;
      }
      d = std::sqrt(d/s);
   }
   else
      ;
   x = y;
   if (type==1 && old>0.)
      d /= old;
   return d;
}


/** \fn void Modulus(const Vect<complex_t> &x, Vect<real_t> &y)
 *  \brief Calculate modulus of complex vector
 *  @param [in] x %Vector with complex value entries
 *  @param [out] y %Vector containing moduli of entries of <tt>x</tt>
 *  \ingroup VectMat
 */
inline void Modulus(const Vect<complex_t>& x,
                    Vect<real_t>&          y)
{
   y.setSize(x.size());
   for (size_t i=1; i<=x.size(); ++i) {
      complex_t z = x(i);
      y.set(i,std::sqrt(z.real()*z.real()+z.imag()*z.imag()));
   }
}


/** \fn void Real(const Vect<complex_t> &x, Vect<real_t> &y)
 *  \brief Calculate real part of complex vector
 *  @param [in] x %Vector with complex value entries
 *  @param [out] y %Vector containing real parts of entries of <tt>x</tt>
 *  \ingroup VectMat
 */
inline void Real(const Vect<complex_t>& x,
                 Vect<real_t>&          y)
{
   y.setSize(x.size());
   for (size_t i=1; i<=x.size(); ++i)
      y.set(i,x(i).real());
}


/** \fn void Imag(const Vect<complex_t> &x, Vect<real_t> &y)
 *  \brief Calculate imaginary part of complex vector
 *  @param [in] x %Vector with complex value entries
 *  @param [out] y %Vector containing imaginary parts of entries of <tt>x</tt>
 *  \ingroup VectMat
 */
inline void Imag(const Vect<complex_t>& x,
                 Vect<real_t>&          y)
{
   y.setSize(x.size());
   for (size_t i=1; i<=x.size(); ++i)
      y.set(i,x(i).imag());
}


/** \fn istream& operator>>(istream& s, Vect<T_> & a)
 *  Read vector from input stream
 *  \ingroup VectMat
 */
template<class T_>
istream& operator>>(istream&  s,
                    Vect<T_>& v)
{
   T_ z;
   for (size_t i=1; i<=v.size(); ++i) {
      s >> i >> z;
      v.set(i,z);
   }
   return s;
}


/** \fn ostream &operator<<(ostream &s, const Vect<T_> &v)
 *  \brief Output vector in output stream
 *  \ingroup VectMat
 */
template<class T_>
ostream &operator<<(ostream&        s,
                    const Vect<T_>& v)
{
   s.setf(ios::scientific);
   if (v.getNz()>1) {
      for (size_t i=1; i<=v.getNx(); i++) {
         if (v.getNy()==1)
            s << setw(6) << i << "  " << setprecision(8) << setw(18) << v(i) << endl;
         else
            s << "\n[[ i = " << i << " ]]" << endl;
         for (size_t j=1; j<=v.getNy(); j++) {
            s << "\n[ j = " << j << " ]" << endl;
            for (size_t k=1; k<=v.getNz(); k++)
               s << setw(6) << setprecision(8) << setw(18) << v(i,j,k);
            s << endl;
         }
      }
   }
   else {
      if (v.getName() != "#")
         s << v.getName() << " at time = " << v.getTime() << endl << endl;
      s.setf(ios::scientific);
      for (size_t i=1; i<=v.getNx(); i++) {
         s << setw(6) << i << "   ";
         for (size_t j=1; j<=v.getNy(); j++)
            s << setprecision(8) << setw(18) << v(i,j);
         s << endl;
      }
      s << endl;
   }
   return s;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

/*! @} */

#endif
