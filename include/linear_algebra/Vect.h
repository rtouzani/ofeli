/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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
#include <cstdlib>

#include <iostream>
#include <string>
#include <map>
using std::string;

using std::ostream;
using std::istream;
using std::endl;

#include <iomanip>
using std::setw;
using std::ios;
using std::setprecision;

#include "OFELI_Config.h"
#include "io/Fct.h"

#if defined (USE_EIGEN)
#include <Eigen/Dense>
using Eigen::Matrix;
#endif

/** \defgroup VectMat Vector and Matrix
 *  \brief Vector and matrix classes
 */

/*! \file Vect.h
 *  \brief Definition file for class Vect.
 */


//! A namespace to group all library classes, functions, ...
namespace OFELI {

/*! \enum NormType
 * Choose type of vector norm to compute
 */
enum NormType {
   NORM1,     /*!< 1-norm                              */
   WNORM1,    /*!< Weighted 1-norm (Discrete L1-Norm)  */
   NORM2,     /*!< 2-norm                              */
   WNORM2,    /*!< Weighted 2-norm (Discrete L2-Norm)  */
   NORM_MAX   /*!< Max-norm (Infinity norm)            */
};


/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! 
 *  \class Vect
 *  \ingroup VectMat
 *  \brief To handle general purpose vectors.
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
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

class Mesh;
class Element;
class Side;
class Grid;
template<class T_,size_t N_> class LocalVect;
template<class T_,size_t NR_,size_t NC_> class LocalMatrix;
template<class T_> struct Point;

#if defined (USE_PETSC)
template<class T_> class PETScVect;
#endif

template<class T_>
class Vect
#if !defined (USE_EIGEN)
          : public vector<T_>
#endif
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

/** \brief Constructor with a Grid instance
 *  \details The constructed vector has as size the total number of grid nodes
 *  @param [in] g Grid instance
 */
    Vect(Grid& g);

/** \brief Constructor with a mesh instance
 *  @param [in] m Mesh instance
 *  @param [in] dof_type Type of degrees of freedom. To be given among the enumerated
 *  values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>, <tt>SIDE_DOF</tt> or <tt>EDGE_DOF</tt>
 *  (Default: <tt>NODE_DOF</tt>)
 *  @param [in] nb_dof Number of degrees of freedom per node, element or side
 *  If <tt>nb_dof</tt> is set to <tt>0</tt> (default value) the constructor picks this number from
 *  the Mesh instance
 */
    Vect(Mesh&      m,
         DOFSupport dof_type=NODE_DOF,
         int        nb_dof=0);

/** \brief Constructor with a mesh instance giving name and time for vector
 *  @param [in] m Mesh instance
 *  @param [in] dof_type Type of degrees of freedom. To be given among the enumerated
 *  values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>, <tt>SIDE_DOF</tt> or <tt>EDGE_DOF</tt>
 *  @param [in] name Name of the vector
 *  @param [in] nb_dof Number of degrees of freedom per node, element or side
 *  If <tt>nb_dof</tt> is set to <tt>0</tt> the constructor picks this number from the Mesh instance
 *  @param [in] t Time value for the vector [Default <tt>0.0</tt>]
 */
    Vect(Mesh&      m,
         DOFSupport dof_type,
         string     name,
         int        nb_dof=0,
         real_t     t=0.0);

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
 *  which are values of vector. This expression must use the variable \c x as coordinate of
 *  vector.
 *  @warning If the time variable <tt>t</tt> is involved in the expression, the time value
 *  associated to the vector instance must be defined (Default value is 0) either by using
 *  the appropriate constructor or by the member function setTime.
 *  @param [in] x Vector that defines coordinates
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
 *  which are coordinates.
 *  Consider for instance that we want to initialize the Vect instance with the values
 *  v[i] = exp(1+x[i]);
 *  then, we use this member function as follows
 *  v.set("exp("1+x",x);
 */
    void set(const Vect<real_t>& x,
             const string&       exp);

/** \brief Define mesh class to size vector
 *  @param [in] m Mesh instance
 *  @param [in] dof_type Parameter to precise the type of degrees of freedom. To be chosen
 *  among the enumerated values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>,
 *  <tt>SIDE_DOF</tt>, <tt>EDGE_DOF</tt> [Default: <tt>NODE_DOF</tt>]
 *  @param [in] nb_dof Number of degrees of freedom per node, element or side.
 *  If <tt>nb_dof</tt> is set to <tt>0</tt> the constructor picks this number from the Mesh instance
 *  [Default: <tt>0</tt>]
 */
    void setMesh(Mesh&      m,
                 DOFSupport dof_type=NODE_DOF,
                 size_t     nb_dof=0);

/** \brief Define grid class to size vector
 *  @param [in] g Grid instance
 */
    void setGrid(Grid& g);

/// \brief Return vector (global) size
    size_t size() const;

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
    void setDOFType(DOFSupport dof_type);

/** \brief Set Discontinuous Galerkin type vector
 *  \details When the vector is associated to a mesh, this one is sized differently if the
 *  DG method is used.
 *  @param [in] degree Polynomial degree of the DG method [Default: <tt>1</tt>]
 */
    void setDG(int degree=1);

/** \brief Say if vector is constructed for a grid
 *  \details Vectors constructed for grids are defined with the help of a Grid instance
 *  @return true if vector is constructed with a Grid instance
 */
    bool isGrid() const;

/// \brief Return vector number of degrees of freedom
    size_t getNbDOF() const;

/// \brief Return vector number of entities (nodes, elements or sides)
    size_t getNb() const;

/// \brief Return Mesh instance
    Mesh &getMesh() const;

/// \brief Return <tt>true</tt> if vector contains a Mesh pointer, <tt>false</tt> if not
/// \details A Vect instance can be constructed using mesh information 
    bool WithMesh() const;

/** Return DOF type of vector
 *  @return dof_type Type of degrees of freedom. Value among the enumerated
 *  values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>, <tt>SIDE_DOF</tt> or <tt>EDGE_DOF</tt>
 */
    DOFSupport getDOFType() const;

/// \brief Set time value for vector
    void setTime(real_t t);

/// \brief Get time value for vector
    real_t getTime() const;

/// \brief Set name of vector
    void setName(string name);

/// \brief Get name of vector
    string getName() const;

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
    real_t getWNorm1() const;

/** \brief Calculate weighted 2-norm of vector
 *  \details The weighted 2-norm is the 2-Norm of the vector divided by the
 *  square root of its size
 */
    real_t getWNorm2() const;

/// \brief Calculate Min value of vector entries
    T_ getMin() const;

/// \brief Calculate Max value of vector entries
    T_ getMax() const;

/// \brief Return number of grid points in the <tt>x</tt>-direction if grid indexing is set
    size_t getNx() const;

/// \brief Return number of grid points in the <tt>y</tt>-direction if grid indexing is set
    size_t getNy() const;

/// \brief Return number of grid points in the <tt>z</tt>-direction if grid indexing is set
    size_t getNz() const;

/** \brief Assign a given function (given by an interpretable algebraic expression) of indices 
 *  components of vector.
 *  \details This function enable assigning a value to vector entries as function of indices 
 *  @param [in] exp  Regular algebraic expression to assign. It must involve the variables <tt>i</tt>,
 *  <tt>j</tt> and/or <tt>k</tt>.
 */
    void setIJK(const string& exp);

/** \brief Assign a given value to components of vector with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] m Mesh instance
 *  @param [in] code The value is assigned if the node has this code
 *  @param [in] val Value to assign
 *  @param [in] dof Degree of freedom to assign
 */
    void setNodeBC(Mesh&  m,
                   int    code,
                   T_     val,
                   size_t dof);

/** \brief Assign a given value to components of vector with given code
 *  \details Vector components are assumed nodewise. Here all dofs of nodes
 *  with given code will be assigned
 *  @param [in] m Mesh instance
 *  @param [in] code The value is assigned if the node has this code
 *  @param [in] val Value to assign
 */
    void setNodeBC(Mesh& m,
                   int   code,
                   T_    val);

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

/** \brief Assign a given function (given by an interpretable algebraic expression) to 
 *  components of vector with given code.
 *  \details Vector components are assumed nodewise
 *  @param [in] m    Instance of mesh
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned
 */
    void setNodeBC(Mesh&         m,
                   int           code,
                   const string& exp,
                   size_t        dof);

/** \brief Assign a given function (given by an interpretable algebraic expression) to 
 *  components of vector with given code.
 *  \details Vector components are assumed nodewise. Case of 1-DOF problem
 *  @param [in] m    Instance of mesh
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 */
    void setNodeBC(Mesh&         m,
                   int           code,
                   const string& exp);

/** \brief Assign a given function (given by an interpretable algebraic expression) to 
 *  components of vector corresponding to sides with given code.
 *  \details Vector components are assumed nodewise
 *  @param [in] m    Instance of mesh
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned
 */
    void setSideBC(Mesh&         m,
                   int           code,
                   const string& exp,
                   size_t        dof);

/** \brief Assign a given function (given by an interpretable algebraic expression) to 
 *  components of vector corresponding to sides with given code.
 *  \details Vector components are assumed nodewise. Case of 1-DOF problem
 *  @param [in] m    Instance of mesh
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 */
    void setSideBC(Mesh&         m,
                   int           code,
                   const string& exp);

/** \brief Assign a given value to components of vector with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] val  Value to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned [default: <tt>1</tt>]
 */
    void setNodeBC(int    code,
                   T_     val,
                   size_t dof);

/** \brief Assign a given value to components of vector with given code
 *  \details Vector components are assumed nodewise. Concerns 1-DOF problems
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] val  Value to prescribe
 */
    void setNodeBC(int    code,
                   T_     val);

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
                   size_t        dof);

/** \brief Assign a given function (given by an interpretable algebraic expression) to 
 *  components of vector with given code
 *  \details Vector components are assumed nodewise. Concerns 1-DOF problems
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 *  @warning This member function is to be used in the case where a constructor with a Mesh 
 *  has been used
 */
    void setNodeBC(int           code,
                   const string& exp);

/** \brief Assign a given function (given by an interpre<table algebraic expression) to 
 *  components of vector with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 *  @param [in] dof  Degree of Freedom for which the value is assigned
 *  @warning This member function is to be used in the case where a constructor with a Mesh 
 *  has been used
 */
    void setSideBC(int           code,
                   const string& exp,
                   size_t        dof);

/** \brief Assign a given function (given by an interpre<table algebraic expression) to 
 *  components of vector with given code
 *  \details Vector components are assumed nodewise. Case of 1-DOF problem
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] exp  Regular algebraic expression to prescribe
 *  @warning This member function is to be used in the case where a constructor with a Mesh 
 *  has been used
 */
    void setSideBC(int           code,
                   const string& exp);

/** \brief Assign a given value to components of vector with given code
 *  \details Vector components are assumed nodewise
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] val Value to prescribe
 *  @param [in] dof Degree of Freedom for which the value is assigned
 *  @warning This member function is to be used in the case where a constructor with a Mesh 
 *  has been used
 */
    void setSideBC(int    code,
                   T_     val,
                   size_t dof);

/** \brief Assign a given value to components of vector with given code
 *  \details Vector components are assumed nodewise. Concerns 1-DOF problems
 *  @param [in] code Code for which nodes will be assigned prescribed value
 *  @param [in] val Value to prescribe
 *  @warning This member function is to be used in the case where a constructor with a Mesh 
 *  has been used
 */
    void setSideBC(int    code,
                   T_     val);

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
/** \brief Define a regular expression that initializes vector
 *  \details This member function enables defining a regular expression to assign a value
 *  to all vector entries
 *  @param [in] dof Label of dof for which expression is defined [Default: <tt>1</tt>] 
 */
    void setRegex(int dof=1);

/** \brief Returns <tt>true</tt> or <tt>false</tt> if a given dof is defined by 
 *  a regular expression
 *  @param [in] dof Label of dof which information is returned
 */
    bool withRegex(int dof) const;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

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
                  const T_*   b);

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
                      int            type) const;
   
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
    void Axpy(T_              a,
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

/** \brief Add a value to an entry for a 1-index vector
 *  @param [in] i Rank index in vector (starts at <tt>1</tt>)
 *  @param [in] val Value to assign
 */
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
 *  @param [in] vmax Maximal value to assign to the lase entry
 *  @param [in] n Number of points (including extremities)
 *  @remark The vector has a size of \c n. It is sized in this function
 */
    void setUniform(T_     vmin,
                    T_     vmax,
                    size_t n);

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
    const Mesh &getMeshPtr() const;
    
/** \brief Return Dot (scalar) product of two vectors
 *  \details A typical use of this operator is <tt>double a = (v,w)</tt>
 *  where <tt>v</tt> and <tt>w</tt> are 2 instances of <tt>Vect<double></tt>
 *  @param [in] v Vect instance by which the current instance is multiplied
 */
    T_ operator,(const Vect<T_>& v) const;

/** \brief Compute FFT transform of vector
 *  \details This member function computes the FFT (Fast Fourier Transform) of the vector
 *  contained in the instance and returns it
 *  @return Vect<complex<double> > instance containing the FFT
 *  @remark The size of Vect instance must be a power of two and must not exceed
 *  the value of 2^MAX_FFT_SIZE (This value is set in the header "constants.h")
 *  @remark The Vect instance can be either a Vect<double> or Vec<complex<double> >
 */
    Vect<complex_t> getFFT();


/** \brief Compute Inverse FFT transform of vector
 *  \details This member function computes the inverse FFT (Fast Fourier Transform) of the vector
 *  contained in the instance and returns it
 *  @return Vect<complex<double> > instance containing the FFT
 *  @remark The size of Vect instance must be a power of two and must not exceed
 *  the value of 2^MAX_FFT_SIZE (This value is set in the header "constants.h")
 *  @remark The Vect instance can be either a Vect<double> or Vec<complex<double> >
 */
    Vect<complex_t> getInvFFT();

#if defined (USE_EIGEN)
/** \brief Casting operator
 *  @warning This constructor is available only if the library <tt>eigen</tt> is used in 
 *  conjunction with OFELI
 */
    operator VectorX() const;
#endif


 private:

    DOFSupport _dof_type;
    size_t _size, _nx, _ny, _nz, _nb, _nb_dof;
    int    _dg_degree;
    bool   _grid, _with_mesh, _with_regex[10];
    Mesh   *_theMesh;
    string _name, _regex[10];
    real_t _time;
    const vector<string> _var {"x","y","z","t"};
    const vector<string> _var_xit {"x","i","t"};
    const vector<string> _var_ijkt {"i","j","k","t"};
    Fct _theFct;
    void dof_select(size_t d, vector<size_t> &dof_list);
    size_t ijk(size_t i, size_t j)           const { return _ny*(i-1)+j-1; }
    size_t ijk(size_t i, size_t j, size_t k) const { return _ny*_nz*(i-1)+_nz*(j-1)+k-1; }
#if defined (USE_EIGEN)
    VectorX     _v;
#endif

};


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
                   const Vect<T_>& y);

/** \fn Vect<T_> operator-(const Vect<T_>& x, const Vect<T_>& y)
 *  \brief Operator - (Difference between two vectors of class Vect)
 *  \ingroup VectMat
 *  \return <tt>x - y</tt>
 */
template<class T_>
Vect<T_> operator-(const Vect<T_>& x,
                   const Vect<T_>& y);


/** \fn Vect<T_> operator*(const T_ &a, const Vect<T_> &x)
 *  \brief Operator * (Premultiplication of vector by constant)
 *  \ingroup VectMat
 *  \return <tt>a*x</tt>
 */
template<class T_>
Vect<T_> operator*(const T_&       a,
                   const Vect<T_>& x);


/** \fn Vect<T_> operator*(const Vect<T_>& x, const T_& a)
 *  \brief Operator * (Postmultiplication of vector by constant)
 *  \ingroup VectMat
 *  \return <tt>x*a</tt>
 */
template<class T_>
Vect<T_> operator*(const Vect<T_>& x,
                   const T_&       a);


/** \fn Vect<T_> operator/(const Vect<T_> &x, const T_ &a)
 *  \brief Operator / (Divide vector entries by constant)
 *  \ingroup VectMat
 *  \return <tt>x/a</tt>
 */
template<class T_>
Vect<T_> operator/(const Vect<T_>& x,
                   const T_&       a);


/** \fn T_ Dot(const Vect<T_> &x, const Vect<T_> &y)
 *  \brief Calculate dot product of two vectors
 *  \ingroup VectMat
 *  \return Dot (inner or scalar) product
 *  Calculate dot (scalar) product of two vectors
 */
template<class T_>
T_ Dot(const Vect<T_>& x,
       const Vect<T_>& y);


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
real_t Discrepancy(Vect<real_t>&       x,
                   const Vect<real_t>& y,
                   int                 n,
                   int                 type=1);


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
real_t Discrepancy(Vect<complex_t>&       x,
                   const Vect<complex_t>& y,
                   int                    n,
                   int                    type=1);


/** \fn void Modulus(const Vect<complex_t> &x, Vect<real_t> &y)
 *  \brief Calculate modulus of complex vector
 *  @param [in] x %Vector with complex value entries
 *  @param [out] y %Vector containing moduli of entries of <tt>x</tt>
 *  \ingroup VectMat
 */
void Modulus(const Vect<complex_t>& x,
             Vect<real_t>&          y);


/** \fn void Real(const Vect<complex_t> &x, Vect<real_t> &y)
 *  \brief Calculate real part of complex vector
 *  @param [in] x %Vector with complex value entries
 *  @param [out] y %Vector containing real parts of entries of <tt>x</tt>
 *  \ingroup VectMat
 */
void Real(const Vect<complex_t>& x,
          Vect<real_t>&          y);


/** \fn void Imag(const Vect<complex_t> &x, Vect<real_t> &y)
 *  \brief Calculate imaginary part of complex vector
 *  @param [in] x %Vector with complex value entries
 *  @param [out] y %Vector containing imaginary parts of entries of <tt>x</tt>
 *  \ingroup VectMat
 */
void Imag(const Vect<complex_t>& x,
          Vect<real_t>&          y);


/** \fn istream& operator>>(istream& s, Vect<T_> & a)
 *  Read vector from input stream
 *  \ingroup VectMat
 */
template<class T_>
istream& operator>>(istream&  s,
                    Vect<T_>& v);


/** \fn ostream &operator<<(ostream &s, const Vect<T_> &v)
 *  \brief Output vector in output stream
 *  \details Level of vector output depends on the global variable \c Verbosity
 *  <ul>
 *     <li> If Verbosity=0, this function outputs vector size only.
 *     <li> If Verbosity>0, this function outputs vector size, vector name, value
 *          of time, and number of components
 *     <li> If Verbosity>1, this function outputs in addition the first 10 entries in vector 
 *     <li> If Verbosity>2, this function outputs in addition the first 50 entries in vector 
 *     <li> If Verbosity>3, this function outputs in addition the first 100 entries in vector 
 *     <li> If Verbosity>4, this function outputs all vector entries 
 *  </ul>
 *  \ingroup VectMat
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
ostream &operator<<(ostream&        s,
                    const Vect<T_>& v);


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
void OddEven(vector<T_>& x,
             vector<T_>& odd,
             vector<T_>& even);


template<class T_>
void OddEven(Vect<T_>&   x,
             vector<T_>& odd,
             vector<T_>& even);


void fft(vector<complex_t>& x);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
