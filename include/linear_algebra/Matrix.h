/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

                     Definition of abstract class 'Matrix'

  ==============================================================================*/


#ifndef __MATRIX_H
#define __MATRIX_H

#include <iostream>
using std::ostream;

#include <algorithm>

#include "mesh/Mesh.h"
#include "linear_algebra/GraphOfMatrix.h"
#include "linear_algebra/Vect.h"
#include "util/util.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

template<class T_> class DMatrix;

class Mesh;
class Element;
class Side;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
typedef std::pair<size_t,size_t> RC;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/*! \file Matrix.h
 *  \brief Definition file for abstract class Matrix.
 */

/*! \enum MatrixType
 * Choose matrix storage and type
 */
enum MatrixType {
   SKYLINE        = 0x00800000,   /*!< Skyline storage         */
   SPARSE         = 0x01000000,   /*!< Sparse storage          */
   DIAGONAL       = 0x02000000,   /*!< Diagonal storage        */
   TRIDIAGONAL    = 0x04000000,   /*!< Tridiagonal storage     */
   SYMMETRIC      = 0x08000000,   /*!< Symmetric matrix        */
   UNSYMMETRIC    = 0x10000000,   /*!< Unsymmetric matrix      */
   IDENTITY       = 0x12000000    /*!< Identity matrix         */
};

/*! \enum Iteration
 * \brief Choose iterative solver for the linear system
 */
enum Iteration {
   DIRECT_SOLVER    = 0,          /*!< Direct solver           */
   CG_SOLVER        = 1,          /*!< CG Method               */
   CGS_SOLVER       = 2,          /*!< CGS Metod               */
   BICG_SOLVER      = 3,          /*!< BiCG Method             */
   BICG_STAB_SOLVER = 4,          /*!< BiCGStab Method         */
   GMRES_SOLVER     = 5           /*!< GMRes Method            */ 
};

/*! \enum Preconditioner
 * \brief Choose preconditioner for the linear system
 */
enum Preconditioner {
   IDENT_PREC       = 0,          /*!< Identity (No preconditioning)                 */
   DIAG_PREC        = 1,          /*!< Diagonal preconditioner                       */
   DILU_PREC        = 2,          /*!< ILU (Incomplete factorization) preconditioner */
   ILU_PREC         = 3,          /*!< DILU (Diagonal Incomplete factorization) preconditioner */
   SSOR_PREC        = 4           /*!< SSOR preconditioner                           */
};

template <class T_> class Vect;
class Mesh;
class Element;
class Side;

 
/*! \class Matrix
 *  \brief Virtual class to handle matrices for all storage formats.
 * \details
 *
 * \tparam <T_> Data type (real_t, float, complex<real_t>, ...)
 */
template<class T_>
class Matrix
{

 public:

/// \brief Default constructor.
/// \details Initializes a zero-size matrix.
    Matrix();

/// \brief Copy Constructor
    Matrix(const Matrix<T_> &m);

/// \brief Destructor
    virtual ~Matrix() { }

/// \brief Return number of rows.
    size_t getNbRows() const { return _nb_rows; }

/// \brief Return number of columns.
    size_t getNbColumns() const { return _nb_cols; }

/// \brief Set Penalty Parameter (For boundary condition prescription).
    void setPenal(real_t p) { _penal = p; }

/// \brief Set the matrix as diagonal
    void setDiagonal();

/// \brief Return <tt>k</tt>-th diagonal entry of matrix.
/// \details First entry is given by \b getDiag(1).
    T_ getDiag(size_t k) const { return _diag[k-1]; }

/// \brief Return matrix dimension (Number of rows and columns).
    size_t size() const { return _size; }

/// \brief Multiply matrix by vector <tt>x</tt> and add to <tt>y</tt>
    virtual void MultAdd(const Vect<T_>& x,
                               Vect<T_>& y) const = 0;

/// \brief Multiply matrix by vector <tt>a*x</tt> and add to <tt>y</tt>
    virtual void MultAdd(      T_        a,
                         const Vect<T_>& x,
                               Vect<T_>& y) const = 0;

/// \brief Multiply matrix by vector <tt>x</tt> and save in <tt>y</tt>
    virtual void Mult(const Vect<T_>& x,
                            Vect<T_>& y) const = 0;

/// \brief Multiply transpose of matrix by vector <tt>x</tt> and save in <tt>y</tt>
    virtual void TMult(const Vect<T_>& v,
                             Vect<T_>& w) const = 0;

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] x %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    virtual void Axpy(      T_          a,
                      const Matrix<T_>* x) = 0;

/// \brief Initialize matrix storage in the case where only diagonal terms are stored.
/// \details This member function is to be used for explicit time integration schemes
    void setDiagonal(Mesh& mesh);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    virtual void setMesh(Mesh&  mesh,
                         size_t dof=0)=0;
    void init_set_mesh(Mesh&  mesh,
                       size_t dof=0);
    void init_set_mesh(Mesh&  mesh,
                       size_t type,
                       size_t dof1,
                       size_t dof2);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    virtual void setMesh(size_t dof,
                         Mesh&  mesh,
                         int    code=0)=0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    virtual void setMesh(size_t dof,
                         size_t nb_eq,
                         Mesh&  mesh)=0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Assembly of element matrix into global matrix.
 *  \details Case where element matrix is given by a C-array.
 *  @param [in] el Pointer to element instance
 *  @param [in] a Element matrix as a C-array
 */
    void Assembly(const Element& el,
                        T_*      a);

/** \brief Assembly of element matrix into global matrix.
 *  \details Case where element matrix is given by a DMatrix instance.
 *  @param [in] el Pointer to element instance
 *  @param [in] a Element matrix as a DMatrix instance
 */
    void Assembly(const Element&     el,
                  const DMatrix<T_>& a);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Assembly of element matrix into global matrix for a Discontinuous Galerkin approximation
 *  \details Case where element matrix is given by a C-array.
 *  @param [in] el Pointer to element instance
 *  @param [in] a Element matrix as a C-array
 */
    void DGAssembly(const Element& el,
                          T_*      a);

    void DGAssembly(const Element&                                               el,
                    const LocalMatrix<T_,MAX_NB_ELEMENT_DOF,MAX_NB_ELEMENT_DOF>& a);
    void DGAssembly(const Side&                                                  sd,
                    const LocalMatrix<T_,MAX_NB_SIDE_DOF,MAX_NB_SIDE_DOF>&       a);

/** \brief Assembly of element matrix into global matrix for a Discontinuous Galerkin approximation
 *  \details Case where element matrix is given by a DMatrix instance.
 *  @param [in] el Pointer to element instance
 *  @param [in] a Element matrix as a DMatrix instance
 */
    void DGAssembly(const Element&     el,
                    const DMatrix<T_>& a);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Assembly of side matrix into global matrix.
 *  \details Case where side matrix is given by a C-array.
 *  @param [in] sd Pointer to side instance
 *  @param [in] a Side matrix as a C-array instance
 */
    void Assembly(const Side& sd, 
                              T_*   a);

/** \brief Assembly of side matrix into global matrix.
 *  \details Case where side matrix is given by a DMatrix instance.
 *  @param [in] sd Pointer to side instance
 *  @param [in] a Side matrix as a DMatrix instance
 */
    void Assembly(const Side&        sd,
                  const DMatrix<T_>& a);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Prescribe(      Mesh& mesh,
                         Vect<T_>& b,
                   const Vect<T_>& u,
                         int       flag=0);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Impose by a penalty method an essential boundary condition, using the Mesh instance
 *  provided by the constructor
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  The penalty parameter is defined by default equal to 1.e20.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in,out] b Vect instance that contains right-hand side.
 *  @param [in] u Vect instance that contains imposed valued at DOFs where they are to be imposed.
 *  @param [in] flag Parameter to determine whether only the right-hand side is to
 *  be modified (<tt>dof>0</tt>)\n
 *  or both matrix and right-hand side (<tt>dof=0</tt>, default value).
 */
    void Prescribe(      Vect<T_>& b,
                   const Vect<T_>& u,
                         int       flag=0)
       { Prescribe(*_theMesh,b,u,flag); }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Prescribe(      int       dof,
                         int       code,
                         Mesh&     mesh,
                         Vect<T_>& b,
                   const Vect<T_>& u,
                         int       flag=0);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Impose by a penalty method an essential boundary condition to a given
 *  degree of freedom for a given code
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  The penalty parameter is defined by default equal to 1.e20.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in] dof Degree of freedom for which a boundary condition is to be enforced
 *  @param [in] code Code for which a boundary condition is to be enforced
 *  @param [in,out] b Vect instance that contains right-hand side.
 *  @param [in] u Vect instance that contains imposed valued at DOFs where they are to be imposed.
 *  @param [in] flag Parameter to determine whether only the right-hand side is to be modified\n
 *  (<tt>dof>0</tt>) or both matrix and right-hand side (<tt>dof=0</tt>, default value).
 */
    void Prescribe(      int       dof,
                         int       code,
                         Vect<T_>& b,
                   const Vect<T_>& u,
                         int       flag=0)
       { Prescribe(dof,code,*_theMesh,b,u,flag); }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Prescribe(Mesh&     mesh,
                   Vect<T_>& b,
                   int       flag=0);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Impose by a penalty method a homegeneous (=0) essential boundary condition.
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  The penalty parameter is defined by default equal to 1.e20.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in,out] b Vect instance that contains right-hand side.
 *  @param [in] flag Parameter to determine whether only the right-hand side is to
 *  be modified (<tt>dof>0</tt>)\n
 *  or both matrix and right-hand side (<tt>dof=0</tt>, default value).
 */
    void Prescribe(Vect<T_>& b,
                   int       flag=0)
       { Prescribe(*_theMesh,b,flag); }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Prescribe(      size_t    dof,
                         Mesh&     mesh,
                         Vect<T_>& b,
                   const Vect<T_>& u,
                         int       flag=0);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Impose by a penalty method an essential boundary condition when only one DOF is treated.
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  This gunction is to be used if only one DOF per node is treated in the linear system.
 *  The penalty parameter is by default equal to 1.e20.
 *  It can be modified by member function setPenal.
 *  @param [in] dof Label of the concerned degree of freedom (DOF).
 *  @param [in,out] b Vect instance that contains right-hand side.
 *  @param [in] u Vect instance that conatins imposed valued at DOFs where they are to be imposed.
 *  @param [in] flag Parameter to determine whether only the right-hand side is to
 *  be modified (<tt>dof>0</tt>)\n
 *  or both matrix and right-hand side (<tt>dof=0</tt>, default value).
 */
    void Prescribe(      size_t    dof,
                         Vect<T_>& b,
                   const Vect<T_>& u,
                         int       flag=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Prescribe1(      Mesh&     mesh,
                          Vect<T_>& b,
                    const Vect<T_>& u,
                          int       flag=0);
    void Prescribe1(      Vect<T_>& b,
                    const Vect<T_>& u,
                          int       flag=0)
       { Prescribe1(*_theMesh,b,u,flag); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Prescribe(Mesh& mesh);
    void Prescribe() { Prescribe(*_theMesh); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void PrescribeSide(Mesh& mesh);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Impose by a penalty method an essential boundary condition when
 *  DOFs are supported by sides.
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  The penalty parameter is defined by default equal to 1.e20.
 *  It can be modified by member function \b setPenal(..).
 */
    void PrescribeSide() { PrescribeSide(*_theMesh); }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Constraint(const Mesh& mesh) { Prescribe(mesh); }
    void Constraint() { Prescribe(); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Add <tt>val</tt> to entry <tt>(i,j)</tt>.
    virtual void add(      size_t  i,
                           size_t  j,
                     const T_&     val) = 0;

/// \brief Factorize matrix. Available only if the storage class enables it.
    virtual int Factor() = 0;

/** \brief Solve the linear system by a direct method.
 *  \details This is available only if the storage class enables it
 *  and if matrix has been primarily factorized (See \b isFactorized).
 */
    virtual int solve(Vect<T_>& b) = 0;

/** \brief Solve system with factorized matrix (forward and back substitution).
 *  @param [in] b Vect instance that contains right-hand side
 *  @param [out] x Vect instance that contains solution
 *  @return - \a 
 *    <ul>
 *       <li><tt>0</tt> if solution was normally performed
 *       <li><tt>n</tt> if the <tt>n</tt>-th pivot is null\n
 *                      Solution is performed only is factorization has previouly been invoked.
 *    </uL>
 */
    int solve(const Vect<T_>& b,
                    Vect<T_>& x);

/** \brief Factorize matrix and solve the linear system.
 *  \details This is available only if the storage cass enables it.
 *  @param [in,out] b Vect instance that contains right-hand side on input and
 *  solution on output
 */
    int FactorAndSolve(Vect<T_>& b);

/** \brief Factorize matrix and solve the linear system.
 *  \details This is available only if the storage class enables it.
 *  @param [in] b Vect instance that contains right-hand side
 *  @param [out] x Vect instance that contains solution
 *  @return
 *    <ul>
 *     <li><tt>0</tt> if solution was normally performed
 *     <li><tt>n</tt> if the n-th pivot is nul
 *    </ul>
 */
    int FactorAndSolve(const Vect<T_>& b,
                             Vect<T_>& x);

/// \brief Return number of stored terms in matrix.
    size_t getLength() const { return _length; }

/// \brief Say if matrix is diagonal or not
    int isDiagonal() const { return _is_diagonal; }

/// \brief Say if matrix is factorized or not.
/// \details If the matrix was not factorized, the class does not allow
/// solving by a direct solver.
    int isFactorized() const { return _fact; }

/// \brief Return Column index for column <tt>i</tt> (See the description for class SpMatrix).
    virtual size_t getColInd(size_t i) const;

/// \brief Return Row pointer for row <tt>i</tt> (See the description for class SpMatrix).
    virtual size_t getRowPtr(size_t i) const;

/** \brief Assign a value to an entry of the matrix
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] val Value to assign
 */
    virtual void set(      size_t  i,
                           size_t  j,
                     const T_&     val) = 0;

/** \brief Operator () (Non constant version).
 *  \details Returns the <tt>(i,j)</tt> entry of the matrix.
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    virtual T_ & operator()(size_t i,
                            size_t j) = 0;

/** \brief Operator () (Non constant version).
 *  \details Returns the <tt>(i,j)</tt> entry of the matrix.
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    virtual T_ operator()(size_t i,
                          size_t j) const = 0;

/** \brief Operator () with one argument (Constant version).
 *  \details Returns <tt>i</tt>-th position in the array storing matrix entries.
 *  The first entry is at location 1.
 *  Entries are stored row by row.
 *  @param [in] i entry index
 */
    T_ operator()(size_t i) const;

/** \brief Operator () with one argument (Non Constant version).
 *  \details Returns <tt>i</tt>-th position in the array storing matrix entries.
 *  The first entry is at location 1.
 *  Entries are stored row by row.
 *  @param [in] i entry index
 */
    T_ & operator()(size_t i);

/** \brief Operator [] (Non constant version).
 *  \details Returns <tt>k</tt>-th stored element in matrix
 *  Index <tt>k</tt> starts at <tt>0</tt>.
 */
    T_ & operator[](size_t k);

/** \brief Operator [] (Constant version).
 *  \details Returns <tt>k</tt>-th stored element in matrix
 *  Index <tt>k</tt> starts at <tt>0</tt>.
 */
    T_ operator[](size_t k) const;

/// \brief Operator =.
/// \details Copy matrix <tt>m</tt> to current matrix instance.
    Matrix & operator=(Matrix<T_>& m);

/// \brief Operator +=.
/// \details Add matrix <tt>m</tt> to current matrix instance.
    Matrix & operator+=(const Matrix<T_>& m);

/// \brief Operator -=.
/// \details Subtract matrix <tt>m</tt> from current matrix instance.
    Matrix & operator-=(const Matrix<T_>& m);

/// \brief Operator =.
/// \details Assign constant value <tt>x</tt> to all matrix entries.
    Matrix & operator=(const T_ &x);

/// \brief Operator *=.
/// \details Premultiply matrix entries by constant value <tt>x</tt>
    Matrix & operator*=(const T_& x);

/// \brief Operator +=.
/// \details Add constant value <tt>x</tt> to all matrix entries.
    Matrix & operator+=(const T_& x);

/// \brief Operator -=.
/// \details Subtract constant value <tt>x</tt> from all matrix entries.
    Matrix & operator-=(const T_& x);

/// \brief Return entry <tt>(i,j)</tt> of matrix if this one is stored, <tt>0</tt> else
    virtual T_ get(size_t i,
                   size_t j) const = 0;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    vector<T_>    _a, _aU, _diag;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   size_t         _dof_type, _nb_rows, _nb_cols, _size, _dof, _length;
   T_             _zero, _temp;
   vector<size_t> _row_ptr, _col_ind, _ch;
   real_t         _penal;
   int            _fact, _set_nodes, _set_elements, _set_sides, _is_diagonal;
   Mesh           *_theMesh;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};


///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////


#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_>
Matrix<T_>::Matrix()
{
   _penal = 1.e20;
   _zero = 0;
   _nb_rows = _nb_cols = _size = _length = 0;
   _is_diagonal = false;
}


template<class T_>
Matrix<T_>::Matrix(const Matrix<T_>& m)
{
   _zero = 0;
   _size = m._size;
   _nb_rows = m._nb_rows;
   _nb_cols = m._nb_cols;
   _length = m._length;
   _penal = m._penal;
   _ch.resize(_size);
   _diag.setSize(_size);
   _ch = m._ch;
   _diag = m._diag;
   _theMesh = m._theMesh;
   _is_diagonal = m._is_diagonal;
}


template<class T_>
void Matrix<T_>::setDiagonal(Mesh& mesh)
{
   init_set_mesh(mesh);
   _size = _theMesh->getNbEq();
   _ch.resize(_size);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] = i+1;
   _a.resize(_size);
   Clear(_a);
   _fact = false;
   _dof = 0;
   _length = _nb_rows = _nb_cols = _size;
   _is_diagonal = true;
}


template<class T_>
void Matrix<T_>::setDiagonal()
{
   _size = _theMesh->getNbEq();
   _ch.resize(_size);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] = i+1;
   _a.resize(_size);
   Clear(_a);
   _fact = false;
   _dof = 0;
   _length = _nb_rows = _nb_cols = _size;
   _is_diagonal = true;
}


template<class T_>
void Matrix<T_>::init_set_mesh(Mesh&  mesh,
                               size_t dof)
{
   _theMesh = &mesh;
   _zero = T_(0);
   _dof_type = 0;
   if (_theMesh->NodesAreDOF())
      _dof_type = NODE_DOF;
   else if (_theMesh->SidesAreDOF())
      _dof_type = SIDE_DOF;
   else if (_theMesh->ElementsAreDOF())
      _dof_type = ELEMENT_DOF;
   _dof = dof;
   _size = _nb_rows = _nb_cols = _theMesh->getNbEq();
   if (_dof_type==NODE_DOF)
      if (_dof)
         _size = _nb_rows = _nb_cols = _theMesh->getNbNodes();
      else
         _size = _nb_rows = _nb_cols = _theMesh->getNbEq();
   else if (_dof_type==SIDE_DOF)
      if (_dof)
         _size = _nb_rows = _nb_cols = _theMesh->getNbSides();
      else
         _size = _nb_rows = _nb_cols = _theMesh->getNbEq();
   else if (_dof_type==ELEMENT_DOF)
      _size = _nb_rows = _nb_cols = _theMesh->getNbElements();
   else;
}


template<class T_>
void Matrix<T_>::Assembly(const Element& el,
                                T_*      a)
{
   size_t kk=0;
   if (_is_diagonal) {
      for (size_t i=1; i<=el.getNbNodes(); ++i) {
         Node *nd=el(i);
         size_t nb_dof = nd->getNbDOF();
         for (size_t k=1; k<=nb_dof; ++k) {
            size_t n=nb_dof*(nd->n()-1) + k;
            add(n,n,a[kk]);
            kk += nb_dof*el.getNbNodes() + 1;
	 }
      }
      return;
   }
   for (size_t i=1; i<=el.getNbNodes(); ++i) {
      Node *nd1=el(i);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k) {
         for (size_t j=1; j<=el.getNbNodes(); ++j) {
            Node *nd2=el(j);
            for (size_t l=1; l<=nd2->getNbDOF(); ++l) {
               if (nd1->getDOF(k)!=0 && nd2->getDOF(l)!=0)
                  add(nd1->getDOF(k),nd2->getDOF(l),a[kk]);
               kk++;
            }
         }
      }
   }
}


template<class T_>
void Matrix<T_>::Assembly(const Element&     el,
                          const DMatrix<T_>& a)
{
   size_t i = 1;
   for (size_t in=1; in<=el.getNbNodes(); ++in) {
      Node *nd1=el(in);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k) {
         size_t j=1;
         for (size_t jn=1; jn<=el.getNbNodes(); ++jn) {
            Node *nd2=el(jn);
            for (size_t l=1; l<=nd2->getNbDOF(); ++l) {
               if (nd1->getDOF(k)!=0 && nd2->getDOF(l)!=0)
                  add(nd1->getDOF(k),nd2->getDOF(l),a(i,j));
               j++;
            }
         }
         i++;
      }
   }
}


template<class T_>
void Matrix<T_>::DGAssembly(const Element& el,
                                  T_*      a)
{
   size_t kk=0;
   if (_is_diagonal) {
      for (size_t i=1; i<=el.getNbDOF(); ++i) {
         Node *nd=el(i);
         size_t nb_dof = nd->getNbDOF();
         for (size_t k=1; k<=nb_dof; ++k) {
            size_t n=nb_dof*(nd->n()-1) + k;
            add(n,n,a[kk]);
            kk += nb_dof*el.getNbNodes() + 1;
         }
      }
      return;
   }

   for (size_t i=1; i<=el.getNbNodes(); ++i) {
      Node *nd1 = el(i);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k) {
         for (size_t j=1; j<=el.getNbNodes(); ++j) {
            Node *nd2=el(j);
            for (size_t l=1; l<=nd2->getNbDOF(); ++l) {
               if (nd1->getDOF(k)!=0 && nd2->getDOF(l)!=0)
                  add(nd1->getDOF(k),nd2->getDOF(l),a[kk]);
               kk++;
            }
         }
      }
   }
}


template<class T_>
void Matrix<T_>::DGAssembly(const Element&                                               el,
                            const LocalMatrix<T_,MAX_NB_ELEMENT_DOF,MAX_NB_ELEMENT_DOF>& a)
{
   for (size_t i=1; i<=el.getNbDOF(); ++i) {
      for (size_t j=1; j<=el.getNbDOF(); ++j) {
         if (el.getDOF(i)!=0 && el.getDOF(j)!=0)
            add(el.getDOF(i),el.getDOF(j),a(i,j));
      }
   }
}


template<class T_>
void Matrix<T_>::DGAssembly(const Side&                                            sd,
                            const LocalMatrix<T_,MAX_NB_SIDE_DOF,MAX_NB_SIDE_DOF>& a)
{
   for (size_t i=1; i<=sd.getNbDOF(); ++i) {
      for (size_t j=1; j<=sd.getNbDOF(); ++j) {
         if (sd.getDOF(i)!=0 && sd.getDOF(j)!=0)
            add(sd.getDOF(i),sd.getDOF(j),a(i,j));
      }
   }
}


template<class T_>
void Matrix<T_>::Assembly(const Side& sd,
                                      T_*   a)
{
   size_t kk = 0;
   for (size_t in=1; in<=sd.getNbNodes(); ++in) {
      Node *nd1=sd(in);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k) {
         for (size_t jn=1; jn<=sd.getNbNodes(); ++jn) {
            Node *nd2=sd(jn);
            for (size_t l=1; l<=nd2->getNbDOF(); ++l) {
               if (nd1->getDOF(k)!=0 && nd2->getDOF(l)!=0)
                  add(nd1->getDOF(k),nd2->getDOF(l),a[kk]);
               kk++;
            }
         }
      }
   }
}


template<class T_>
void Matrix<T_>::Assembly(const Side&        sd,
                          const DMatrix<T_>& a)
{
   size_t i=1;
   for (size_t in=1; in<=sd.getNbNodes(); ++in) {
      Node *nd1=sd(in);
      for (size_t k=1; k<=nd1->getNbDOF(); ++k) {
         size_t j=1;
         for (size_t jn=1; jn<=sd.getNbNodes(); ++jn) {
            Node *nd2=sd(jn);
            for (size_t l=1; l<=nd2->getNbDOF(); ++l) {
               if (nd1->getDOF(k)!=0 && nd2->getDOF(l)!=0)
                  add(nd1->getDOF(k),nd2->getDOF(l),a(i,j));
               j++;
            }
         }
         i++;
      }
   }
}


template<class T_>
void Matrix<T_>::Prescribe(      Mesh&     mesh,
                                 Vect<T_>& b,
                           const Vect<T_>& u,
                                 int       flag)
{
   MeshNodes(mesh) {
      for (size_t i=1; i<=theNode->getNbDOF(); ++i) {
         if (theNode->getCode(i)>0) {
            size_t k=theNode->getDOF(i);
            if (flag==0) {
               _diag[k-1] = get(k,k)*_penal;
               set(k,k,_diag[k-1]);
            }
            b.set(k,u(k)*_diag[k-1]);
         }
      }
   }
}


template<class T_>
void Matrix<T_>::Prescribe(      int       dof,
                                 int       code,
                                 Mesh&     mesh,
                                 Vect<T_>& b,
                           const Vect<T_>& u,
                                 int       flag)
{
   MeshNodes(mesh) {
      if (theNode->getCode(dof)==code) {
         size_t k=theNode->getDOF(dof);
         if (flag==0) {
            _diag[k-1] = get(k,k)*_penal;
            set(k,k,_diag[k-1]);
         }
         b.set(k,u(k)*_diag[k-1]);
      }
   }
}


template<class T_>
void Matrix<T_>::Prescribe(Mesh&     mesh,
                           Vect<T_>& b,
                           int       flag)
{
   MeshNodes(mesh) {
      for (size_t j=1; j<=theNode->getNbDOF(); ++j)
         if (theNode->getCode(j)>0) {
            size_t k=theNode->getDOF(j);
            if (!flag) {
               _diag[k-1] = get(k,k)*_penal;
               set(k,k,_diag[k-1]);
            }
            b.set(k,0);
         }
   }
}


template<class T_>
void Matrix<T_>::Prescribe(      size_t    dof,
                                 Mesh&     mesh,
                                 Vect<T_>& b,
                           const Vect<T_>& u,
                                 int       flag)
{
   mesh_nodes(mesh) {
      if (The_node.getCode(dof)>0) {
         size_t k=node_label;
         if (!flag) {
            _diag[k-1] = get(k,k)*_penal;
            set(k,k,_diag[k-1]);
         }
         b.set(k,u(k)*_diag[k-1]);
      }
   }
}


template<class T_>
void Matrix<T_>::Prescribe1(      Mesh&     mesh,
                                  Vect<T_>& b,
                            const Vect<T_>& u,
                                  int       flag)
{
   mesh_nodes(mesh) {
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
         if (The_node.getCode(i)>0) {
            size_t k=The_node.getDOF(i);
            if (!flag)
               add(k,k,_penal);
            b.set(k,u(k)*_penal);
         }
      }
   }
}


template<class T_>
void Matrix<T_>::Prescribe(Mesh& mesh)
{
   mesh_nodes(mesh) {
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
         if (The_node.getCode(i)>0) {
            size_t k=The_node.getDOF(i);
            set(k,k,get(k,k)*_penal);
         }
      }
   }
}


template<class T_>
void Matrix<T_>::PrescribeSide(Mesh& mesh)
{
   mesh_sides(mesh) {
      for (size_t i=1; i<=The_side.getNbDOF(); i++) {
         if (The_side.getCode(i)>0) {
            size_t k=The_side.getDOF(i);
            set(k,k,get(k,k)*_penal);
         }
      }
   }
}


template<class T_>
int Matrix<T_>::solve(const Vect<T_>& b,
                            Vect<T_>& x)
{
   x = b;
   return solve(x);
}


template<class T_>
int Matrix<T_>::FactorAndSolve(Vect<T_>& b)
{
   Factor();
   return solve(b);
}


template<class T_>
int Matrix<T_>::FactorAndSolve(const Vect<T_>& b,
                               Vect<T_>&       x)
{
   int ret = Factor();
   solve(b,x);
   return ret;
}


template<class T_>
size_t Matrix<T_>::getColInd(size_t i) const
{
   i = 0;
   return 0;
}


template<class T_>
size_t Matrix<T_>::getRowPtr(size_t i) const
{
   i = 0;
   return 0;
}


template<class T_>
Matrix<T_> & Matrix<T_>::operator=(Matrix<T_>& m)
{
   _a = m._a;
   return *this;
}


template<class T_>
Matrix<T_> & Matrix<T_>::operator+=(const Matrix<T_>& m)
{
   for (size_t i=1; i<=_length; ++i)
      _a.add(i,m._a[i-1]);
   return *this;
}


template<class T_>
Matrix<T_> & Matrix<T_>::operator-=(const Matrix<T_>& m)
{
   for (size_t i=0; i<_length; ++i)
      _a[i] -= m._a[i];
   return *this;
}


template<class T_>
Matrix<T_> & Matrix<T_>::operator=(const T_& x)
{
   for (size_t i=0; i<_length; ++i)
      _a[i] = x;
   return *this;
}


template<class T_>
Matrix<T_> & Matrix<T_>::operator*=(const T_& x)
{
   for (size_t i=0; i<_length; ++i)
      _a[i] *= x;
   return *this;
}


template<class T_>
Matrix<T_> & Matrix<T_>::operator+=(const T_& x)
{
   for (size_t i=0; i<_length; ++i)
      _a[i] += x;
   return *this;
}


template<class T_>
Matrix<T_> & Matrix<T_>::operator-=(const T_& x)
{
   for (size_t i=0; i<_length; ++i)
      _a[i] = -x;
   return *this;
}


template<class T_>
T_ Matrix<T_>::operator()(size_t i) const
{
   return _a[i-1];
}


template<class T_>
T_ & Matrix<T_>::operator()(size_t i)
{
   return _a[i-1];
}


template<class T_>
T_ & Matrix<T_>::operator[](size_t k)
{
   return _a[k];
}


template<class T_>
T_ Matrix<T_>::operator[](size_t k) const
{
   return _a[k];
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
