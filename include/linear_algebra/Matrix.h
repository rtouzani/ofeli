/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani
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
#include "OFELI_Config.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"
using std::to_string;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

template<class T_,size_t NR_,size_t NC_> class LocalMatrix;
 
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
   DENSE          =   1,   /*!< Dense storage           */
   SKYLINE        =   2,   /*!< Skyline storage         */
   SPARSE         =   4,   /*!< Sparse storage          */
   DIAGONAL       =   8,   /*!< Diagonal storage        */
   TRIDIAGONAL    =  16,   /*!< Tridiagonal storage     */
   BAND           =  32,   /*!< Band storage            */
   SYMMETRIC      =  64,   /*!< Symmetric matrix        */
   UNSYMMETRIC    = 128,   /*!< Unsymmetric matrix      */
   IDENTITY       = 256    /*!< Identity matrix         */
};

/*! \enum Iteration
 * \brief Choose iterative solver for the linear system
 */
enum Iteration {
   DIRECT_SOLVER    = 0,   /*!< Direct solver           */
   CG_SOLVER        = 1,   /*!< CG Method               */
   CGS_SOLVER       = 2,   /*!< CGS Metod               */
   BICG_SOLVER      = 3,   /*!< BiCG Method             */
   BICG_STAB_SOLVER = 4,   /*!< BiCGStab Method         */
   GMRES_SOLVER     = 5    /*!< GMRes Method            */ 
};

/*! \enum Preconditioner
 * \brief Choose preconditioner for the linear system
 */
enum Preconditioner {
   IDENT_PREC       = 0,   /*!< Identity (No preconditioning)                           */
   DIAG_PREC        = 1,   /*!< Diagonal preconditioner                                 */
   DILU_PREC        = 2,   /*!< ILU (Incomplete factorization) preconditioner           */
   ILU_PREC         = 3,   /*!< DILU (Diagonal Incomplete factorization) preconditioner */
   SSOR_PREC        = 4    /*!< SSOR preconditioner                                     */
};


//template<class T_> class DMatrix;
 
/*! \class Matrix
 *  \brief Virtual class to handle matrices for all storage formats.
 * \details
 * This class enables storing and manipulating dense matrices.
 * The template parameter is the type of matrix entries.
 * Any matrix entry can be accessed by the () operator: For instance,
 * if \c A is an instance of this class, \c A(i,j) stands for the entry
 * at the i-th row and j-th column, \c i and \c j starting from 1.
 * Entries of \c A can be assigned a value by the same operator.
 *
 * \tparam <T_> Data type (real_t, float, complex<real_t>, ...)
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
class Matrix
{

 public:

/// \brief Default constructor
/// \details Initializes a zero-size matrix.
  Matrix();

/// \brief Copy Constructor
    Matrix(const Matrix<T_> &m);

/// \brief Destructor
    virtual ~Matrix();

/** \brief Set matrix to 0 and reset factorization parameter
 *  @warning This function must be used if after a factorization, the matrix has
 *  been modified
 */
    virtual void reset();

/// \brief Return number of rows.
    size_t getNbRows() const;

/// \brief Return number of columns.
    size_t getNbColumns() const;

/// \brief Set Penalty Parameter (For boundary condition prescription).
    void setPenal(real_t p);

/// \brief Set the matrix as diagonal
    void setDiagonal();

/// \brief Return <tt>k</tt>-th diagonal entry of matrix.
/// \details First entry is given by \b getDiag(1).
    T_ getDiag(size_t k) const;

/// \brief Return matrix dimension (Number of rows and columns).
    size_t size() const;

/// \brief Multiply matrix by vector <tt>x</tt> and add to <tt>y</tt>
    virtual void MultAdd(const Vect<T_>& x,
                         Vect<T_>&       y) const = 0;

/// \brief Multiply matrix by vector <tt>a*x</tt> and add to <tt>y</tt>
    virtual void MultAdd(T_              a,
                         const Vect<T_>& x,
                         Vect<T_>&       y) const = 0;

/// \brief Multiply matrix by vector <tt>x</tt> and save in <tt>y</tt>
    virtual void Mult(const Vect<T_>& x,
                      Vect<T_>&       y) const = 0;

/// \brief Multiply transpose of matrix by vector <tt>x</tt> and save in <tt>y</tt>
    virtual void TMult(const Vect<T_>& v,
                       Vect<T_>&       w) const = 0;

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] x %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    virtual void Axpy(T_                a,
                      const Matrix<T_>* x) = 0;

/// \brief Initialize matrix storage in the case where only diagonal terms are stored.
/// \details This member function is to be used for explicit time integration schemes
    void setDiagonal(Mesh& mesh);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    virtual void setGraph(const Vect<RC>& I,
                          int             opt=1) = 0;

    virtual void setMesh(Mesh&  mesh,
                         size_t dof=0) = 0;

    void init_set_mesh(Mesh&  mesh,
                       size_t dof=0);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    virtual void setMesh(size_t dof,
                         Mesh&  mesh,
                         int    code=0) = 0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    virtual void setMesh(size_t dof,
                         size_t nb_eq,
                         Mesh&  mesh) = 0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// brief Set all matrix entries to zero
    virtual void clear();

/** \brief Assembly of element matrix into global matrix.
 *  \details Case where element matrix is given by a C-array.
 *  @param [in] el Pointer to element instance
 *  @param [in] a Element matrix as a C-array
 */
    void Assembly(const Element& el,
                  T_*            a);

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/** \brief Assembly of element matrix into global matrix.
 *  \details Case where element matrix is given by a DMatrix instance.
 *  @param [in] el Pointer to element instance
 *  @param [in] A Element matrix as a DMatrix instance
 */
//    void Assembly(const Element&     el,
//                  const DMatrix<T_>& A);


/** \brief Assembly of side matrix into global matrix.
 *  \details Case where side matrix is given by a DMatrix instance.
 *  @param [in] sd Pointer to side instance
 *  @param [in] a Side matrix as a DMatrix instance
 */
//    void Assembly(const Side&        sd,
//                  const DMatrix<T_>& a);

/** \brief Assembly of element matrix into global matrix for a Discontinuous Galerkin approximation
 *  \details Case where element matrix is given by a C-array.
 *  @param [in] el Pointer to element instance
 *  @param [in] a Element matrix as a C-array
 */
    void DGAssembly(const Element& el,
                    T_*            a);

    void DGAssembly(const Element&                                               el,
                    const LocalMatrix<T_,MAX_NB_ELEMENT_DOF,MAX_NB_ELEMENT_DOF>& a);

    void DGAssembly(const Side&                                            sd,
                    const LocalMatrix<T_,MAX_NB_SIDE_DOF,MAX_NB_SIDE_DOF>& a);

/** \brief Assembly of element matrix into global matrix for a Discontinuous Galerkin approximation
 *  \details Case where element matrix is given by a DMatrix instance.
 *  @param [in] el Pointer to element instance
 *  @param [in] a Element matrix as a DMatrix instance
 */
//    void DGAssembly(const Element&     el,
//                    const DMatrix<T_>& a);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Assembly of side matrix into global matrix.
 *  \details Case where side matrix is given by a C-array.
 *  @param [in] sd Pointer to side instance
 *  @param [in] a Side matrix as a C-array instance
 */
    void Assembly(const Side& sd,
                  T_*         a);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Prescribe(Mesh&           mesh,
                   Vect<T_>&       b,
                   const Vect<T_>& u,
                   int             flag=0);
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
    void Prescribe(Vect<T_>&       b,
                   const Vect<T_>& u,
                   int             flag=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Prescribe(int             dof,
                   int             code,
                   Mesh&           mesh,
                   Vect<T_>&       b,
                   const Vect<T_>& u,
                   int             flag=0);
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
    void Prescribe(int             dof,
                   int             code,
                   Vect<T_>&       b,
                   const Vect<T_>& u,
                   int             flag=0);

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
                   int       flag=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Prescribe(size_t          dof,
                   Mesh&           mesh,
                   Vect<T_>&       b,
                   const Vect<T_>& u,
                   int             flag=0);
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
    void Prescribe(size_t          dof,
                   Vect<T_>&       b,
                   const Vect<T_>& u,
                   int             flag=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Prescribe1(Mesh&           mesh,
                    Vect<T_>&       b,
                    const Vect<T_>& u,
                    int             flag=0);

    void Prescribe1(Vect<T_>&       b,
                    const Vect<T_>& u,
                    int             flag=0);

    void Prescribe(Mesh& mesh);

    void Prescribe();

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
    void PrescribeSide();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Constraint(const Mesh& mesh);
    void Constraint();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Add <tt>val</tt> to entry <tt>(i,j)</tt>.
    virtual void add(size_t    i,
                     size_t    j,
                     const T_& val) = 0;

/// \brief Factorize matrix. Available only if the storage class enables it.
    virtual int Factor() = 0;

/** \brief Solve the linear system.
 *  \details If the inherited class is SpMatrix, the function uses an iterative method
 *  once this one has been chosen. Otherwise, the method solves the linear system
 *  by factorization.
 */
    virtual int solve(Vect<T_>& b,
                      bool      fact=true) = 0;

/** \brief Solve the linear system.
 *  \details If the inherited class is SpMatrix, the function uses an iterative method
 *  once this one has been chosen. Otherwise, the method solves the linear system
 *  by factorization.
 *  @param [in] b Vect instance that contains right-hand side
 *  @param [out] x Vect instance that contains solution
 *  @param [in] fact Set to \c true if factorization is to be performed, \c false if not.
 *  [Default: <tt>true</tt>]
 *  @return
 *    <ul>
 *       <li><tt>0</tt> if solution was normally performed
 *       <li><tt>n</tt> if the <tt>n</tt>-th pivot is null\n
 *                      Solution is performed only is factorization has previouly been invoked.
 *    </ul>
 */
    virtual int solve(const Vect<T_>& b,
                      Vect<T_>&       x,
                      bool            fact=true) = 0;

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
                       Vect<T_>&       x);

/// \brief Return number of stored terms in matrix.
    size_t getLength() const;

/// \brief Say if matrix is diagonal or not
    int isDiagonal() const;

/// \brief Say if matrix is factorized or not.
/// \details If the matrix was not factorized, the class does not allow
/// solving by a direct solver.
    int isFactorized() const;

/// \brief Return Column index for column <tt>i</tt> (See the description for class SpMatrix).
    virtual size_t getColInd(size_t i) const;

/// \brief Return Row pointer for row <tt>i</tt> (See the description for class SpMatrix).
    virtual size_t getRowPtr(size_t i) const;

/** \brief Assign a value to an entry of the matrix
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] val Value to assign
 */
    virtual void set(size_t    i,
                     size_t    j,
                     const T_& val) = 0;

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
    std::vector<T_> _a, _aU, _diag;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   size_t              _dof_type, _nb_rows, _nb_cols, _size, _dof, _length;
   T_                  _zero, _temp;
   std::vector<size_t> _row_ptr, _col_ind, _ch;
   real_t              _penal;
   int                 _set_nodes, _set_elements, _set_sides, _is_diagonal;
   Mesh                *_theMesh;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
