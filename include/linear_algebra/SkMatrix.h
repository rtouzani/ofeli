/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

            Definition of class 'SkMatrix' for general skyline matrix

  ==============================================================================*/


#ifndef __SKMATRIX_H
#define __SKMATRIX_H

#include "linear_algebra/Matrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file SkMatrix.h
 *  \brief Definition file for class SkMatrix.
 */

class Mesh;
  
/*! \class SkMatrix
 *  \ingroup VectMat
 *  \brief To handle square matrices in skyline storage format.
 *
 *  \details This template class allows storing and manipulating a matrix in
 *  skyline storage format.
 *
 * The matrix entries are stored in 2 vectors column by column as in the following 
 * example:
 * 
 * @verbatim
 /                    \        /                       \ 
 | l0            .    |        | u0   u1    0   0   u7 |
 | l1  l2        .    |        |      u2   u3   0   u8 |
 |  0  l3  l4    .    |        | ...       u4  u5   u9 |
 |  0   0  l5  l6     |        |               u6  u10 |
 | l7  l8  l9 l10 l11 |        |                   u11 |
 \                    /        \                       /
 @endverbatim
 * 
 * \tparam T_ Data type (double, float, complex<double>, ...)
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

template<class T_>
class SkMatrix : public Matrix<T_>
{

 public:

   using Matrix<T_>::_size;
   using Matrix<T_>::_length;
   using Matrix<T_>::_nb_rows;
   using Matrix<T_>::_nb_cols;
   using Matrix<T_>::_zero;
   using Matrix<T_>::_temp;
   using Matrix<T_>::_fact;
   using Matrix<T_>::_a;
   using Matrix<T_>::_aU;
   using Matrix<T_>::_diag;
   using Matrix<T_>::_ch;
   using Matrix<T_>::_dof_type;
   using Matrix<T_>::_is_diagonal;
   using Matrix<T_>::_theMesh;
   using Matrix<T_>::operator();

/// \brief Default constructor.
/// \details Initializes a zero-dimension matrix
    SkMatrix()
    {
       _dof = 0;
       _fact = false;
       _is_diagonal = false;
    }

/** \brief Constructor that initializes a dense symmetric matrix.
 *  \details Normally, for a dense matrix this is not the right class.
 *  @param [in] size Number of matrix rows (and columns).
 *  @param [in] is_diagonal Boolean to select if the matrix is diagonal or not [Default: false]
 */
    SkMatrix(size_t size,
             int    is_diagonal=false)
    {
       _dof = 0;
       _zero = 0;
       _fact = false;
       _is_diagonal = is_diagonal;
       _dof_type = NODE_DOF;
       _nb_rows = _nb_cols = _size = size;
       _ch.resize(size);
       _diag.resize(_size);
       _ch[0] = 0;
       for (size_t i=1; i<_size; i++)
          _ch[i] = _ch[i-1] + i + 1;
       if (_is_diagonal) {
          for (size_t i=1; i<_size; i++)
             _ch[i] = _ch[i-1] + 1;
       }
       _length = _ch[_size-1] + 1;
       _a.resize(_length);
       _aU.resize(_length);
    }

/** \brief Constructor using mesh to initialize skyline structure of matrix.
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix
 *  structure is determined using all DOFs.
 *  @param [in] is_diagonal Boolean argument to say is the matrix is actually a diagonal matrix or not.
 */
    SkMatrix(Mesh&  mesh,
             size_t dof=0,
             int    is_diagonal=false)
    {
       _is_diagonal = is_diagonal;
       _fact = false;
       setMesh(mesh,dof);
    }

/** \brief Constructor that initializes skyline structure of matrix using vector of column heights.
 *  @param [in] ColHt Vect instance that contains rows lengths of matrix.
 */
    SkMatrix(const Vect<size_t> &ColHt) : _dof(0)
    {
       _is_diagonal = false;
       _zero = 0;
       _size = ColHt.size();
       _ch.resize(_size,0);
       for (size_t i=1; i<_size; i++)
          _ch[i] = _ch[i-1] + ColHt[i];
       _length = _ch[_size-1] + 1;
       _a.resize(_length);
       _aU.resize(_length);
       _diag.resize(_size);
       _fact = false;
    }

/// \brief Copy Constructor
    SkMatrix(const SkMatrix<T_>& m)
    {
       _is_diagonal = m._is_diagonal;
       _size = m._size;
       _length = m._length;
       _ch.resize(_size);
       _ch = m._ch;
       _diag.resize(_size);
       _diag = m._diag;
       _a.resize(_length);
       _aU.resize(_length);
       _a = m._a;
       _aU = m._aU;
       _fact = m._fact;
       _dof = m._dof;
       _zero = static_cast<T_>(0);
    }

/// \brief Destructor
    ~SkMatrix() { }

/** \brief Determine mesh graph and initialize matrix.
 *  \details This member function is called by constructor with the same arguments
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix
 *  structure is determined using all DOFs.
 */
    void setMesh(Mesh&  mesh,
                 size_t dof=0)
    {
       _dof_type = mesh.getDOFSupport();
       Matrix<T_>::init_set_mesh(mesh,dof);
       if (_dof_type==NODE_DOF) {
          if (dof)
             _length = NodeSkyline(mesh,_ch,dof);
          else
             _length = NodeSkyline(mesh,_ch);
       }
       else if (_dof_type==SIDE_DOF) {
          if (dof)
             _length = SideSkyline(mesh,_ch,dof);
          else
             _length = SideSkyline(mesh,_ch);
       }
       else if (_dof_type==ELEMENT_DOF) {
          if (dof)
             _length = ElementSkyline(mesh,_ch,dof);
          else
             _length = ElementSkyline(mesh,_ch);
       }
       else
          ;
       _diag.resize(_size);
       _a.resize(_length);
       _aU.resize(_length);
       _ch[0] = 0;
       for (size_t i=1; i<_size; i++)
           _ch[i] += _ch[i-1];
        _zero = T_(0);
        _fact = false;
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(size_t dof,
                 Mesh&  mesh,
                 int    code=0)
    {
//     This is just to avoid warning on unused variable
       dof = 0;
       if (mesh.getDim()==0) { }
       code = 0;
       _fact = false;
       _theMesh = &mesh;
    }

    void setMesh(size_t dof,
                 size_t nb_eq,
                 Mesh&  mesh)
    {
//     This is just to avoid warning on unused variable
       dof = 0;
       nb_eq = 0;
       if (mesh.getDim()==0) { }
       _fact = false;
       _theMesh = &mesh;
    }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Determine matrix structure.
 *  \details This member function calculates matrix structure using a Mesh instance.
 *  @param [in] mesh Mesh instance
 */
    void setSkyline(Mesh& mesh)
    {
       _zero = 0;
       int set_sides = mesh.SidesAreDOF();
       _size = mesh.getNbEq();
       _theMesh = &mesh;
       if (_dof)
          _size = mesh.getNbNodes();
       if (set_sides) {
          if (_dof) {
             _size = mesh.getNbSides();
             _length = SideSkyline(mesh,_ch,_dof);
          }
          _length = SideSkyline(mesh,_ch);
       }
       else {
          if (_dof) {
             _size = mesh.getNbNodes();
             _length = SideSkyline(mesh,_ch,_dof);
          }
          _length = SideSkyline(mesh,_ch);
       }
       _diag.resize(_size);
       _ch[0] = 0;
       for (size_t i=1; i<_size; i++)
          _ch[i] += _ch[i-1];
       _a.resize(_length);
       _aU.resize(_length);
       _fact = false;
    }

/// \brief Store diagonal entries in a separate internal vector.
    void setDiag()
    {
       for (size_t i=0; i<_size; i++)
          _diag[i] = _aU[_ch[i]];
    }

/** \brief Choose DOF to activate.
 *  \details This function is available only if variable <tt>dof</tt> is equal to 1 in the constructor
 *  @param [in] i Index of the DOF
 */
    void setDOF(size_t i) { _dof = i; }

/** \brief Assign a value to an entry ofthe matrix.
 *  @param [in] i Row index (starting at <tt>i=1</tt>)
 *  @param [in] j Column index (starting at <tt>i=1</tt>)
 *  @param [in] val Value to assign to entry <tt>a(i,j)</tt>
 */
    void set(size_t    i,
             size_t    j,
             const T_& val)
    {
       int k=0, l=0;
       if (i>1)
          k = int(j-i+_ch[i-1]-_ch[i-2]-1);
       if (j>1)
          l = int(i-j+_ch[j-1]-_ch[j-2]-1);
       if (k>=0 && i>j)
          _a[_ch[i-1]+j-i] = val;
       else if (l>=0 && i<=j)
          _aU[_ch[j-1]+i-j] = val;
       else
          throw OFELIException("In SkMatrix::Set(i,j,x): Index pair: ("+itos(int(i))+"," +
                               itos(int(j))+") is not " + "compatible with skyline symmeric storage.");
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void SSet(size_t    i,
              size_t    j,
              const T_& val)
    {
       int k=0, l=0;
       if (i>1)
          k = j-i+_ch[i-1]-_ch[i-2]-1;
       if (j>1)
          l = i-j+_ch[j-1]-_ch[j-2]-1;
       if (k>=0 && i>j)
          _a[_ch[i-1]+j-i] = val;
       else if (l>=0 && i<=j)
          _aU[_ch[j-1]+i-j] = val;
    }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Add to matrix the product of a matrix by a scalar
  *  @param [in] a Scalar to premultiply
  *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
  *  to current instance
  */
    void Axpy(T_                  a,
              const SkMatrix<T_>& m)
    {
       _a  += a * m._a;
       _aU += a * m._aU;
    }

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                a,
              const Matrix<T_>* m)
    {
       for (size_t i=0; i<_length; i++) {
          _a[i]  += a * m->_a[i];
          _aU[i] += a * m->_aU[i];
       }
    }

/** \brief Multiply matrix by vector <tt>x</tt> and add to <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [in,out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(const Vect<T_>& x,
                 Vect<T_>&       y) const
    {
       for (size_t i=0; i<_size; i++)
          for (size_t j=0; j<_size; j++)
             y[i] += operator()(i+1,j+1)*x[j];
    }

/** \brief Multiply transpose of matrix by vector <tt>x</tt> and add to <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [in,out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void TMultAdd(const Vect<T_>& x,
                  Vect<T_>&       y) const
    {
       cerr << "TMultAdd is not implemented for class SkMatrix" << endl;
    }

/** \brief Multiply matrix by a vector and add to another one.
 *  @param [in] a Constant to multiply by matrix
 *  @param [in] x Vector to multiply by matrix
 *  @param [in,out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(T_              a,
                 const Vect<T_>& x,
                 Vect<T_>&       y) const
    {
       for (size_t i=0; i<_size; i++)
          for (size_t j=0; j<_size; j++)
             y[i] += a * operator()(i+1,j+1)*x[j];
    }

/** \brief Multiply matrix by vector <tt>x</tt> and save in <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void Mult(const Vect<T_>& x,
              Vect<T_>&       y) const
    {
       y = static_cast<T_>(0);
       MultAdd(x,y);
    }

/** \brief Multiply transpose of matrix by vector <tt>x</tt> and save in <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void TMult(const Vect<T_>& x,
               Vect<T_>&       y) const
    {
       y = static_cast<T_>(0);
       TMultAdd(x,y);
    }

/** \brief Add a constant value to an entry ofthe matrix.
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] val Constant value to add to <tt>a(i,j)</tt>
 */
    void add(size_t    i,
             size_t    j,
             const T_& val)
    {
       if (i>j)
          _a[_ch[i-1]+j-i] += val;
       else if (i<=j)
          _aU[_ch[j-1]+i-j] += val;
    }

/// \brief Return column height.
/// \details Column height at entry <tt>i</tt> is returned.
    size_t getColHeight(size_t i) const
    {
       if (i==1)
          return 1;
       else
          return _ch[i-1]-_ch[i-2];
    }

/** \brief Operator () (Constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ operator()(size_t i,
                  size_t j) const
    {
       int k=0, l=0;
       if (i>1)
          k = int(j-i+_ch[i-1]-_ch[i-2]-1);
       if (j>1)
          l = int(i-j+_ch[j-1]-_ch[j-2]-1);
       if (k>=0 && i>j)
          return _a[_ch[i-1]+j-i];
       else if (l>=0 && i<=j)
          return _aU[_ch[j-1]+i-j];
       else
          return _zero;
    }

/** \brief Operator () (Non constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ & operator()(size_t i,
                    size_t j)
    {
       int k=0, l=0;
       if (i>1)
          k = int(j-i+_ch[i-1]-_ch[i-2]-1);
       if (j>1)
          l = int(i-j+_ch[j-1]-_ch[j-2]-1);
       if (k>=0 && i>j)
          return _a[_ch[i-1]+j-i];
       else if (l>=0 && i<=j)
          return _aU[_ch[j-1]+i-j];
       else
          throw OFELIException("In SkMatrix::Operator(): Index pair (" + itos(int(i)) + "," +
                               itos(int(j)) + ") is not compatible with skyline structure");
       return _temp;
    }

/** \brief Impose an essential boundary condition.
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in] mesh Mesh instance from which information is extracted.
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [in] u Vect instance that conatins imposed valued at DOFs where they are to be imposed.
 *  @param [in] flag Parameter to determine whether only the right-hand side is to
 *  be modified (<tt>dof>0</tt>)\n
 *  or both matrix and right-hand side (<tt>dof=0</tt>, default value).
 */
    void DiagPrescribe(Mesh&           mesh,
                       Vect<T_>&       b,
                       const Vect<T_>& u,
                       int             flag=0)
    {
       real_t p=0;
       for (size_t l=0; l<_size; l++)
          p = std::max(p,_aU[_ch[l]]);
       MeshNodes(mesh) {
          for (size_t i=1; i<=TheNode.getNbDOF(); ++i) {
             if (TheNode.getCode(i)>0) {
                size_t ii = TheNode.getDOF(i)-1;
                for (size_t j=ii+1+_ch[ii-1]-_ch[ii]; j<=ii; j++) {
                   b[ii] = p*u[ii];
                   _a[_ch[ii]+j-ii] = _aU[_ch[ii]+j-ii] = 0;
                }
                _diag[ii] = _aU[_ch[ii]] = p;
             }
          }
       }
    }

/** \brief Impose an essential boundary condition using the Mesh instance provided by the constructor.
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [in] u Vect instance that conatins imposed valued at DOFs where they are to be imposed.
 *  @param [in] flag Parameter to determine whether only the right-hand side is to
 *  be modified (<tt>dof>0</tt>)\n
 *  or both matrix and right-hand side (<tt>dof=0</tt>, default value).
 */
    void DiagPrescribe(Vect<T_>&       b,
                       const Vect<T_>& u,
                       int             flag=0)
    {
       real_t p=0;
       for (size_t l=0; l<_size; l++)
          p = std::max(p,_aU[_ch[l]]);
       MESH_ND {
          for (size_t i=1; i<=TheNode.getNbDOF(); ++i) {
             if (TheNode.getCode(i)>0) {
                size_t ii = TheNode.getDOF(i)-1;
                for (size_t j=ii+1+_ch[ii-1]-_ch[ii]; j<=ii; j++) {
                   b[ii] = p*u[ii];
                   _a[_ch[ii]+j-ii] = _aU[_ch[ii]+j-ii] = 0;
                }
                _diag[ii] = _aU[_ch[ii]] = p;
             }
          }
       }
    }

/// \brief Operator =.
/// \details Copy matrix <tt>m</tt> to current matrix instance.
    SkMatrix<T_> & operator=(const SkMatrix<T_>& m)
    {
       _a = m._a;
       _aU = m._aU;
       return *this;
    }

/** \brief Operator =.
 *  \details define the matrix as a diagonal one with all diagonal entries equal
 *  to <tt>x</tt>.
 */
    SkMatrix<T_> & operator=(const T_& x)
    {
       _fact = false;
       for (size_t i=0; i<_length; i++)
          _a[i]  = _aU[i] = 0;
       for (size_t i=0; i<_nb_rows; i++) {
          _diag[i] = x;
          set(i+1,i+1,x);
       }
       return *this;
    }

/// \brief Operator +=.
/// \details Add matrix <tt>m</tt> to current matrix instance.
    SkMatrix<T_> & operator+=(const SkMatrix<T_>& m)
    {
       _fact = false;
       _a  += m._a;
       _aU += m._aU;
       return *this;
    }

/// \brief Operator +=.
/// \details Add constant value <tt>x</tt> to matrix entries.
    SkMatrix<T_> & operator+=(const T_& x)
    {
       _fact = false;
       _a += x;
       _aU += x;
       return *this;
    }

/// \brief Operator *=.
/// \details Premultiply matrix entries by constant value <tt>x</tt>.
    SkMatrix<T_> & operator*=(const T_& x)
    {
       _fact = false;
       _a  *= x;
       _aU *= x;
       return *this;
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int Factor() { return setLU(); }
    int Solve(Vect<T_>& b) { return solve(b); }
    int Solve(const Vect<T_>& b, Vect<T_>& x) { return solve(b,x); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Factorize the matrix (LU factorization)
 *  \details LU factorization of the matrix is realized. Note that since this
 *  is an in place factorization, the contents of the matrix are modified.
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if factorization was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 *  @remark A flag in this class indicates after factorization that this one
 *  has been realized, so that, if the member function solve is called after this
 *  no further factorization is done.
 */
    int setLU()
    {
       if (_is_diagonal)
          return 0;
       size_t k, di, dij, i=0;
       if (Abs(_aU[_ch[0]]) < OFELI_EPSMCH)
          throw OFELIException("In SkMatrix::Factor(): The first pivot is null.");
       for (i=1; i<_size; i++) {
          size_t dj = 0;
          for (size_t j=di=i+1+_ch[i-1]-_ch[i]; j<i; j++) {
             if (j>0)
                dj = j+1+_ch[j-1]-_ch[j];
             dij = std::max(di,dj);
             for (k=0; k<j-dij; k++)
                _a[_ch[i]+j-i] -= _a[_ch[i]+dij+k-i]*_aU[_ch[j]+dij+k-j];
             _a[_ch[i]+j-i] /= _aU[_ch[j]];
             for (k=0; k<j-dij; k++)
                _aU[_ch[i]+j-i] -= _a[_ch[j]+dij+k-j]*_aU[_ch[i]+dij+k-i];
          }
          for (k=0; k<i-di; k++)
             _aU[_ch[i]] -= _a[_ch[i]+k+di-i]*_aU[_ch[i]+k+di-i];
          if (Abs(_aU[_ch[i]]) < OFELI_EPSMCH)
             throw OFELIException("In SkMatrix::Factor(): The " + itos(i+1) + "-th pivot is null.");
       }
       _fact = true;
       return 0;
    }

/** \brief Solve linear system.
 *  \details The linear system having the current instance as a matrix is solved by using the LU decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLU.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLU realizes the factorization step only.
 *  @param [in,out] b Vect instance that contains right-hand side on input and solution on output.
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(Vect<T_>& b)
    {
       int ret = 0;
       if (_is_diagonal) {
          for (size_t i=0; i<_size; i++) {
             if (Abs(_aU[i]) < OFELI_EPSMCH)
                throw OFELIException("In SkMatrix::solve(b): The " + itos(i+1) + "-th diagonal is null.");
             b[i] /= _aU[i];
          }
          return 0;
       }
       size_t di;
       if (!_fact)
          ret = setLU();
       size_t i, j;
       for (i=1; i<_size; i++) {
          di = i+1+_ch[i-1]-_ch[i];
          T_ s = 0;
          for (j=0; j<i-di; j++)
             s += _a[_ch[i]+di+j-i] * b[di+j];
          b[i] -= s;
       }
       for (int k=int(_size-1); k>0; k--) {
          if (Abs(_aU[_ch[k]]) < OFELI_EPSMCH)
             throw OFELIException("In SkMatrix::solve(b): The " + itos(k+1) + "-th pivot is null.");
          b[k] /= _aU[_ch[k]];
          di = k+1+_ch[k-1]-_ch[k];
          for (j=0; j<k-di; j++)
             b[j+di] -= b[k] * _aU[_ch[k]+di+j-k];
       }
       b[0] /= _aU[_ch[0]];
       return ret;
    }

/** \brief Solve linear system.
 *  \details The linear system having the current instance as a matrix is solved by using the LU decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLU.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLU realizes the factorization step only.
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [out] x Vect instance that contains solution
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(const Vect<T_>& b,
              Vect<T_>&       x)
    {
       x = b;
       return solve(x);
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int solveLU(const Vect<T_>& b,
                Vect<T_>&       x)
    {
       int ret = 0;
       if (!_fact)
          ret = setLU();
       x = b;
       if (_is_diagonal) {
          for (size_t i=0; i<_size; i++)
             x[i] /= _aU[i];
          return 0;
       }
       size_t di, i, j;
       for (i=1; i<_size; i++) {
          di = i+1+_ch[i-1]-_ch[i];
          T_ s = 0;
          for (j=0; j<i-di; j++)
             s += _a[_ch[i]+di+j-i] * x[di+j];
          x[i] -= s;
       }
       for (int k=int(_size-1); k>0; k--) {
          if (Abs(_aU[_ch[k]]) < OFELI_EPSMCH)
             throw OFELIException("In SkMatrix::solveLU(b,x): The " + itos(k+1) + "-th pivot is null.");
          x[k] /= _aU[_ch[k]];
          di = k+1+_ch[k-1]-_ch[k];
          for (j=0; j<k-di; j++)
             x[j+di] -= b[k] * _aU[_ch[k]+di+j-k];
       }
       x[0] /= _aU[_ch[0]];
       return ret;
    }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return C-Array.
/// \details Skyline of matrix is stored row by row.
    T_ *get() const { return _a; }
   
/// \brief Return entry <tt>(i,j)</tt> of matrix if this one is stored, 0 else
    T_ get(size_t i,
           size_t j) const
    {
       if (i>j)
          return _a[_ch[i-1]+j-i];
       else
          return _aU[_ch[j-1]+i-j];
    }

 private:
    int _dof;
};

///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////


/** \fn Vect<T_> operator*(const SkMatrix<T_> &A, const Vect<T_> &b)
 *  \brief Operator * (Multiply vector by matrix and return resulting vector
 *  \ingroup VectMat
 *  @param [in] A SkMatrix instance to multiply by vector
 *  @param [in] b Vect instance 
 *  \return Vect instance containing <tt>A*b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
Vect<T_> operator*(const SkMatrix<T_>& A,
                   const Vect<T_>&     b)
{
   Vect<T_> v(b.size());
   A.Mult(b,v);
}


/** \fn ostream & operator<<(ostream &s, const SkMatrix<T_> &a)
 *  \ingroup VectMat
 *  \brief Output matrix in output stream
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
ostream& operator<<(ostream&            s,
                    const SkMatrix<T_>& a)
{
   s.setf(ios::right|ios::scientific);
   s << endl;
   for (size_t i=1; i<=a.getNbRows(); i++) {
      s << "\nRow  " << setw(6) << i << endl;
      for (size_t j=1; j<=a.getNbRows(); j++)
          s << "  " << setprecision(8) << std::setfill(' ') << setw(18) << a(i,j);
      s << endl;
   }
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
