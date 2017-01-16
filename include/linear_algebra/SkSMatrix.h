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

          Definition of class 'SkSMatrix' for Symmetric Skyline Matrix

  ==============================================================================*/


#ifndef __SKSMATRIX_H
#define __SKSMATRIX_H

#include "linear_algebra/Matrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file SkSMatrix.h
 *  \brief Definition file for class SkSMatrix.
 */

/*! \class SkSMatrix
 *  \ingroup VectMat
 *  \brief To handle symmetric matrices in skyline storage format.
 *
 * \details This template class allows storing and manipulating a symmetric matrix in skyline storage format.
 *
 * The matrix entries are stored column by column as in the following example:
 * 
 * @verbatim
       /                       \ 
       | a0   a1    0   0   a7 |
       |      a2   a3   0   a8 |
       | ...       a4  a5   a9 |
       |               a6  a10 |
       |                   a11 |
       \                       /
   @endverbatim
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 */

template<class T_> class SkSMatrix;
class Mesh;

template<class T_>
class SkSMatrix : public Matrix<T_>
{
   using Matrix<T_>::_nb_rows;
   using Matrix<T_>::_nb_cols;
   using Matrix<T_>::_size;
   using Matrix<T_>::_length;
   using Matrix<T_>::_zero;
   using Matrix<T_>::_dof_type;
   using Matrix<T_>::_temp;
   using Matrix<T_>::_fact;
   using Matrix<T_>::_a;
   using Matrix<T_>::_diag;
   using Matrix<T_>::_ch;
   using Matrix<T_>::_is_diagonal;
   using Matrix<T_>::_theMesh;

public:

    using Matrix<T_>::operator();

/// \brief Default constructor.
/// \details Initializes a zero-dimension matrix
    SkSMatrix();

/** \brief Constructor that initializes a dense symmetric matrix.
 *  \details Normally, for a dense matrix this is not the right class.
 *  @param [in] size Number of matrix rows (and columns).
 *  @param [in] is_diagonal Boolean to select if the matrix is diagonal or not [Default: false]
 */
    SkSMatrix(size_t size,
              int    is_diagonal=false);

/** \brief Constructor using mesh to initialize skyline structure of matrix.
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix
 *  structure is determined using all DOFs.
 *  @param [in] is_diagonal Boolean argument to say is the matrix is actually a diagonal matrix or not.
 */
    SkSMatrix(Mesh&  mesh,
              size_t dof=0,
              int    is_diagonal=false);

/** \brief Constructor that initializes skyline structure of matrix using vector of column height.
 *  @param [in] ColHt Vect instance that contains rows lengths of matrix.
 */
    SkSMatrix(const Vect<size_t>& ColHt);

/** \brief Constructor for a square matrix using non zero row and column indices
 *  @param [in] I Vector containing row indices
 *  @param [in] J Vector containing column indices
 *  @param [in] opt Flag indicating if vectors I and J are cleaned and ordered
 *  (opt=1) or not (opt=0).\n In the latter case, these vectors can contain
 *  the same contents more than once and are not necessarily ordered.
 */
    SkSMatrix(const Vect<size_t>& I,
              const Vect<size_t>& J,
              int                 opt=1);

/** \brief Constructor for a square matrix using non zero row and column indices
 *  @param [in] I Vector containing row indices
 *  @param [in] J Vector containing column indices
 *  @param [in] a Vector containing matrix entries in the same order than
 *  the one given by <tt>I</tt> and <tt>J</tt>
 *  @param [in] opt Flag indicating if vectors <tt>I</tt> and <tt>J</tt> are cleaned and ordered
 *  (<tt>opt=1</tt>) or not (<tt>opt=0</tt>).\n In the latter case, these vectors can contain
 *  the same contents more than once and are not necessarily ordered              
 */
    SkSMatrix(const Vect<size_t>& I,
              const Vect<size_t>& J,
              const Vect<T_>&     a,
              int                 opt=1);

/// \brief Copy Constructor
    SkSMatrix(const SkSMatrix<T_>& m);

/// \brief Destructor
    ~SkSMatrix() { }

/** \brief Determine mesh graph and initialize matrix.
 *  \details This member function is called by constructor with the same arguments
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix
 *  structure is determined using all DOFs.
 */
    void setMesh(Mesh&  mesh,
                 size_t dof=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Determine mesh graph and initialize matrix.
 *  \details This member function is called by constructor with the same arguments
 *  @param [in] mesh Mesh instance for which matrix graph is determined.
 *  @param [in] dof_type Type of support of dof. To be given among the enumerated
 *  values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>, <tt>SIDE_DOF</tt> 
 *  @param [in] dof1 Label of first degree of freedom for which numbering is performed.
 *  @param [in] dof2 Label of second degree of freedom for which numbering is performed.
 */
    void setGraph(Mesh&  mesh,
                  int    dof_type,
                  size_t dof1,
                  size_t dof2);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(size_t dof, 
                 Mesh&  mesh,
                 int    code=0);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(size_t dof,
                 size_t nb_eq,
                 Mesh&  mesh);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Determine matrix structure.
/// \details This member function calculates matrix structure using Mesh instance <tt>mesh</tt>.
    void setSkyline(Mesh& mesh);

/// \brief Store diagonal entries in a separate internal vector
    void setDiag();

/** \brief Assign a value to an entry ofthe matrix.
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] val Value to assign to <tt>a(i,j)</tt>
 */
    void set(size_t    i,
             size_t    j,
             const T_& val);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void SSet(size_t    i,
              size_t    j,
              const T_& val);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                   a,
              const SkSMatrix<T_>& m);

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m Pointer to %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                a,
              const Matrix<T_>* m);

/** \brief Multiply matrix by vector <tt>x</tt> and add to <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [in,out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(const Vect<T_>& x,
                 Vect<T_>&       y) const;

/** \brief Multiply matrix by vector <tt>a*x</tt> and add to <tt>y</tt>.
 *  @param [in] a Constant to multiply by matrix
 *  @param [in] x Vector to multiply by matrix
 *  @param [in,out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(T_              a,
                 const Vect<T_>& x,
                 Vect<T_>&       y) const;

/** \brief Multiply matrix by vector <tt>x</tt> and save in <tt>y</tt>
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void Mult(const Vect<T_>& x,
              Vect<T_>&       y) const;

/** \brief Multiply transpose of matrix by vector x and save in y
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void TMult(const Vect<T_>& x,
               Vect<T_>&       y) const;

/** \brief Add a constant to an entry of the matrix
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] val Constant value to add to <tt>a(i,j)</tt>
 */
    void add(size_t    i,
             size_t    j,
             const T_& val);

/// \brief Return column height.
/// \details Column height at entry <tt>i</tt> is returned.
    size_t getColHeight(size_t i) const;

/// \brief Get <tt>j</tt>-th column vector.
    Vect<T_> getColumn(size_t j) const;

/// \brief Get <tt>i</tt>-th row vector.
    Vect<T_> getRow(size_t i) const;

/** \brief Operator () (Non constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @warning To modify a value of an entry of the matrix it is safer not to
 *  modify both lower and upper triangles. Otherwise, wrong values will be 
 *  assigned. If not sure, use the member functions set or add.
 */
    T_ & operator()(size_t i,
                    size_t j);

/** \brief Operator () (Constant version).
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ operator()(size_t i,
                  size_t j) const;

/// \brief Operator =.
/// \details Copy matrix <tt>m</tt> to current matrix instance.
    SkSMatrix<T_> & operator=(const SkSMatrix<T_>& m);

/** \brief Operator =.
 *  \details define the matrix as a diagonal one with all diagonal entries equal
 *  to <tt>x</tt>.
 */
    SkSMatrix<T_> & operator=(const T_& x);

/// \brief Operator +=.
/// \details Add matrix <tt>m</tt> to current matrix instance.
    SkSMatrix<T_> & operator+=(const SkSMatrix<T_>& m);

/// \brief Operator *=.
/// \details Premultiply matrix entries by constant value <tt>x</tt>.
    SkSMatrix<T_> & operator*=(const T_& x);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int Factor() { return setLDLt(); }
    int Solve(Vect<T_>& b) { return solve(b); }
    int Solve(const Vect<T_>& b, Vect<T_>& x) { return solve(b,x); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Factorize matrix (LDLt (Crout) factorization).
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if factorization was normally performed
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null
 *  </ul>
 */
    int setLDLt();

/** \brief Solve a linear system using the LDLt (Crout) factorization
 *  \details This function solves a linear system. The LDLt factorization is 
 *  performed if this was not already done using the function setLU.
 *  @param [in] b Vect instance that contains right-hand side
 *  @param [out] x Vect instance that contains solution
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null
 *  </ul>
 *  Solution is performed only is factorization has previouly been invoked.
 */
    int solveLDLt(const Vect<T_>& b,
                  Vect<T_>&       x);

/** \brief Solve linear system.
 *  \details The linear system having the current instance as a matrix is solved by using the LDLt decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLU.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLDLt realizes the factorization step only.
 *  @param [in,out] b Vect instance that contains right-hand side on input and solution on output.
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(Vect<T_>& b);

/** \brief Solve linear system.
 *  \details The linear system having the current instance as a matrix is solved by using the LDLt decomposition.
 *  Solution is thus realized after a factorization step and a forward/backward substitution step.
 *  The factorization step is realized only if this was not already done.\n
 *  Note that this function modifies the matrix contents is a factorization is performed. Naturally, if the
 *  the matrix has been modified after using this function, the user has to refactorize it using the function setLDLt.
 *  This is because the class has no non-expensive way to detect if the matrix has been modified.
 *  The function setLDLt realizes the factorization step only.
 *  @param [in] b Vect instance that contains right-hand side.
 *  @param [out] x Vect instance that contains solution
 *  @return
 *  <ul>
 *     <li><tt>0</tt> if solution was normally performed,
 *     <li><tt>n</tt> if the <tt>n</tt>-th pivot is null.
 *  </ul>
 */
    int solve(const Vect<T_>& b,
              Vect<T_>&       x);

/// \brief Return C-Array
/// \details Skyline of matrix is stored row by row.
    T_ *get() const { return _a; }

/// \brief Assign a value to the i-th entry of C-array containing matrix
    void set(size_t i,
             T_     x)
    { _a[i] = x; }

/// \brief Return entry <tt>(i,j)</tt> of matrix if this one is stored, <tt>0</tt> else
    T_ get(size_t i,
           size_t j) const;

 private:

   int _dof;
};

///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////


template<class T_>
SkSMatrix<T_>::SkSMatrix() : _dof(0)
{
   _fact = false;
   _is_diagonal = false;
   _dof_type = NODE_DOF;
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(size_t size,
                         int    is_diagonal) : _dof(0)
{
   _zero = 0;
   _fact = false;
   _is_diagonal = is_diagonal;
   _dof_type = NODE_DOF;
   _nb_rows = _nb_cols = _size = size;
   _ch.resize(size);
   _diag.resize(_size);
   _ch[0] = 0;
   if (_is_diagonal) {
      for (size_t i=1; i<_size; i++)
         _ch[i] = _ch[i-1] + 1;
   }
   else {
      for (size_t i=1; i<_size; i++)
         _ch[i] = _ch[i-1] + i + 1;
   }
   _length = _ch[_size-1] + 1;
   _a.resize(_length);
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(Mesh&  mesh,
                         size_t dof,
                         int    is_diagonal)
{
   _fact = false;
   _is_diagonal = is_diagonal;
   setMesh(mesh,dof);
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(const Vect<size_t>& ColHt) : _dof(0)
{
   _fact = false;
   _is_diagonal = false;
   _dof_type = NODE_DOF;
   _size = ColHt.size();
   _zero = 0;
   _ch.resize(_size);
   _ch[0] = 0;
   for (size_t i=1; i<_size; ++i)
      _ch[i] = _ch[i-1] + ColHt[i];
   _length = _ch[_size-1] + 1;
   _a.resize(_length);
   _diag.resize(_size);
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(const Vect<size_t>& I,
                         const Vect<size_t>& J,
                         int                 opt) : _dof(0)
{
   _is_diagonal = false;
   _size = 0;
   _fact = false;
   _is_diagonal = false;
   _dof_type = NODE_DOF;
   size_t i;
   size_t n = I.size();
   std::vector<RC> pp(n);
   for (i=0; i<n; i++) {
      _size = std::max(_size,I[i]);
      pp[i] = RC(I[i]-1,J[i]-1);
   }
   _ch.resize(_size);
   _nb_rows = _nb_cols = _size;
   if (opt==0) {
      sort(pp.begin(),pp.end());
      vector<RC>::iterator new_end = std::unique(pp.begin(),pp.end());
      pp.erase(new_end,pp.end());
   }
   for (i=0; i<n; i++) {
      if (I[i]>J[i])
         _ch[I[i]-1] = std::max(static_cast<unsigned>(abs(int(I[i])-int(J[i]))),_ch[I[i]-1]);
   }
   _ch[0] = 0;
   for (i=1; i<_size; ++i)
      _ch[i] += _ch[i-1] + 1;
   _length = _ch[_size-1]+1;
   _a.resize(_length);
   _diag.resize(_size);
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(const Vect<size_t>& I,
                         const Vect<size_t>& J,
                         const Vect<T_>&     a,
                         int                 opt) : _dof(0)
{
   _is_diagonal = false;
   _fact = false;
   _dof_type = NODE_DOF;
   size_t i;
   size_t n = I.size();
   std::vector<RC> pp(n);
   _size = 0;
   for (i=0; i<n; i++) {
      _size = std::max(_size,I[i]);
      pp[i] = RC(I[i]-1,J[i]-1);
   }
   _ch.resize(_size,0);
   _nb_rows = _nb_cols = _size;
   if (opt==0) {
      sort(pp.begin(),pp.end());
      vector<RC>::iterator new_end = std::unique(pp.begin(),pp.end());
      pp.erase(new_end,pp.end());
   }
   for (i=0; i<n; i++)
      if (I[i]>J[i])
         _ch[I[i]-1] = std::max(static_cast<unsigned>(abs(int(I[i])-int(J[i]))),_ch[I[i]-1]);
   _ch[0] = 0;
   for (i=1; i<_size; ++i)
      _ch[i] += _ch[i-1] + 1;
   _length = _ch[_size-1]+1;
   _a.resize(_length);
   size_t k=0;
   for (i=0; i<n; i++)
      set(I[i],J[i],a[k++]);
   _diag.resize(_size);
}


template<class T_>
SkSMatrix<T_>::SkSMatrix(const SkSMatrix<T_>& m) : _dof(m._dof)
{
   _length = m._length;
   _size = m._size;
   _ch.resize(_size);
   _ch = m._ch;
   _a.resize(_length);
   _a = m._a;
   _diag.resize(_size);
   _diag = m._diag;
   _fact = m._fact;
   _zero = T_(0);
   _theMesh = m._theMesh;
   _is_diagonal = m._is_diagonal;
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
void SkSMatrix<T_>::setGraph(Mesh&  mesh,
                             int    dof_type,
                             size_t dof1,
                             size_t dof2)
{
   _theMesh = &mesh;
   _theMesh->selectDOF(dof_type,dof1,dof2);
   if (dof_type==NODE_DOF)
      _length = NodeSkyline(mesh,_ch,dof1,dof2);
   else if (dof_type==SIDE_DOF) {
      mesh.getAllSides();
      _length = SideSkyline(mesh,_ch,dof1,dof2);
   }
   else if (dof_type==ELEMENT_DOF)
      _length = ElementSkyline(mesh,_ch);
   else;
   _size = _ch.size();
   _diag.resize(_size);
   _a.resize(_length);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] += _ch[i-1];
   _fact = false;
   _dof = 0;
   _nb_rows = _nb_cols = _size;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


template<class T_>
void SkSMatrix<T_>::setMesh(Mesh&  mesh,
                            size_t dof)
{
   _dof_type = mesh.getDOFSupport();
   Matrix<T_>::init_set_mesh(mesh,dof);
   if (_dof_type==NODE_DOF)
      if (dof)
         _length = NodeSkyline(*_theMesh,_ch,dof);
      else
         _length = NodeSkyline(*_theMesh,_ch);
   else if (_dof_type==SIDE_DOF)
      if (dof)
         _length = SideSkyline(*_theMesh,_ch,dof);
      else
         _length = SideSkyline(*_theMesh,_ch);
   else if (_dof_type==ELEMENT_DOF)
      if (dof)
         _length = ElementSkyline(*_theMesh,_ch,dof);
      else
         _length = ElementSkyline(*_theMesh,_ch);
   else;
   _diag.resize(_size);
   _a.resize(_length);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] += _ch[i-1];
   _fact = false;
   _dof = 0;
   _nb_rows = _nb_cols = _size;
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
void SkSMatrix<T_>::setMesh(size_t dof,
                            Mesh&  mesh,
                            int    code)
{
// This is just to avoid warning on unused variable
   dof = 0;
   if (mesh.getDim()==0) { }
   code = 0;
   _fact = false;
   _theMesh = &mesh;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
void SkSMatrix<T_>::setMesh(size_t dof,
                            size_t nb_eq,
                            Mesh&  mesh)
{
// This is just to avoid warning on unused variable
   dof = 0;
   nb_eq = 0;
   if (mesh.getDim()==0) { }
   _fact = false;
   _theMesh = &mesh;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


template<class T_>
void SkSMatrix<T_>::setSkyline(Mesh& mesh)
{
   _zero = 0;
   int set_sides = mesh.SidesAreDOF();
   _size = mesh.getNbEq();
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
         _length = NodeSkyline(mesh,_ch,_dof);
      }
      _length = NodeSkyline(mesh,_ch);
   }
   _diag.resize(_size);
   _ch[0] = 0;
   for (size_t i=1; i<_size; i++)
      _ch[i] += _ch[i-1];
   _a.resize(_length);
   _fact = false;
   _theMesh = &mesh;
}


template<class T_>
void SkSMatrix<T_>::setDiag()
{
   for (size_t i=0; i<_size; i++)
      _diag[i] = _a[_ch[i]];
}


template<class T_>
void SkSMatrix<T_>::set(size_t    i,
                        size_t    j,
                        const T_& val)
{
   try {
      if (i>=j)
         _a[_ch[i-1]+j-i] = val;
      else
         THROW_RT("set(i,j,x): Index pair (" + itos(int(i)) + "," + itos(int(j)) + ") is not " +
                  "compatible with skyline symmeric storage.");
   }
   CATCH("SkSMatrix");
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
void SkSMatrix<T_>::SSet(size_t    i,
                         size_t    j,
                         const T_& val)
{
   _a[_ch[i-1]+j-i] = val;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


template<class T_>
void SkSMatrix<T_>::MultAdd(const Vect<T_>& x,
                            Vect<T_>&       y) const
{
   for (size_t i=0; i<_size; i++) {
      for (size_t j=0; j<_size; j++)
         y[i] += operator()(i+1,j+1)*x[j];
   }
}


template<class T_>
void SkSMatrix<T_>::MultAdd(T_              a,
                            const Vect<T_>& x,
                            Vect<T_>&       y) const
{
  for (size_t i=0; i<_size; i++)
     for (size_t j=0; j<_size; j++)
        y[i] += a * operator()(i+1,j+1)*x[j];
}


template<class T_>
void SkSMatrix<T_>::Mult(const Vect<T_>& x,
                         Vect<T_>&       y) const
{
   y = static_cast<T_>(0);
   MultAdd(x,y);
}


template<class T_>
void SkSMatrix<T_>::TMult(const Vect<T_>& x,
                          Vect<T_>&       y) const
{
   Mult(x,y);
}


template<class T_>
void SkSMatrix<T_>::add(size_t    i,
                        size_t    j,
                        const T_& val)
{
   if (i>=j)
      _a[_ch[i-1]+j-i] += val;
}


template<class T_>
size_t SkSMatrix<T_>::getColHeight(size_t i) const
{
   if (i==1)
      return 1;
   else
      return _ch[i-1]-_ch[i-2];
}


template<class T_>
Vect<T_> SkSMatrix<T_>::getColumn(size_t j) const
{
   Vect<T_> v(_nb_rows);
   for (size_t i=1; i<=_nb_rows; i++)
      v(i) = (*this)(i,j);
   return v;
}


template<class T_>
Vect<T_> SkSMatrix<T_>::getRow(size_t i) const
{
   return getColumn(i);
}


template<class T_>
T_ & SkSMatrix<T_>::operator()(size_t i,
                               size_t j)
{
   int k=0, l=0;
   if (i>1)
      k = int(j-i+_ch[i-1]-_ch[i-2]-1);
   if (j>1)
      l = int(i-j+_ch[j-1]-_ch[j-2]-1);
   if (k>=0 && i>=j)
      return _a[_ch[i-1]+j-i];
   else if (l>=0 && i<j)
      return _a[_ch[j-1]+i-j];
   else
      return _zero;
}


template<class T_>
T_ SkSMatrix<T_>::operator()(size_t i,
                             size_t j) const
{
   int k=0, l=0;
   if (i>1)
      k = int(j-i+_ch[i-1]-_ch[i-2]-1);
   if (j>1)
      l = int(i-j+_ch[j-1]-_ch[j-2]-1);
   if (k>=0 && i>=j)
      return _a[_ch[i-1]+j-i];
   else if (l>=0 && i<j)
      return _a[_ch[j-1]+i-j];
   else
      return _zero;
}


template<class T_>
SkSMatrix<T_> & SkSMatrix<T_>::operator=(const SkSMatrix<T_>& m)
{
   _a = m._a;
   return *this;
}


template<class T_>
SkSMatrix<T_> & SkSMatrix<T_>::operator=(const T_& x)
{
   _fact = false;
   for (size_t i=0; i<_length; i++)
      _a[i] = 0;
   for (size_t i=0; i<_nb_rows; i++) {
      _diag[i] = x;
      set(i+1,i+1,x);
   }
   return *this;
}


template<class T_>
SkSMatrix<T_> & SkSMatrix<T_>::operator+=(const SkSMatrix<T_>& m)
{
   _fact = false;
   for (size_t i=0; i<_length; i++)
      _a[i] += m._a[i];
   return *this;
}


template<class T_>
SkSMatrix<T_> & SkSMatrix<T_>::operator*=(const T_& x)
{
   _fact = false;
   for (size_t i=0; i<_length; i++)
      _a[i] *= x;
   return *this;
}


template<class T_>
int SkSMatrix<T_>::setLDLt()
{
   if (_is_diagonal)
      return 0;
   size_t di=0, dij;
   T_ s;
   T_ pivot = _a[_ch[0]];
   try {
      if (Abs(pivot) < OFELI_EPSMCH)
         THROW_RT("Factor(): The first pivot is null.");
      else
         _a[_ch[0]] = T_(1.)/pivot;
   }
   CATCH_EXIT("SkSMatrix");

   for (size_t i=1; i<_size; i++) {
      size_t dj = 0;
      for (size_t j=di=i+1+_ch[i-1]-_ch[i]; j<i; j++) {
         if (j>0)
            dj = j+1+_ch[j-1]-_ch[j];
         dij = std::max(di,dj);
         for (size_t k=0; k<j-dij; k++)
            _a[_ch[i]+j-i] -= _a[_ch[i]+dij+k-i]*_a[_ch[j]+dij+k-j];
      }

      pivot = _a[_ch[i]];
      for (size_t k=di; k<i; k++) {
         s = _a[_ch[i]+k-i]*_a[_ch[k]];
         pivot -= s*_a[_ch[i]+k-i];
         _a[_ch[i]+k-i] = s;
      }
      try {
         if (Abs(pivot) < OFELI_EPSMCH)
            THROW_RT("Factor(): The " + itos(int(i)+1) + "-th pivot is null.");
         else
            _a[_ch[i]] = T_(1.)/pivot;
      }
      CATCH_EXIT("SkSMatrix");
   }
   _fact = true;
   return 0;
}


template<class T_>
int SkSMatrix<T_>::solve(Vect<T_>& b)
{
   int ret = 0;
   if (_is_diagonal) {
      for (size_t i=0; i<_size; i++) {
         try {
            if (Abs(_a[i]) < OFELI_EPSMCH)
               THROW_RT("solve(b): The " + itos(i+1) + "-th diagonal is null.");
         }
         CATCH_EXIT("SkSMatrix");
         b[i] /= _a[i];
      }
      return 0;
   }
   if (!_fact)
      ret = setLDLt();
   for (size_t i=1; i<_size; i++) {
      size_t di = i+1+_ch[i-1]-_ch[i];
      T_ s = 0;
      for (size_t j=0; j<i-di; j++)
         s += _a[_ch[i]+di+j-i] * b[di+j];
      b[i] -= s;
   }
   for (size_t i=0; i<_size; i++)
      b[i] *= _a[_ch[i]];
   for (int k=int(_size-1); k>0; k--) {
      size_t di = k+1+_ch[k-1]-_ch[k];
      for (size_t j=0; j<k-di; j++)
         b[j+di] -= b[k] * _a[_ch[k]+di+j-k];
   }
   return ret;
}


template<class T_>
int SkSMatrix<T_>::solve(const Vect<T_>& b,
                         Vect<T_>&       x)
{
   x = b;
   return Solve(x);
}


template<class T_>
int SkSMatrix<T_>::solveLDLt(const Vect<T_>& b,
                             Vect<T_>&       x)
{
   int ret = 0;
   if (!_fact)
      ret = setLDLt();
   x = b;
   if (_is_diagonal) {
      for (size_t i=0; i<_size; i++)
         x[i] /= _a[i];
      return 0;
   }
   for (size_t i=1; i<_size; i++) {
      size_t di = i+1+_ch[i-1]-_ch[i];
      T_ s = 0;
      for (size_t j=0; j<i-di; j++)
         s += _a[_ch[i]+di+j-i] * x[di+j];
      x[i] -= s;
   }
   for (size_t i=0; i<_size; i++)
      x[i] *= _a[_ch[i]];
   for (int k=int(_size-1); k>0; k--) {
      size_t di = k+1+_ch[k-1]-_ch[k];
      for (size_t j=0; j<k-di; j++)
         x[j+di] -= x[k] * _a[_ch[k]+di+j-k];
   }
   return ret;
}


template<class T_>
T_ SkSMatrix<T_>::get(size_t i,
                      size_t j) const
{
   if (i>=j)
      return _a[_ch[i-1]+j-i];
   else
      return _a[_ch[j-1]+i-j];
}


template<class T_>
void SkSMatrix<T_>::Axpy(T_                   a,
                         const SkSMatrix<T_>& m)
{
   _a += a * m._a;
}


template<class T_>
void SkSMatrix<T_>::Axpy(T_                a,
                         const Matrix<T_>* m)
{
   for (size_t i=0; i<_length; i++)
      _a[i] += a * m->_a[i];
}


///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////


/** \fn Vect<T_> operator*(const SkSMatrix<T_> &A, const Vect<T_> &b)
 *  \brief Operator * (Multiply vector by matrix and return resulting vector
 *  \ingroup VectMat
 *  @param [in] A SkSMatrix instance to multiply by vector
 *  @param [in] b Vect instance 
 *  \return Vect instance containing <tt>A*b</tt>
 */
template<class T_>
Vect<T_> operator*(const SkSMatrix<T_>& A,
                   const Vect<T_>&      b)
{
   Vect<T_> v(b.size());
   A.Mult(b,v);
}


/** \fn ostream & operator<<(ostream& s, const SkSMatrix<T_> &a)
 *  \ingroup VectMat
 *  \brief Output matrix in output stream
 */
template<class T_>
ostream& operator<<(ostream&             s,
                    const SkSMatrix<T_>& a)
{
   s.setf(ios::right|ios::scientific);
   s << endl;
   for (size_t i=1; i<=a.getNbRows(); i++) {
      s << "\nRow:  " << setw(6) << i << endl;
      for (size_t j=1; j<=a.getNbColumns(); j++)
          s << "  " << setprecision(8) << std::setfill(' ') << setw(18) << a(i,j);
      s << endl;
   }
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
