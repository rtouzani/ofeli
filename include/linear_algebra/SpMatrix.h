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

                  Definition of class 'SpMatrix' for Sparse matrix

  ==============================================================================*/


#ifndef __SPMATRIX_H
#define __SPMATRIX_H

#include "linear_algebra/Matrix.h"
#ifdef USE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace Eigen;
#else
#include "solvers/LinearSolver.h"
#endif

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file SpMatrix.h
 *  \brief Definition file for class SpMatrix.
 */

template<class T_> class Vect;
template<class T_> class LinearSolver;
class Mesh;

using std::pair;
using std::unique;

/*! \class SpMatrix
 *  \ingroup VectMat
 *  \brief To handle matrices in sparse storage format.
 *
 * This template class enables storing and manipulating a sparse matrix,
 * i.e. only nonzero terms are stored. Internally, the matrix is stored as a
 * vector instance and uses for the definition of its graph a <tt>Vect<size_t></tt>
 * instance row_ptr and a Vect<size_t> instance <tt>col_ind</tt> that contains
 * respectively addresses of first element of each row and column indices.
 *
 * To illustrate this, consider the matrix
 * \verbatim
            1   2   0
            3   4   0
            0   5   0
   \endverbatim
 *
 * Such a matrix is stored in the vector<real_t> instance {1,2,3,4,5}.
 * The vectors <tt>row_ptr</tt> and <tt>col_ind</tt> are respectively:
 * <tt>{0,2,4,5}</tt>, <tt>{1,2,1,2,2}</tt>
 *
 * When the library <tt>eigen</tt> is used in conjunction with <tt>OFELI</tt>,
 * the class uses the sparse matrix class of <tt>eigen</tt> and enables then 
 * access to specific solvers (see class LinearSolver)
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

template<class T_> class SpMatrix : public Matrix<T_>
{
    using Matrix<T_>::_nb_rows;
    using Matrix<T_>::_nb_cols;
    using Matrix<T_>::_size;
    using Matrix<T_>::_length;
    using Matrix<T_>::_zero;
    using Matrix<T_>::_temp;
    using Matrix<T_>::_a;
    using Matrix<T_>::_diag;
    using Matrix<T_>::_ch;
    using Matrix<T_>::_dof_type;
    using Matrix<T_>::_is_diagonal;
    using Matrix<T_>::_theMesh;
    using Matrix<T_>::_row_ptr;
    using Matrix<T_>::_col_ind;

#ifdef USE_EIGEN
    typedef Eigen::Matrix<T_,Eigen::Dynamic,1> VectorX;
    typedef SparseMatrix<T_>                   SpMat;
    typedef Triplet<real_t>                    Tr;
#endif

 public:

/// \brief Default constructor.
/// \details Initialize a zero-dimension matrix
    SpMatrix() : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
    {
       _max_it = 1000;
       _toler = 1.e-8;
       _is_diagonal = 0;
    }

/** \brief Constructor that initializes current instance as a dense matrix.
 *  \details Normally, for a dense matrix this is not the right class.
 *  @param [in] nr Number of matrix rows.
 *  @param [in] nc Number of matrix columns.
 */
    SpMatrix(size_t nr,
             size_t nc) : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
    {
       _is_diagonal = 0;
       _nb_rows = nr;
       _nb_cols = nc;
       _size = 0;
       if (_nb_rows==_nb_cols)
          _size = _nb_rows;
       _length = _nb_rows*_nb_cols;
#ifdef USE_EIGEN
       _A.resize(_nb_rows,_nb_cols);
#else
       _row_ptr.resize(_nb_rows+1);
       _col_ind.resize(_length);
       _row_ptr[0] = 0;
       for (size_t i=0; i<_nb_rows; i++)
          _row_ptr[i+1] = _row_ptr[i] + _nb_cols;
       size_t l=0;
       for (size_t i=0; i<_nb_rows; i++)
          for (size_t j=0; j<_nb_cols; j++)
             _col_ind[l++] = j+1;
       _a.resize(_length,T_(0.));
#endif
       _max_it = 1000;
       _toler = 1.e-8;
    }

/** \brief Constructor that initializes current instance as a dense matrix.
 *  \details Normally, for a dense matrix this is not the right class.
 *  @param [in] size Number of matrix rows (and columns).
 *  @param [in] is_diagonal Boolean argument to say is the matrix is actually a diagonal matrix or not.
 */
    SpMatrix(size_t size,
             int    is_diagonal=false) : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
   {
       _is_diagonal = is_diagonal;
       _size = _nb_rows = _nb_cols = size;
#ifdef USE_EIGEN
       _A.resize(_nb_rows,_nb_cols);
#else
       _length = _size*_size;
       _row_ptr.resize(_nb_rows+1);
       _col_ind.resize(_length);
       _row_ptr[0] = 0;
       for (size_t i=0; i<_nb_rows; i++)
          _row_ptr[i+1] = _row_ptr[i] + _nb_cols;
       size_t l = 0;
       for (size_t i=0; i<_nb_rows; i++)
          for (size_t j=0; j<_nb_cols; j++)
             _col_ind[l++] = j+1;
       _a.resize(_length,T_(0.));
#endif
       _max_it = 1000;
       _toler = 1.e-8;
    }

/** \brief Constructor using a Mesh instance.
 *  @param [in] mesh Mesh instance from which matrix graph is extracted.
 *  @param [in] dof Option parameter, with default value <tt>0</tt>.\n
 *  <tt>dof=1</tt> means that only one degree of freedom for each node (or element or side)
 *  is taken to determine matrix structure. The value <tt>dof=0</tt> means that matrix 
 *  structure is determined using all DOFs.
 *  @param [in] is_diagonal Boolean argument to say is the matrix is actually a diagonal matrix or not.
 */
    SpMatrix(Mesh&  mesh,
             size_t dof=0,
             int    is_diagonal=false) : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
    {
       _is_diagonal = is_diagonal;
       setMesh(mesh,dof);
       _max_it = 1000;
       _toler = 1.e-8;
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    SpMatrix(size_t dof,
             Mesh&  mesh,
             int    code=0) : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
    {
       _is_diagonal = false;
       setMesh(dof,mesh,mesh.getDOFSupport());
       _max_it = 1000;
       _toler = 1.e-8;
    }

#ifndef USE_EIGEN
    SpMatrix(size_t dof,
             size_t nb_eq,
             Mesh&  mesh) : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
{
   _is_diagonal = false;
   setMesh(dof,nb_eq,mesh);
   _max_it = 1000;
   _toler = 1.e-8;
}
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef USE_EIGEN
/** \brief Constructor for a square matrix using non zero row and column indices
 *  @param [in] I Vector containing pairs of row and column indices
 *  @param [in] opt Flag indicating if vectors I is cleaned and ordered
 *  (opt=1) or not (opt=0). In the latter case, this vector can have
 *  the same contents more than once and are not necessarily ordered
 */
    SpMatrix(const Vect<RC>& I,
             int             opt=1) : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
    {
       _is_diagonal = false;
       setGraph(I,opt);
       _max_it = 1000;
       _toler = 1.e-8;
}
#endif

#ifndef USE_EIGEN
/** \brief Constructor for a square matrix using non zero row and column indices
 *  @param [in] I Vector containing pairs of row and column indices
 *  @param [in] a Vector containing matrix entries in the same order than
 *  the one given by <tt>I</tt>
 *  @param [in] opt Flag indicating if vector <tt>I</tt> is cleaned and ordered
 *  (<tt>opt=1</tt>: default) or not (<tt>opt=0</tt>). In the latter case, this vector can have
 *  the same contents more than once and are not necessarily ordered
 */
    SpMatrix(const Vect<RC>& I,
             const Vect<T_>& a,
             int             opt=1) : _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
    {
       _is_diagonal = false;
       size_t n=I.size();
       _nb_rows = _nb_cols = 0;
       for (size_t i=0; i<n; i++) {
          _nb_rows = std::max(_nb_rows,I[i].first);
          _nb_cols = std::max(_nb_cols,I[i].second);
          _IJ.push_back(RC(I[i].first-1,I[i].second-1));
       }
       if (_nb_rows==_nb_cols)
          _size = _nb_rows;
       if (opt==0) {
          sort(_IJ.begin(),_IJ.end());
          vector<RC>::iterator new_end=unique(_IJ.begin(),_IJ.end());
          _IJ.erase(new_end,_IJ.end());
       }
       _row_ptr.resize(_size+1);
       _col_ind.resize(n);
       StoreGraph(_size,_IJ,_row_ptr,_col_ind);
       _length = _IJ.size();
       _a.resize(_length);
       for (size_t j=0; j<n; j++)
          _a[_row_ptr[I[j].first-1]+_col_index(I[j].first,I[j].second)] = a[j];
       _max_it = 1000;
       _toler = 1.e-8;
    }
#endif

#ifndef USE_EIGEN
/** \brief Constructor for a rectangle matrix.
 *  @param [in] nr Number of rows
 *  @param [in] nc Number of columns
 *  @param [in] row_ptr Vector of row pointers (See the above description of this class).
 *  @param [in] col_ind Vector of column indices (See the above description of this class).\n
 */
    SpMatrix(size_t                nr,
             size_t                nc,
             const vector<size_t>& row_ptr,
             const vector<size_t>& col_ind)
             : _type(0), _dof(0), _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
    {
       _is_diagonal = false;
       _nb_rows = nr;
       _nb_cols = nc;
       _size = 0;
       _length = col_ind.size();
       _row_ptr.resize(_nb_rows+1);
       _row_ptr = row_ptr;
       _col_ind.resize(_length);
       _col_ind = col_ind;
       _a.resize(_length,T_(0.));
       _diag.resize(_size);
       _max_it = 1000;
       _toler = 1.e-8;
    }
#endif

#ifndef USE_EIGEN
/** \brief Constructor for a rectangle matrix
 *  @param [in] nr Number of rows
 *  @param [in] nc Number of columns
 *  @param [in] row_ptr Vector of row pointers (See the above description of this class).
 *  @param [in] col_ind Vector of column indices (See the above description of this class).\n
 *  @param [in] a vector instance containing matrix entries stored columnwise
 */
    SpMatrix(size_t                nr,
             size_t                nc,
             const vector<size_t>& row_ptr,
             const vector<size_t>& col_ind,
             const vector<T_>&     a)
             : _type(0), _dof(0), _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
    {
       _is_diagonal = false;
       _size = 0;
       _nb_rows = nr; _nb_cols = nc;
       _length = col_ind.size();
       _row_ptr.resize(_nb_rows+1);
       _row_ptr = row_ptr;
       _col_ind.resize(_length);
       _col_ind = col_ind;
       _a.resize(_length);
       _a = a;
       _diag.resize(_size);
       _max_it = 1000;
       _toler = 1.e-8;
    }
#endif

#ifndef USE_EIGEN
/** \brief Constructor for a rectangle matrix
 *  @param [in] row_ptr Vector of row pointers (See the above description of this class).
 *  @param [in] col_ind Vector of column indices (See the above description of this class).
 */
    SpMatrix(const vector<size_t>& row_ptr,
             const vector<size_t>& col_ind)
             : _type(0), _dof(0), _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
    {
       _is_diagonal = false;
       _size = _nb_rows = _nb_cols = row_ptr.size()-1;
       _length = col_ind.size();
       _row_ptr.resize(_size+1);
       _row_ptr = row_ptr;
       _col_ind.resize(_length);
       _col_ind = col_ind;
       _a.resize(_length,T_(0.));
       _diag.resize(_size);
       _max_it = 1000;
       _toler = 1.e-8;
    }
#endif

#ifndef USE_EIGEN
/** \brief Constructor for a rectangle matrix
 *  @param [in] row_ptr Vector of row pointers (See the above description of this class).
 *  @param [in] col_ind Vector of column indices (See the above description of this class).
 *  @param [in] a vector instance that contain matrix entries stored row by row.\n
 *  Number of rows is extracted from vector <tt>row_ptr</tt>.
 */
    SpMatrix(const vector<size_t>& row_ptr,
             const vector<size_t>& col_ind,
             const vector<T_>&     a)
             : _type(0), _dof(0), _is_dense(0), _extended(0), _solver(DIRECT_SOLVER), _prec(IDENT_PREC)
    {
       _is_diagonal = false;
       _size = _nb_rows = _nb_cols = row_ptr.size()-1;
       _length = col_ind.size();
       _row_ptr.resize(_size+1);
       _row_ptr = row_ptr;
       _col_ind.resize(_length);
       _col_ind = col_ind;
       _a.resize(_length);
       _a = a;
       _diag.resize(_size);
       _max_it = 1000;
       _toler = 1.e-8;
    }
#endif

/// \brief Copy constructor.
    SpMatrix(const SpMatrix& m)
    {
       _is_diagonal = m._is_diagonal;
       _dof = m._dof;
       _size = m._size;
       _is_dense = m._is_dense;
       _extended = m._extended;
#ifdef USE_EIGEN
       if (_nb_rows)
          _A.resize(_nb_rows,_nb_cols);
#else
       _length = m._length;
       _row_ptr.resize(_size+1);
       _col_ind.resize(_length);
       _col_ind = m._col_ind;
       _row_ptr = m._row_ptr;
       _a.resize(_length);
       _a = m._a;
#endif
       _diag.resize(_size);
       _diag = m._diag;
       _solver = -1;
       _prec = -1;
       _max_it = 1000;
       _toler = 1.e-8;
       _solver = -1;
       _prec = -1;
       _max_it = 1000;
       _toler = 1.e-8;
    }

/// \brief Destructor.
    ~SpMatrix() { }

/// \brief Define matrix as a dense one
    void Dense()
    {
       _is_dense = 1;
#ifdef USE_EIGEN
       _length = _nb_rows*_nb_cols;
       _A.reserve(_length);
       clear();
#endif
    }

/// \brief Define matrix as identity matrix
    void Identity()
    {
#ifdef USE_EIGEN
       _A.reserve(_length);
       for (size_t i=0; i<_nb_rows; i++)
          _A.coeffRef(i,i) = static_cast<T_>(1.);
#else
       for (size_t i=0; i<_nb_rows; ++i) {
          for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; ++j)
             _a[_row_ptr[i]+j-1] = static_cast<T_>(0.);
          _a[_row_ptr[i]+i-1] = static_cast<T_>(1.);
       }
#endif
    }

/// \brief Define matrix as a diagonal one
    void Diagonal()
    {
       _is_dense = 0;
       _is_diagonal = 1;
#ifdef USE_EIGEN
       _length = _nb_rows;
       _A.reserve(_length);
       for (size_t i=0; i<_nb_rows; i++)
          _A.insert(i,i) = 0;
#else
       for (size_t i=0; i<_length; i++)
          _a[i] = static_cast<T_>(0);
#endif
    }

/// \brief Define matrix as a diagonal one
/// with diagonal entries equal to <tt>a</tt>
    void Diagonal(const T_& a)
    {
#ifdef USE_EIGEN
       _A.reserve(_length);
       for (size_t i=0; i<_nb_rows; i++)
          _A.coeffRef(i,i) = a;
#else
       for (size_t i=0; i<_nb_rows; ++i) {
          for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; ++j)
             _a[_row_ptr[i]+j-1] = 0;
          _a[_row_ptr[i]+i-1] = a;
       }
#endif
    }

/** \brief Sets the matrix as the one for the Laplace equation in 1-D
 *  \details The matrix is initialized as the one resulting from P<sub>1</sub>
 *  finite element discretization of the classical elliptic operator
 *     -u'' = f
 *  with homogeneous Dirichlet boundary conditions
 *  \remark This function is available for real valued matrices only.
 *  @param [in] n Size of matrix (Number of rows)
 *  @param [in] h %Mesh size (assumed constant)
 */
    void Laplace1D(size_t n,
                   real_t h);


/** \brief Sets the matrix as the one for the Laplace equation in 2-D
 *  \details The matrix is initialized as the one resulting from P<sub>1</sub>
 *  finite element discretization of the classical elliptic operator
 *     -Delta u = f
 *  with homogeneous Dirichlet boundary conditions
 *  \remark This function is available for real valued matrices only.
 *  @param [in] nx Number of unknowns in the <tt>x</tt>-direction
 *  @param [in] ny Number of unknowns in the <tt>y</tt>-direction
 *  @remark The number of rows is equal to <tt>nx*ny</tt>
 */
    void Laplace2D(size_t nx,
                   size_t ny);

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
          if (_extended)
             _length = XGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
          else if (dof)
             _length = NodeGraphScal(mesh,_row_ptr,_col_ind,_IJ,_nbc);
          else
             _length = NodeGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
       }
       else if (_dof_type==SIDE_DOF)
          _length = SideGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
       else if (_dof_type==ELEMENT_DOF)
          _length = ElementGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
       else
          ;
#ifdef USE_EIGEN
       _A.resize(_size,_size);
       _A.reserve(_nbc);
       clear();
#else
       _a.resize(_length,static_cast<T_>(0));
#endif
       _diag.resize(_size);
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setMesh(size_t dof,
                 Mesh&  mesh,
                 int    code=0)
    {
       _size = _nb_rows = _nb_cols = mesh.getNbEq();
       if (dof)
          _size = _nb_rows = _nb_cols = mesh.getNbNodes();
       if (code!=0)
          _length = XGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
       else
          _length = NodeGraphScal(mesh,_row_ptr,_col_ind,_IJ,_nbc);
#ifdef USE_EIGEN
       _A.resize(_size,_size);
       _A.reserve(_nbc);
       clear();
#else
       _a.resize(_length,T_(0.));
#endif
       _diag.resize(_size);
    }

    void setMesh(size_t dof, 
                 size_t nb_eq,
                 Mesh&  mesh)
    {
       _type = 0;
       _dof = 0;
       _size = _nb_rows = _nb_cols = nb_eq;
       _length = NodeGraphScal(mesh,dof,nb_eq,_row_ptr,_col_ind,_IJ,_nbc);
#ifdef USE_EIGEN
       _A.resize(_size,_size);
       _A.reserve(_nbc);
       clear();
#else
       _a.resize(_length,T_(0.));
#endif
       _diag.resize(_size);
    }

    void setMesh(Mesh&  mesh,
                 size_t dof,
                 size_t type)
    {
       _dof_type = mesh.getDOFSupport();
       Matrix<T_>::init_set_mesh(mesh,dof);
       if (_dof_type==NODE_DOF) {
          if (type && dof)
             _length = XGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
          else if (dof)
             _length = NodeGraphScal(mesh,_row_ptr,_col_ind,_IJ,_nbc);
          else if (type==0 && dof==0)
             _length = NodeGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
       }
       else if (_dof_type==SIDE_DOF)
          _length = SideGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
       else if (_dof_type==ELEMENT_DOF)
          _length = ElementGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
       else
          ;
#ifdef USE_EIGEN
       _A.resize(_size,_size);
       _A.reserve(_nbc);
       clear();
#else
       _a.resize(_length,T_(0.));
#endif
       _diag.resize(_size);
    }

/** \brief Activate extended graph option.
 *  \details An extended graph is the one for which a node is linked not only to
 *  its neighbor nodes but also to the neighbors of the neighbors.
 */
    void setExtendedGraph() { _extended = true; }
   
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Activate 1-DOF per node option.
    void setOneDOF() { _one_dof = true; }

/// \brief Activate Sides option.
    void setSides() { _sides = true; }

/// \brief Store diagonal entries in a separate internal vector.
    void setDiag()
    {
#ifdef USE_EIGEN
       for (size_t i=0; i<_nb_rows; i++)
          _diag[i] = _A.coeff(i,i);
#else
       for (size_t i=0; i<_size; i++) {
          for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; j++) {
             if (i==j) {
                _diag[i] = _a[_row_ptr[i]+i-1];
                break;
             }
          }
       }
#endif
    }

/** \brief Impose by a diagonal method an essential boundary condition.
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  The penalty parameter is defined by default equal to 1.e20.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in] mesh Mesh instance from which information is extracted.
 *  @param [in,out] b Vect instance that contains right-hand side.
 *  @param [in] u Vect instance that conatins imposed valued at DOFs where they are to be imposed.
 */
    void DiagPrescribe(Mesh&           mesh,
                       Vect<T_>&       b,
                       const Vect<T_>& u)
    {
       real_t p = _zero;
       for (size_t j=1; j<=_nb_rows; j++)
          p = std::max(p,Abs(get(j,j)));
#ifdef USE_EIGEN
#else
       size_t k=0;
       mesh_nodes(mesh) {
          for (size_t i=1; i<=The_node.getNbDOF(); ++i) {
             size_t ii=The_node.getDOF(i)-1;
             for (size_t j=0; j<_row_ptr[ii+1]-_row_ptr[ii]; j++,k++)
                if (The_node.getCode(i)>0) {
                   b[ii] = p*u[ii];
                   _a[k] = 0;
                   if (ii+1==_col_ind[k])
                      _a[k] = p;
                }
          }
       }
#endif
    }

/** \brief Impose by a diagonal method an essential boundary condition using the Mesh instance
 *  provided by the constructor
 *  \details This member function modifies diagonal terms in matrix and terms
 *  in vector that correspond to degrees of freedom with nonzero code
 *  in order to impose a boundary condition.
 *  The penalty parameter is defined by default equal to 1.e20.
 *  It can be modified by member function \b setPenal(..).
 *  @param [in,out] b Vect instance that contains right-hand side.
 *  @param [in] u Vect instance that conatins imposed valued at DOFs where they are to be imposed.
 */
    void DiagPrescribe(Vect<T_>&       b,
                       const Vect<T_>& u)
    {
#ifdef USE_EIGEN
#else
       real_t p = 0;
       for (size_t j=1; j<=_nb_rows; j++)
          p = std::max(p,Abs(get(j,j)));
       size_t k=0;
       MESH_ND {
          for (size_t i=1; i<=TheNode.getNbDOF(); ++i) {
             size_t ii=TheNode.getDOF(i)-1;
             for (size_t j=0; j<_row_ptr[ii+1]-_row_ptr[ii]; j++,k++)
                if (TheNode.getCode(i)>0) {
                   b[ii] = p*u[ii];
                   _a[k] = 0;
                   if (ii+1==_col_ind[k])
                      _a[k] = p;
                }
          }
       }
#endif
    }

/// \brief Set size of matrix (case where it's a square matrix).
/// @param [in] size Number of rows and columns.
    void setSize(size_t size)
    {
       _nb_rows = _nb_cols = _size = size;
#ifdef USE_EIGEN
       _A.resize(size,size);
#else
       _length = lsize_t(_nb_rows*_nb_cols);
       _row_ptr.resize(_nb_rows+1);
       _col_ind.resize(_length);
       _row_ptr[0] = 0;
       for (size_t i=1; i<=_nb_rows; i++)
          _row_ptr[i] = _row_ptr[i-1] + _nb_cols;
       size_t l = 0;
       for (size_t i=0; i<_nb_rows; i++)
          for (size_t j=0; j<_nb_cols; j++)
             _col_ind[l++] = j+1;
       _a.resize(_length,T_(0.));
#endif
    }

/** \brief Set size (number of rows) of matrix
 *  @param [in] nr Number of rows
 *  @param [in] nc Number of columns
 */
    void setSize(size_t nr,
                 size_t nc)
    {
       _nb_rows = nr;
       _nb_cols = nc;
       _size = 0;
       if (_nb_rows==_nb_cols)
          _size = _nb_rows;
       _length = lsize_t(_nb_rows*_nb_cols);
#ifdef USE_EIGEN
       _A.resize(nr,nc);
#else
       _row_ptr.resize(_nb_rows+1);
       _col_ind.resize(_length);
       _row_ptr[0] = 1;
       for (size_t i=1; i<=_nb_rows; i++)
          _row_ptr[i] = _row_ptr[i-1] + _nb_cols;
       size_t l = 0;
       for (size_t i=0; i<_nb_rows; i++)
          for (size_t j=0; j<_nb_cols; j++)
             _col_ind[l++] = j + 1;
       _a.resize(_length,T_(0.));
#endif
    }

/** \brief Set graph of matrix by giving a vector of its nonzero entries
 *  @param [in] I Vector containing pairs of row and column indices
 *  @param [in] opt Flag indicating if vector <tt>I</tt> is cleaned and ordered
 *  (<tt>opt=1</tt>: default) or not (<tt>opt=0</tt>). In the latter case, this vector can have
 *  the same contents more than once and are not necessarily ordered
 */
    void setGraph(const Vect<RC>& I,
                  int             opt=1)
    {
       _length = I.size();
       _nb_rows = _nb_cols = 0;
       _IJ.resize(_length);
       for (size_t i=0; i<_length; i++) {
          _nb_rows = std::max(_nb_rows,I[i].first);
          _nb_cols = std::max(_nb_cols,I[i].second);
          _IJ[i] = RC(I[i].first-1,I[i].second-1);
       }
       _size = _nb_rows;
       if (opt==0) {
          sort(_IJ.begin(),_IJ.end());
          vector<RC>::iterator new_end = unique(_IJ.begin(),_IJ.end());
          _IJ.erase(new_end,_IJ.end());
       }
       _length = _IJ.size();
       _row_ptr.resize(_size+1);
       _col_ind.resize(_length);
       StoreGraph(_IJ,_row_ptr,_col_ind);
#ifdef USE_EIGEN
       for (size_t i=0; i<_row_ptr.size()-1; i++)
          _nbc.push_back(_row_ptr[i+1]-_row_ptr[i]);
       _A.reserve(_nbc);
       clear();
#else
       _a.resize(_length,0);
#endif
    }

/// \brief Get <tt>i</tt>-th row vector.
    Vect<T_> getRow(size_t i) const
    {
       Vect<T_> v(_nb_cols);
#ifdef USE_EIGEN
       size_t j=0;
       for (size_t k=0; k<_length; k++) {
          if (_IJ[k].first==i-1)
             v[j++] = _A.coeff(_IJ[k].first,_IJ[k].second);
       }
#else
       for (size_t j=1; j<=_nb_cols; j++)
          v(j) = get(i,j);
#endif
       return v;
    }

/// \brief Get <tt>j</tt>-th column vector.
    Vect<T_> getColumn(size_t j) const
    {
       Vect<T_> v(_nb_rows);
#ifdef USE_EIGEN
       size_t i=0;
       for (size_t k=0; k<_length; k++) {
          if (_IJ[k].second==j-1)
             v[i++] = _A.coeff(_IJ[k].first,_IJ[k].second);
       }
#else
       for (size_t i=1; i<=_nb_rows; i++)
          v(i) = get(i,j);
#endif
       return v;
    }

/** \brief Operator () (Non constant version)
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ & operator()(size_t i,
                    size_t j)
    {
#ifdef USE_EIGEN
       return _A.coeffRef(i-1,j-1);
#else
       int k=_col_index(i,j);
       if (k<0)
          throw OFELIException("In SpMatrix::set(i,j,x): Index pair: (" + itos(int(i)) +
                               "," + itos(int(j)) + ") is not compatible "
                               "with sparse storage.");
       else
          return _a[_row_ptr[i-1]+k];
       return _temp;
#endif
    }

/** \brief Operator <tt>()</tt> (Constant version)
 *  @param [in] i Row index
 *  @param [in] j Column index
 */
    T_ operator()(size_t i,
                  size_t j) const
    {
#ifdef USE_EIGEN
       return _A.coeff(i-1,j-1);
#else
       int k=_col_index(i,j);
       if (k<0)
          return _zero;
       else
          return _a[_row_ptr[i-1]+k];
#endif
    }

/** \brief Operator <tt>()</tt> with one argument (Constant version)
 *  \details Returns <tt>i</tt>-th position in the array storing matrix entries.
 *  The first entry is at location <tt>1</tt>.
 *  Entries are stored row by row.
 */
    const T_ operator()(size_t i) const { return _a[i-1]; }

/** \brief Operator <tt>[]</tt> (Constant version).
 *  \details Returns <tt>i</tt>-th position in the array storing matrix entries.
 *  The first entry is at location <tt>0</tt>.
 *  Entries are stored row by row.
 */
    const T_ operator[](size_t i) const { return _a[i]; }

/** \brief Operator <tt>*</tt> to multiply matrix by a vector
 *  @param [in] x Vect instance to multiply by
 *  @return Vector product of matrix by <tt>x</tt>
 */
    Vect<T_> operator*(const Vect<T_>& x) const
    {
       Vect<T_> y(_nb_rows);
       y.clear();
#ifdef USE_EIGEN
      for (size_t i=0; i<_length; i++)
         y[_IJ[i].first] += _A.coeff(_IJ[i].first,_IJ[i].second)*x[_IJ[i].second];
#else
      Mult(x,y);
#endif
      return y;
    }

/** \brief Operator <tt>*=</tt> to premultiply matrix by a constant
 *  @param [in] a Constant to multiply matrix by
 *  @return Resulting matrix
 */
    SpMatrix<T_>& operator*=(const T_& a)
    {
#ifdef USE_EIGEN
       for (size_t i=0; i<_length; i++)
          _A.coeffRef(_IJ[i].first,_IJ[i].second) *= a;
#else
       for (size_t k=0; k<_length; ++k)
          _a[k] *= a;
#endif
       return *this;
    }

/// \brief Get mesh instance whose reference will be stored in current instance of SpMatrix.
    void getMesh(Mesh& mesh)
    {
#if USE_EIGEN
#else
       if (_sides)
          SideGraph(mesh,_row_ptr,_col_ind);
       else {
          if (_dof)
             _size = _nb_rows = _nb_cols = mesh.getNbNodes();
          else
             _size = _nb_rows = _nb_cols = mesh.getNbEq();
          if (_type)
            XGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
//             XGraphScal(mesh,_row_ptr,_col_ind);
      else
        NodeGraph(mesh,_row_ptr,_col_ind,_IJ,_nbc);
   }
   _a.resize(_length,T_(0.));
#endif
   _diag.resize(_size);
}

/** \brief Multiply matrix by vector and save in another one.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void Mult(const Vect<T_>& x,
              Vect<T_>&       y) const
    {
#ifdef USE_EIGEN
       for (size_t i=0; i<_length; i++)
          y[_IJ[i].first] = _A.coeff(_IJ[i].first,_IJ[i].second)*x[_IJ[i].second];
#else
       y = static_cast<T_>(0);
       MultAdd(x,y);
#endif
    }

/** \brief Multiply matrix by vector <tt>x</tt> and add to <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(const Vect<T_>& x,
                 Vect<T_>&       y) const
    {
#ifdef USE_EIGEN
       for (size_t i=0; i<_length; i++)
          y[_IJ[i].first] += _A.coeff(_IJ[i].first,_IJ[i].second)*x[_IJ[i].second];
#else
       size_t l=0;
       for (size_t i=0; i<_nb_rows; ++i)
          for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; ++j)
             y[i] += _a[_row_ptr[i]+j] * x[_col_ind[l++]-1];
#endif
    }

/** \brief Multiply matrix by vector <tt>a*x</tt> and add to <tt>y</tt>.
 *  @param [in] a Constant to multiply by matrix
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector to add to the result. <tt>y</tt> contains on output the result.
 */
    void MultAdd(T_              a,
                 const Vect<T_>& x,
                 Vect<T_>&       y) const
{
#ifdef USE_EIGEN
       for (size_t i=0; i<_length; i++)
          y[_IJ[i].first] += a*_A.coeff(_IJ[i].first,_IJ[i].second)*x[_IJ[i].second];
#else
       size_t l=0;
       for (size_t i=0; i<_nb_rows; i++)
          for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; j++)
             y[i] += a * _a[_row_ptr[i]+j] * x[_col_ind[l++]-1];
#endif
    }

/** \brief Multiply transpose of matrix by vector <tt>x</tt> and save in <tt>y</tt>.
 *  @param [in] x Vector to multiply by matrix
 *  @param [out] y Vector that contains on output the result.
 */
    void TMult(const Vect<T_>& x,
               Vect<T_>&       y) const
    {
#ifdef USE_EIGEN
       for (size_t i=0; i<_length; i++)
          y[_IJ[i].first] += _A.coeff(_IJ[i].second,_IJ[i].first)*x[_IJ[i].second];
#else
       for (size_t i=0; i<_nb_rows; i++)
          for (size_t j=0; j<_row_ptr[i+1]-_row_ptr[i]; j++)
             y[_col_ind[j+_row_ptr[i]-1]] += _a[_row_ptr[i]+j] * x[i];
#endif
    }

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m %Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                  a,
              const SpMatrix<T_>& m)
    {
#ifdef USE_EIGEN
       _A += a * m._A;
#else
       _a += a * m._a;
#endif
    }

/** \brief Add to matrix the product of a matrix by a scalar
 *  @param [in] a Scalar to premultiply
 *  @param [in] m Pointer to Matrix by which <tt>a</tt> is multiplied. The result is added
 *  to current instance
 */
    void Axpy(T_                a,
              const Matrix<T_>* m) { }

/** \brief Assign a value to an entry of the matrix
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] val Value to assign to <tt>a(i,j)</tt>
 */
    void set(size_t    i,
             size_t    j,
             const T_& val)
    {
#ifdef USE_EIGEN
       _A.coeffRef(i-1,j-1) = val;
#else
       int k=_col_index(i,j);
       if (k<0)
          throw OFELIException("In SpMatrix::set(i,j,x): Index pair (" + itos(int(i)) +
                               "," + itos(int(j)) + ") is not compatible "
                               "with sparse storage.");
       else
          _a[_row_ptr[i-1]+k] = val;
#endif
    }

/** \brief Add a value to an entry of the matrix
 *  @param [in] i Row index
 *  @param [in] j Column index
 *  @param [in] val Constant value to add to <tt>a(i,j)</tt>
 */
    void add(size_t    i,
             size_t    j,
             const T_& val)
    {
#ifdef USE_EIGEN
       _A.coeffRef(i-1,j-1) += val;
#else
       _a[_row_ptr[i-1]+_col_index(i,j)] += val;
#endif
    }

/// \brief Operator =.
/// \details Assign constant value <tt>x</tt> to all matrix entries.
    void operator=(const T_& x)
    {
#ifdef USE_EIGEN
       for (size_t i=0; i<_length; i++) {
          _A.coeffRef(_IJ[i].first,_IJ[i].second) = static_cast<T_>(0);
          if (_IJ[i].first==_IJ[i].second)
             _A.coeffRef(_IJ[i].first,_IJ[i].first) = x;	
    }
#else
    for (size_t i=0; i<_length; i++)
       _a[i] = x;
#endif
}

/// \brief Return storage information.
/// @return Column index of the <tt>i</tt>-th stored element in matrix
    size_t getColInd(size_t i) const { return _col_ind[i-1]; }

/// \brief Return Row pointer at position <tt>i</tt>.
    size_t getRowPtr(size_t i) const { return _row_ptr[i-1]; }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Factorize matrix
    int Factor() { return -1; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Perform a diagonal incomplete LU factorization of matrix
 *  @param [out] id
 *  @param [out] pivot
 *  \return Return <tt>0</tt> if the factorization was normally achieved,
 *  <tt>n</tt> if the <tt>n</tt>-th pivot is null.
 */
#ifdef USE_EIGEN
    int DILUFactorize(vector<size_t>& id,
                      vector<T_>&     pivot) const;
#else
    int DILUFactorize(vector<size_t>& id,
                      vector<T_>&     pivot) const
    {
       id.resize(_size);
       pivot.resize(_size);
       size_t k=0;
       for (size_t i=0; i<_size; i++) {
          for (size_t j=_row_ptr[i]; j<_row_ptr[i+1]; j++, k++) {
             if (_col_ind[j]==i+1) {
                id[i] = k + 1;
                if (_a[k]==static_cast<T_>(0))
                   throw OFELIException("In SpMatrix::DILUFactorize(...): Zero pivot detected in row "
                                        +itos(i+1));
                pivot[i] = _a[k];
             }
          }
       }
       int found=0;
       T_ c=static_cast<T_>(0);
       for (size_t i=0; i<_size; ++i) {
          pivot[i] = static_cast<T_>(1)/pivot[i];
          for (size_t j=id[i]; j<_row_ptr[i+1]; ++j) {
             found = 0;
             size_t l=_col_ind[j]-1;
             for (k=_row_ptr[l]; k<id[l]-1; ++k) {
                if (_col_ind[k]==i+1)
                   found = 1, c = _a[k];
             }
             if (found)
                pivot[l] -= c*pivot[i]*_a[j];
          }
       }
       return 0;
    }
#endif

/** \brief Perform an Incomplete LU factorization of matrix
 *  \return Return <tt>0</tt> if the factorization was normally achieved,
 *  <tt>n</tt> if the <tt>n</tt>-th pivot is null.
 */
#ifdef USE_EIGEN
    int ILUFactorize(const SpMatrix<T_>& A);
#else
    int ILUFactorize(const SpMatrix<T_>& A)
    {
       _nb_rows = A._nb_rows, _nb_cols = A._nb_cols;
       _size = A._size, _length = A._length;
       _row_ptr = A._row_ptr; _col_ind = A._col_ind;
       _lnnz = (_length-_nb_rows)/2, _unnz = (_length+_nb_rows)/2;
       _aL.resize(_lnnz); _aU.resize(_unnz);
       _l_col_ind.resize(_lnnz);
       _u_col_ind.resize(_unnz);
       _l_row_ptr.resize(_nb_rows+1);
       _u_row_ptr.resize(_nb_rows+1);
       _l_row_ptr[0] = _u_row_ptr[0] = 0;

       for (size_t i=0; i<_nb_rows; i++) {
          _l_row_ptr[i+1] = _l_row_ptr[i];
          _u_row_ptr[i+1] = _u_row_ptr[i];
          for (size_t j=A._row_ptr[i]; j<A._row_ptr[i+1]; j++) {
             if (A._col_ind[j]<i+1) {
                size_t k=_l_row_ptr[i+1]++;
                _aL[k] = A._a[j];
                _l_col_ind[k] = A._col_ind[j];
             }
             else {
                size_t k=_u_row_ptr[i+1]++;
                _aU[k] = A._a[j];
                _u_col_ind[k] = A._col_ind[j];
             }
          }
       }
       for (size_t i=1; i<_nb_rows; i++) {
          for (size_t j=_l_row_ptr[i]; j<_l_row_ptr[i+1]; j++) {
             size_t pn=_u_row_ptr[_col_ind[j]-1], qn=j+1, rn=_u_row_ptr[i];
             T_ p = (_aL[j]/=_aU[pn]);
             for (pn++; pn<_u_row_ptr[_l_col_ind[j]] && _u_col_ind[pn]<i+1; pn++) {
                while (qn<_l_row_ptr[i+1] && _l_col_ind[qn]<_u_col_ind[pn])
                   qn++;
                if (qn<_l_row_ptr[i+1] && _u_col_ind[pn]==_l_col_ind[qn])
                   _aL[qn] -= p*_aU[pn];
             }
             for (; pn<_u_row_ptr[_l_col_ind[j]]; pn++) {
                while (rn<_u_row_ptr[i+1] && _u_col_ind[rn]<_u_col_ind[pn])
                   rn++;
                if (rn<_u_row_ptr[i+1] && _u_col_ind[pn]==_u_col_ind[rn])
                   _aU[rn] -= p*_aU[pn];
             }
          }
       }

       size_t k=0, kl=0, ku=0;
       _a.resize(_length);
       for (size_t i=0; i<_nb_rows; i++) {
          for (size_t j=_row_ptr[i]; j<_row_ptr[i+1]; j++) {
             if (_col_ind[j]<i+1)
                _a[k++] = _aL[kl++];
             else
                _a[k++] = _aU[ku++];
          }
       }
       return 0;
    }
#endif

/** \brief Solve a linear system with an diagonal incompletely factorized matrix
 *  @param [in] b Vect instance containing the right-hand side
 *  @param [out] x Vect instance containing on output the solution
 */
#ifndef USE_EIGEN
    void DILUSolve(const vector<size_t>& id,
                   const vector<T_>&     pivot,
                   const Vect<T_>&       b,
                   Vect<T_>&             x) const
    {
       vector<T_> z(_size);
       for (size_t i=0; i<_size; i++) {
          T_ s = 0;
          for (size_t j=_row_ptr[i]; j<id[i]; ++j)
             s += _a[j] * z[_col_ind[j]-1];
          z[i] = pivot[i] * (b[i]-s);
       }
       for (size_t i=0; i<_size; ++i) {
          T_ s = 0;
          for (size_t j=id[_size-i-1]; j<_row_ptr[_size-i]; ++j)
             s += _a[j] * x(_col_ind[j]);
          x[_size-i-1] = z[_size-i-1] - pivot[_size-i-1] * s;
       }
    }
#endif

/** \brief Solve a linear system with an incompletely factorized matrix
 *  @param [in] b Vect instance containing the right-hand side
 *  @param [out] x Vect instance containing on output the solution
 */
#ifndef USE_EIGEN
    void ILUSolve(const Vect<T_>& b,
                  Vect<T_>&       x) const
    {
       for (size_t i=0; i<_size; i++) {
          T_ s=0;
          for (size_t j=_l_row_ptr[i]; j<_l_row_ptr[i+1]; j++)
             s += _aL[_l_col_ind[j]-1]*b[_l_col_ind[j]-1];
          x[i] = b[i] - s;
       }
       for (int i=int(_size)-1; i>=0; i--) {
          T_ s=0;
          for (size_t j=_l_row_ptr[i]; j<_l_row_ptr[i+1]; j++)
             s += _aU[_u_col_ind[j]-1]*x[_u_col_ind[j]-1];
          x[i] = (x[i]-s)/_aU[_u_col_ind[_l_row_ptr[i]]-1];
       }
    }
#endif

/** \brief Solve a linear system with an incompletely factorized matrix
 *  @param [in] b Vect instance containing the right-hand side
 *  @param [out] x Vect instance containing on output the solution
 */
#ifndef USE_EIGEN
    void SSORSolve(const Vect<T_>& b,
                   Vect<T_>&       x) const
    {
       size_t k=0;
       vector<size_t> id(_size);
       for (size_t i=0; i<_size; i++) {
          for (size_t j=_row_ptr[i]; j<_row_ptr[i+1]; j++, k++) {
             if (_col_ind[j]==i+1) {
                id[i] = k + 1;
                if (_a[k]==static_cast<T_>(0))
                   throw OFELIException("In SpMatrix::setMatrix(SpMatrix<T_>):"
                                        " Zero pivot detected in row " + itos(i+1));
             }
          }
       }

       Vect<T_> z(_size);
       for (size_t i=0; i<_size; i++) {
          T_ s = 0;
          for (size_t j=_row_ptr[i]; j<id[i]-1; j++)
             s += _a[j] * z(_col_ind[j]);
          z[i] = (b[i]-s)/_a[id[i]-1];
       }
       for (int i=int(_size)-1; i>=0; i--) {
          T_ s = 0;
          for (size_t j=id[i]; j<_row_ptr[i+1]; j++)
             s += _a[j] * x(_col_ind[j]);
          x[i] = z[i] - s/_a[id[i]-1];
       }
    }
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int solve(Vect<T_>& b)
    {
       Vect<T_> x(b.size());
       int ret = solve(b,x);
       b = x;
       return ret;
    }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \brief Solve the linear system of equations.
 *  \details The default parameters are:
 *  <ul>
 *     <li><tt>CG_SOLVER</tt> for solver
 *     <li><tt>DIAG_PREC</tt> for preconditioner
 *     <li>Max. Number of iterations is <tt>1000</tt>
 *     <li>Tolerance is <tt>1.e-8</tt>
 *  </ul>
 *  To change these values, call function setSolver before this function
 *  @param [in] b Vector that contains right-hand side
 *  @param [out] x Vector that contains the obtained solution
 *  @return Number of actual performed iterations 
 */
    int solve(const Vect<T_>& b,
              Vect<T_>&       x)
    {
       if (_solver==DIRECT_SOLVER)
          throw OFELIException("In SpMatrix::solve(...): No solver provided.");
       LinearSolver<T_> ls(*this,b,x);
       int ret = ls.solve(_solver,_prec);
       return ret;
    }

/** \brief Choose solver and preconditioner for an iterative procedure
 *  @param [in] solver Option to choose iterative solver in an enumerated variable
 *  <ul>
 *     <li><tt>CG_SOLVER</tt>: Conjugate Gradient [default]
 *     <li><tt>CGS_SOLVER</tt>: Squared conjugate gradient
 *     <li><tt>BICG_SOLVER</tt>: Biconjugate gradient
 *     <li><tt>BICG_STAB_SOLVER</tt>: Biconjugate gradient stabilized
 *     <li><tt>GMRES_SOLVER</tt>: Generalized Minimal Residual
 *  </ul>
 *  Default value is <tt>CG_SOLVER</tt>
 *  @param [in] prec Option to choose preconditioner in an enumerated variable
 *  <ul>
 *     <li><tt>IDENT_PREC</tt>: Identity preconditioner (no preconditioning)
 *     <li><tt>DIAG_PREC</tt>: Diagonal preconditioner [default]
 *     <li><tt>SSOR_PREC</tt>: SSOR (Symmetric Successive Over Relaxation) preconditioner
 *     <li><tt>DILU_PREC</tt>: ILU (Diagonal Incomplete factorization) preconditioner
 *     <li><tt>ILU_PREC</tt>: ILU (Incomplete factorization) preconditioner
 *  </ul>
 *  Default value is <tt>DIAG_PREC</tt>
 *  @param [in] max_it Maximum number of allowed iterations. Default value is <tt>1000</tt>.
 *  @param [in] toler Tolerance for convergence. Default value is <tt>1.e-8</tt>
 */
    void setSolver(Iteration      solver=CG_SOLVER,
                   Preconditioner prec=DIAG_PREC,
                   int            max_it=1000,
                   real_t         toler=1.e-8)
    {
       _solver = solver;
       _prec = prec;
       _max_it = max_it;
       _toler = toler;
    }

/// brief Set all matrix entries to zero
    void clear()
    {
#ifdef USE_EIGEN
       for (size_t i=0; i<_length; i++)
          _A.coeffRef(_IJ[i].first,_IJ[i].second) = static_cast<T_>(0);
#endif
}

/// \brief Return C-Array.
/// \details Non zero terms of matrix is stored row by row.
    T_ *get() const { return &_a[0]; }

/** \brief  Return entry <tt>(i,j)</tt> of matrix if this one is stored, <tt>0</tt> otherwise
 *  @param [in] i Row index (Starting from 1)
 *  @param [in] j Column index (Starting from 1)
 */
    T_ get(size_t i,
           size_t j) const
    {
#ifdef USE_EIGEN
       return _A.coeff(i-1,j-1);
#else
       int k=_col_index(i,j);
       if (k<0)
          return _zero;
       else
          return _a[_row_ptr[i-1]+k-1];
#endif
}

#ifdef USE_EIGEN
/// \brief Return reference to the matrix instance in Eigen library
    SpMat& getEigenMatrix() { return _A; }
#endif

#ifdef USE_EIGEN
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#endif

/** \fn ostream & operator<<(ostream& s, const SpMatrix<T_> &A)
 *  \ingroup VectMat
 *  Output matrix in output stream
 */
    template<class TT_>
    friend ostream& operator<<(ostream& s, const SpMatrix<TT_>& A);

 private:

   size_t _dof;
   int _is_dense, _type, _max_it;
   int _one_dof, _sides, _extended;
   real_t _toler;
   vector<RC> _IJ;
   vector<size_t> _nbc;
   Iteration _solver;
   Preconditioner _prec;
#ifdef USE_EIGEN
   SpMat _A;
#else
   vector<T_> _pivot, _aL, _aU;
   vector<size_t> _id, _l_row_ptr, _u_row_ptr, _l_col_ind, _u_col_ind;
   size_t _lnnz, _unnz;
   int _col_index(size_t i, size_t j) const
   {
      for (int k=0; k<int(_row_ptr[i]-_row_ptr[i-1]); ++k)
         if (_col_ind[_row_ptr[i-1]+k]==j)
            return k;
      return -1;
   }
#endif
};


///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<>
inline void SpMatrix<real_t>::Laplace1D(size_t n,
                                        real_t h)
{
   _is_dense = 0;
   _is_diagonal = 0;
   _size = _nb_rows = _nb_cols = n;
#ifdef USE_EIGEN
   _length = 3*_size - 2;
   _A.reserve(_length);
   _A.insert(0,0) =  2./h;
   _A.insert(0,1) = -1./h;
   for (size_t i=1; i<_nb_rows-1; i++) {
      _A.insert(i,i-1) = -1./h;
      _A.insert(i,i  ) =  2./h;
      _A.insert(i,i+1) = -1./h;
   }
   _A.insert(_nb_rows-1,_nb_rows-2) = -1./h;
   _A.insert(_nb_rows-1,_nb_rows-1) =  2./h;
#else
   _row_ptr.push_back(0);
   _row_ptr.push_back(2);
   _col_ind.push_back(1);
   _col_ind.push_back(2);
   _a.push_back( 2./h);
   _a.push_back(-1./h);
   for (size_t i=1; i<_size-1; i++) {
      _row_ptr.push_back(_row_ptr[i]+3);
      _col_ind.push_back(i);
      _col_ind.push_back(i+1);
      _col_ind.push_back(i+2);
      _a.push_back(-1./h);
      _a.push_back( 2./h);
      _a.push_back(-1./h);
   }
   _col_ind.push_back(_size-1);
   _col_ind.push_back(_size);
   _row_ptr.push_back(_row_ptr[_size-1]+2);
   _a.push_back(-1./h);
   _a.push_back( 2./h);
   _length = _row_ptr[_size];
#endif
}


template<>
inline void SpMatrix<real_t>::Laplace2D(size_t nx,
                                        size_t ny)
{
   _is_dense = 0;
   _is_diagonal = 0;
   _size = _nb_rows = _nb_cols = nx*ny;
#ifdef USE_EIGEN
   _length = 5*_size;
   _A.reserve(_length);
   for (size_t ii=0; ii<_size; ii++) {
      _A.insert(ii,ii) =  4.;
      size_t i=ii/ny, j=ii-i*ny;
      if (i>0)
         _A.insert(ii,ii-ny) = -1.;
      if (i<nx-1)
         _A.insert(ii,ii+ny) = -1.;
      if (j>0)
         _A.insert(ii,ii- 1) = -1.;
      if (j<ny-1)
         _A.insert(ii,ii+ 1) = -1.;
   }
#else
   _row_ptr.push_back(0);
   _length = 0;
   for (size_t ii=0; ii<_size; ii++) {
      size_t i=ii/ny, j=ii-i*ny;
      if (i>0) {
         _a.push_back(-1.);
         _col_ind.push_back(ii-ny+1);
         _length++;
      }
      if (j>0) {
         _a.push_back(-1.);
         _col_ind.push_back(ii);
         _length++;
      }
      _a.push_back(4.);
      _col_ind.push_back(ii+1);
      _length++;
      if (j<ny-1) {
         _a.push_back(-1.);
         _col_ind.push_back(ii+2);
         _length++;
      }
      if (i<nx-1) {
         _a.push_back(-1.);
         _col_ind.push_back(ii+ny+1);
         _length++;
      }
      _row_ptr.push_back(_length);
   }
#endif
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \fn Vect<T_> operator*(const SpMatrix<T_> &A, const Vect<T_> &b)
 *  \brief Operator * (Multiply vector by matrix and return resulting vector
 *  \ingroup VectMat
 *  @param [in] A SpMatrix instance to multiply by vector
 *  @param [in] b Vect instance 
 *  \return Vect instance containing <tt>A*b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
Vect<T_> operator*(const SpMatrix<T_>& A,
                   const Vect<T_>&     b)
{
   Vect<T_> v(b.size());
   A.Mult(b,v);
}


/** \fn ostream & operator<<(ostream& s, const SpMatrix<T_> &A)
 *  \ingroup VectMat
 *  \brief Output matrix in output stream
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template<class T_>
ostream& operator<<(ostream&            s,
                    const SpMatrix<T_>& A)
{
#ifdef USE_EIGEN
   s << endl << A._A << endl;
#else
   s.setf(ios::right|ios::scientific);
   s << endl;
   size_t k = 0;
   for (size_t i=0; i<A._nb_rows; ++i) {
      for (size_t j=A._row_ptr[i]; j<A._row_ptr[i+1]; ++j)
         s << "(" << setw(6) << i+1 << "," << setw(6) << A._col_ind[j] << "): "
           << setprecision(8) << std::setfill(' ') << setw(18) << A._a[k++] << endl;
   }
#endif
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
