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

                    Definition of abstract class 'Equation'
                          for finite element equations

  ==============================================================================*/


#ifndef __EQUATION_H
#define __EQUATION_H

#include "equations/Equa.h"
#include "linear_algebra/LocalMatrix.h"
#include "linear_algebra/LocalVect.h"


namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Equation.h
 *  \brief Definition file for class Equation.
 */

class Mesh;
class Side;
class Node;

/*! \class Equation
 *  \ingroup Equation
 * \brief Abstract class for all equation classes.
 *
 * Template Arguments:
 *
 * \arg \b T_   : data type (real_t, float, ...)
 * \arg \b NEN_ : Number of element nodes
 * \arg \b NEE_ : Number of element equations
 * \arg \b NSN_ : Number of side nodes
 * \arg \b NSN_ : Number of side equations
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
class Equation : virtual public Equa<T_>
{

 public:

  using Equa<T_>::_theMesh;
  using Equa<T_>::_field_type;
  using Equa<T_>::_matrix_type;
  using Equa<T_>::_TimeInt;
  using Equa<T_>::_nb_fields;
  using Equa<T_>::_nb_nodes;
  using Equa<T_>::_nb_el;
  using Equa<T_>::_nb_dof;
  using Equa<T_>::_nb_dof_total;
  using Equa<T_>::_nb_eq;
  using Equa<T_>::_theFct;

/// Default constructor.
/// Constructs an "empty" equation
    Equation();

/// \brief Constructor with mesh instance.
/// @param [in] mesh Mesh instance
    Equation(Mesh &mesh);

/** \brief Constructor with mesh instance and solution vector.
 *  @param [in] mesh Mesh instance
 *  @param [in] u Vect instance containing solution.
 */
    Equation(Mesh&     mesh,
             Vect<T_>& u);

/** \brief Constructor with mesh instance, matrix and right-hand side.
 *  @param [in] mesh Mesh instance
 *  @param [in] u Vect instance containing Right-hand side.
 *  @param [in] init_time Initial Time value
 *  @param [in] final_time Final Time value
 *  @param [in] time_step Time step value
 */
    Equation(Mesh&     mesh,
             Vect<T_>& u,
             real_t&   init_time,
             real_t&   final_time,
             real_t&   time_step);

/// \brief Destructor
    virtual ~Equation();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Initialize equation.
 *  @param [in] mesh Mesh instance
 *  @param [in] ts Time step value
 *  @param [in] t Initial time [Default: <tt>0.</tt>]
 */
    void initEquation(Mesh&  mesh,
                      real_t init_time,
                      real_t final_time,
                      real_t time_step);

/// \brief Return problem matrix: Case of a SpMatrix (sparse).
/// @param [in] A Reference to matrix
    void getMatrix(const SpMatrix<T_>& A) const;

/// \brief Return problem matrix: Case of a SkMatrix (skyline)
/// @param [in] A Reference to matrix
    void getMatrix(const SkMatrix<T_>& A) const;

/// \brief Return problem matrix: Case of a SkSMatrix (symmetric skyline)
/// @param [in] A Reference to matrix
    void getMatrix(const SkSMatrix<T_>& A) const;

/// \brief Return problem matrix (Case of a SpMatrix (sparse))
/// @param [in] A Reference to matrix
    void getSolution(const Vect<T_>& u) const;

/// \brief Return problem right-hand side vector
    void getRHS(const Vect<T_>& b) const;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Update Right-Hand side by taking into account essential boundary conditions
 *  @param [in] el Reference to current element instance
 *  @param [in] bc Vector that contains imposed values at all DOFs
 */
    void updateBC(const Element&  el,
                  const Vect<T_>& bc);

/** \brief Update element matrix to impose bc by diagonalization technique
 *  @param [in] dof_type DOF type option. To choose among the enumerated values:
 *  <ul>
 *    <li><tt>NODE_DOF</tt>, DOFs are supported by nodes [Default]
 *    <li><tt>ELEMENT_DOF</tt>, DOFs are supported by elements
 *    <li><tt>SIDE_DOF</tt>, DOFs are supported by sides
 *  </ul>
 *  @param [in] dof DOF setting:
 *  <ul>
 *     <li><tt> = 0</tt>, All DOFs are taken into account [Default]
 *     <li><tt>!= 0</tt>, Only DOF No. <tt>dof</tt> is handled in the system
 *  </ul>
 */
    void DiagBC(DOFSupport dof_type=NODE_DOF,
                int        dof=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Set Body Force
    void setBodyForce(const Vect<T_>& f);

/// \brief Calculate residue in element
    void setResidue();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Localize Element Vector from a Vect instance.
 *  @param [in] b Reference to global vector to be localized.
 *  The resulting local vector can be accessed by attribute ePrev.
 *  This member function is to be used if a constructor with Element was invoked.
 */
    void LocalNodeVector(Vect<T_>& b);

/** \brief Localize Element Vector from a Vect instance.
 *  @param [in] b Global vector to be localized.
 *  @param [out] be Local vector, the length of which is the total number of element equations. 
 *  @remark All degrees of freedom are transferred to the local vector
 */
    void ElementNodeVector(const Vect<T_>&     b,
                           LocalVect<T_,NEE_>& be);

/** \brief Localize Side Vector from a Vect instance.
 *  @param [in] b Global vector to be localized.
 *  @param [out] bs Local vector, the length of which is the total number of side equations. 
 *  @remark All degrees of freedom are transferred to the local vector
 */
    void SideNodeVector(const Vect<T_>&     b,
                        LocalVect<T_,NSE_>& bs);

/** \brief Localize Side Vector from a Vect instance.
 *  @param [in] b Global vector to be localized.
 *  @param [out] bs Local constant value of vector at given side. 
 *  @remark All degrees of freedom are transferred to the local vector
 */
    void SideSideVector(const Vect<T_>& b,
                        T_*             bs);

/** \brief Localize Element Vector from a Vect instance.
 *  @param [in] b Global vector to be localized.
 *  @param [out] be Local vector, the length of which is the total number of element equations. 
 *  @remark Vector <tt>b</tt> is assumed to contain only one degree of freedom by node.
 */
    void ElementNodeVectorSingleDOF(const Vect<T_>&     b,
                                    LocalVect<T_,NEN_>& be);

/** \brief Localize Element Vector from a Vect instance.
 *  @param [in] b Global vector to be localized.
 *  @param [out] be Local vector, the length of which is the total number of element equations. 
 *  @param [in] dof Degree of freedom to transfer to the local vector
 *  @remark Only yhe dega dof is transferred to the local vector
 */
    void ElementNodeVector(const Vect<T_>&     b,
                           LocalVect<T_,NEN_>& be,
                           int                 dof);

/** \brief Localize Element Vector from a Vect instance.
 *  @param [in] b Global vector to be localized.
 *  @param [out] be Local vector, the length of which is 
 */
    void ElementSideVector(const Vect<T_>&     b,
                           LocalVect<T_,NSE_>& be);

/** \brief Localize Element Vector.
 *  @param [in] b Global vector to be localized
 *  @param [in] dof_type DOF type option. To choose among the enumerated values:
 *  <ul>
 *    <li><tt>NODE_DOF</tt>, DOFs are supported by nodes [Default]
 *    <li><tt>ELEMENT_DOF</tt>, DOFs are supported by elements
 *    <li><tt>SIDE_DOF</tt>, DOFs are supported by sides
 *  </ul>
 *  @param [in] flag Option to set:
 *  <ul>
 *    <li><tt> = 0</tt>, All DOFs are taken into account [Default]
 *    <li><tt>!= 0</tt>, Only DOF number <tt>dof</tt> is handled in the system
 *  </ul>
 *  The resulting local vector can be accessed by attribute <tt>ePrev</tt>.
 *  @remark This member function is to be used if a constructor with Element was invoked.
 *  It uses the Element pointer <tt>_theElement</tt>
 */
    void ElementVector(const Vect<T_>& b,
                       DOFSupport      dof_type=NODE_DOF,
                       int             flag=0);

/** \brief Localize Side Vector.
 *  @param [in] b Global vector to be localized
 *  <ul>
 *     <li><tt>NODE_DOF</tt>, DOFs are supported by nodes [ default ]
 *     <li><tt>ELEMENT_DOF</tt>, DOFs are supported by elements
 *     <li><tt>SIDE_DOF</tt>, DOFs are supported by sides
 *     <li><tt>BOUNDARY_SIDE_DOF</tt>, DOFs are supported by boundary sides
 *  </ul>
 *  @param [out] sb Array in which local vector is stored
 *  The resulting local vector can be accessed by attribute <tt>ePrev</tt>.
 *  @remark This member function is to be used if a constructor with Side was invoked.
 *  It uses the Side pointer <tt>_theSide</tt>
 */
    void SideVector(const Vect<T_>& b,
                    T_*             sb);

/** \brief Localize coordinates of element nodes
 *  \details Coordinates are stored in array <tt>_x[0], _x[1], ...</tt> which
 *  are instances of class Point<real_t>
 *  @remark This member function uses the Side pointer <tt>_theSide</tt>
 */
    void ElementNodeCoordinates();

/** \brief Localize coordinates of side nodes
 *  \details Coordinates are stored in array <tt>_x[0], _x[1], ...</tt> which
 *  are instances of class Point<real_t>
 *  @remark This member function uses the Element pointer <tt>_theElement</tt>
 */
    void SideNodeCoordinates();

/** \brief Assemble element matrix into global one
 *  @param A Pointer to global matrix (abstract class: can be any of classes SkSMatrix, SkMatrix, SpMatrix)
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(Matrix<T_>* A);

#if defined (USE_PETSC)
/** \brief Assemble element matrix into global one
 *  @param A Reference to global matrix
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(PETScMatrix<T_>& A);

/** \brief Assemble side matrix into global one
 *  @param A Reference to global matrix
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(PETScMatrix<T_>& A);

/** \brief Assemble element right-hand side vector into global one
 *  @param b Reference to global right-hand side vector
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(PETScVect<T_>& b);

/** \brief Assemble side right-hand side vector into global one
 *  @param b Reference to global right-hand side vector
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(PETScVect<T_>& b);
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void ElementAssembly(Element& el, Matrix<T_>* A);
    void ElementAssembly(Element& el, SkMatrix<T_>& A);
    void ElementAssembly(Element& el, SkSMatrix<T_>& A);
    void ElementAssembly(Element& el, SpMatrix<T_>& A);
    void ElementAssembly(Element& el, BMatrix<T_>& A);
    void ElementAssembly(Element& el, TrMatrix<T_>& A);
    void SideAssembly(Side& sd, Matrix<T_>* A);
    void SideAssembly(Side& sd, BMatrix<T_>& A);
    void SideAssembly(Side& sd, SkMatrix<T_>& A);
    void SideAssembly(Side& sd, SkSMatrix<T_>& A);
    void SideAssembly(Side& sd, SpMatrix<T_>& A);
    void ElementAssembly(Element& el, Vect<T_>& v);
    void SideAssembly(Side& sd, Vect<T_>& v);

#if defined(USE_PETSC)
    void updateBC(const Element& el, const PETScVect<T_>& bc);
    void updateBC(const PETScVect<T_>& bc);
    void ElementAssembly(Element& el, PETScMatrix<T_>& A);
    void SideAssembly(Side& sd, PETScMatrix<T_>& A);
    void ElementAssembly(Element& el, PETScVect<T_>& v);
    void SideAssembly(Side& sd, PETScVect<T_>& v);
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Assemble element matrix into global one
 *  @param A Global matrix stored as a BMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(BMatrix<T_>& A);

/** \brief Assemble element matrix into global one
 *  @param A Global matrix stored as an SkSMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(SkSMatrix<T_>& A);

/** \brief Assemble element matrix into global one
 *  @param [in] A Global matrix stored as an SkMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(SkMatrix<T_>& A);

/** \brief Assemble element matrix into global one
 *  @param [in] A Global matrix stored as an SpMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(SpMatrix<T_>& A);

/** \brief Assemble element matrix into global one
 *  @param [in] A Global matrix stored as an TrMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(TrMatrix<T_>& A);

/** \brief Assemble element matrix into global one for the Discontinuous Galerkin approximation
 *  @param A Pointer to global matrix (abstract class: can be any of classes SkSMatrix, SkMatrix, SpMatrix)
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void DGElementAssembly(Matrix<T_>* A);

/** \brief Assemble element matrix into global one for the Discontinuous Galerkin approximation
 *  @param A Global matrix stored as an SkSMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void DGElementAssembly(SkSMatrix<T_>& A);

/** \brief Assemble element matrix into global one for the Discontinuous Galerkin approximation
 *  @param [in] A Global matrix stored as an SkMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void DGElementAssembly(SkMatrix<T_>& A);

/** \brief Assemble element matrix into global one for the Discontinuous Galerkin approximation
 *  @param [in] A Global matrix stored as an SpMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void DGElementAssembly(SpMatrix<T_>& A);

/** \brief Assemble element matrix into global one for the Discontinuous Galerkin approximation
 *  @param [in] A Global matrix stored as an TrMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void DGElementAssembly(TrMatrix<T_>& A);

/** \brief Assemble side (edge or face) matrix into global one
 *  @param A Pointer to global matrix (abstract class: can be any of classes SkSMatrix, SkMatrix, SpMatrix)
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(Matrix<T_>* A);

/** \brief Assemble side (edge or face) matrix into global one
 *  @param [in] A Global matrix stored as an SkSMatrix instance
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(SkSMatrix<T_>& A);

/** \brief Assemble side (edge or face) matrix into global one
 *  @param [in] A Global matrix stored as an SkMatrix instance
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(SkMatrix<T_>& A);

/** \brief Assemble side (edge or face) matrix into global one
 *  @param [in] A Global matrix stored as an SpMatrix instance
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(SpMatrix<T_>& A);

/** \brief Assemble element vector into global one
 *  @param [in] v Global vector (Vect instance)
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(Vect<T_>& v);

/** \brief Assemble side (edge or face) vector into global one
 *  @param [in] v Global vector (Vect instance)
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(Vect<T_>& v);

/** \brief Assemble product of element matrix by element vector into global vector
 *  @param [in] el Reference to Element instance
 *  @param [in] x Global vector to multiply by (Vect instance)
 *  @param [out] b Global vector to add (Vect instance)
 */
    void AxbAssembly(const Element&  el,
                     const Vect<T_>& x,
                     Vect<T_>&       b);

/** \brief Assemble product of side matrix by side vector into global vector
 *  @param [in] sd Reference to Side instance
 *  @param [in] x Global vector to multiply by (Vect instance)
 *  @param [out] b Global vector (Vect instance)
 */
    void AxbAssembly(const Side&     sd,
                     const Vect<T_>& x,
                     Vect<T_>&       b);

/// \brief Return number of element nodes.
    size_t getNbNodes() const;

/// \brief Return number of element equations
    size_t getNbEq() const;

/** \brief Define a material property by an algebraic expression.
 *  @param [in] exp Algebraic expression
 *  @param [in] prop Property name
 *  @return Return value in expression evaluation:
 *  <ul>
 *    <li><tt>=0</tt>, Normal evaluation
 *    <li><tt>!=0</tt>, An error message is displayed
 *  </ul>
 */
    real_t setMaterialProperty(const string& exp,
                               const string& prop);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief LocalMatrix instance containing local matrix associated to current element
    LocalMatrix<T_,NEE_,NEE_> eMat, eA0, eA1, eA2;

/// \brief LocalMatrix instance containing local matrix associated to current side
    LocalMatrix<T_,NSE_,NSE_> sMat, sA0, sA1;

/// \brief LocalVect instance containing local right-hand side vector associated to current element
    LocalVect<T_,NEE_> eRHS;

/// \brief LocalVect instance containing local residual vector associated to current element
    LocalVect<T_,NEE_> eRes;

/// \brief LocalVect instance containing local right-hand side vector associated to current side
    LocalVect<T_,NSE_> sRHS;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    size_t                        _label;
    const Element*                _theElement;
    const Side*                   _theSide;
    ElementGeom                   _el_geo;
    LocalVect<Point<real_t>,NEN_> _x;
    LocalVect<T_,NEE_>            _eu, _ebf;
    LocalVect<T_,NSE_>            _su;
    T_                            _ssf[MAX_NBDOF_NODE];
    void updateBC(const Element&  el,
                  Vect<T_>*       bc=nullptr);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
