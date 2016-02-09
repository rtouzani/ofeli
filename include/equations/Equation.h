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

                 Definition of abstract class 'Equations'
                       for finite element equations

  ==============================================================================*/


#ifndef __EQUATION_H
#define __EQUATION_H

#include "equations/AbsEqua.h"

extern FunctionParser theParser;


namespace OFELI {

/*! \file Equation.h
 *  \brief Definition file for class Equation.
 */


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_> class Equation;
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
 */

template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
class Equation : virtual public AbsEqua<T_>
{

 public:

  using AbsEqua<T_>::_theMesh;
  using AbsEqua<T_>::_field_type;
  using AbsEqua<T_>::_matrix_type;
  using AbsEqua<T_>::_time_scheme;
  using AbsEqua<T_>::_verbose;
  using AbsEqua<T_>::_step;
  using AbsEqua<T_>::_nb_fields;
  using AbsEqua<T_>::_theta;
  using AbsEqua<T_>::_alpha;
  using AbsEqua<T_>::_beta;
  using AbsEqua<T_>::_time_step;
  using AbsEqua<T_>::_time;
  using AbsEqua<T_>::_A;
  using AbsEqua<T_>::_b;
  using AbsEqua<T_>::_bc;
  using AbsEqua<T_>::_bf;
  using AbsEqua<T_>::_sf;
  using AbsEqua<T_>::_u;
  using AbsEqua<T_>::_sol_given;
  using AbsEqua<T_>::_init_given;
  using AbsEqua<T_>::_bc_given;
  using AbsEqua<T_>::_bf_given;
  using AbsEqua<T_>::_sf_given;

/// Default constructor.
/// Constructs an "empty" equation
    Equation() { _b = NULL; }

/// \brief Constructor with mesh instance.
/// @param [in] mesh Mesh instance
    Equation(Mesh &mesh);

/** \brief Constructor with mesh instance, matrix and right-hand side.
 *  @param [in] mesh Mesh instance
 *  @param [in] b Vect instance containing Right-hand side.
 *  @param [in] t Time value
 *  @param [in] ts Time step
 */
    Equation(Mesh&     mesh,
             Vect<T_>& b,
             real_t&   t,
             real_t&   ts);

/// \brief Constructor using Element data.
/// @param [in] el Pointer to Element
    Equation(const Element* el);

/// Constructor using Side data.
    Equation(const Side* sd);

/** \brief Constructor using element data, solution at previous time step and time value.
 *  @param [in] el Pointer to element
 *  @param [in] u Vect instance containing solution at previous time step
 *  @param [in] time Time value (Default value is 0)
 */
    Equation(const Element*  el,
             const Vect<T_>& u,
                   real_t    time=0);

/** \brief Constructor using side data, solution at previous time step and time value.
 *  @param [in] sd Pointer to side
 *  @param [in] u Vect instance containing solution at previous time step
 *  @param [in] time Time value (Default value is 0)
 */
    Equation(const Side*     sd,
             const Vect<T_>& u,
                   real_t    time=0);

/// \brief Destructor
    virtual ~Equation() { }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Initialize equation.
 *  @param [in] mesh Mesh instance
 *  @param [in] ts Time step value
 */
    virtual void initEquation(Mesh&  mesh,
                              real_t ts);

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Return problem matrix: Case of a SpMatrix (sparse).
/// @param [in] A 
    void getMatrix(const SpMatrix<T_>& A) const { eMat.Localize(_theElement,A); }

/// \brief Return problem matrix: Case of a SkMatrix (skyline)
    void getMatrix(const SkMatrix<T_>& A) const { eMat.Localize(_theElement,A); }

/// \brief Return problem matrix: Case of a SkSMatrix (symmetric skyline)
    void getMatrix(const SkSMatrix<T_>& A) const { eMat.Localize(_theElement,A); }

/// \brief Return problem matrix (Case of a SpMatrix (sparse))
    void getSolution(const Vect<T_>& u) const { ePrev.Localize(_theElement,u); }

/// \brief Return problem right-hand side vector
    void getRHS(const Vect<T_>& b) const { eRHS.Localize(_theElement,b); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Update Right-Hand side by taking into account essential boundary conditions
 *  @param [in] el Reference to current element instance
 *  @param [in] bc Vector that contains imposed values at all DOFs
 */
    void updateBC(const Element&  el,
                  const Vect<T_>& bc);

/** \brief Update Right-Hand side by taking into account essential boundary conditions
 *  @param [in] bc Vector that contains imposed values at all DOFs
 *  @remark The current element is pointed by <tt>_theElement</tt>
 */
    void updateBC(const Vect<T_>& bc) { updateBC(TheElement,bc); }

/** \brief Update element matrix to impose bc by diagonalization technique
 *  @param [in] dof_type DOF type option. To choose among the enumerated values:
 *  <ul>
 *    <li><tt>NODE_FIELD</tt>, DOFs are supported by nodes [Default]
 *    <li><tt>ELEMENT_FIELD</tt>, DOFs are supported by elements
 *    <li><tt>SIDE_FIELD</tt>, DOFs are supported by sides
 *  </ul>
 *  @param [in] dof DOF setting:
 *  <ul>
 *     <li><tt> = 0</tt>, All DOFs are taken into account [Default]
 *     <li><tt>!= 0</tt>, Only DOF No. <tt>dof</tt> is handled in the system
 *  </ul>
 */
    void DiagBC(int dof_type=NODE_DOF,
                int dof=0);

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
    void ElementNodeVector(const Vect<T_>&           b,
                                 LocalVect<T_,NEE_>& be);

/** \brief Localize Element Vector from a Vect instance.
 *  @param [in] b Global vector to be localized.
 *  @param [out] be Local vector, the length of which is the total number of element equations. 
 *  @remark Vector <tt>b</tt> is assumed to contain only one degree of freedom by node.
 */
    void ElementNodeVectorSingleDOF(const Vect<T_>&           b,
                                          LocalVect<T_,NEN_>& be);

/** \brief Localize Element Vector from a Vect instance.
 *  @param [in] b Global vector to be localized.
 *  @param [out] be Local vector, the length of which is the total number of element equations. 
 *  @param [in] dof Degree of freedom to transfer to the local vector
 *  @remark Only yhe dega dof is transferred to the local vector
 */
    void ElementNodeVector(const Vect<T_>&           b,
                                 LocalVect<T_,NEN_>& be,
                                 int                 dof);

/** \brief Localize Element Vector from a Vect instance.
 *  @param [in] b Global vector to be localized.
 *  @param [out] be Local vector, the length of which is 
 */
    void ElementSideVector(const Vect<T_>&           b,
                                 LocalVect<T_,NSE_>& be);

/** \brief Localize Element Vector.
 *  @param [in] b Global vector to be localized
 *  @param [in] dof_type DOF type option. To choose among the enumerated values:
 *  <ul>
 *    <li><tt>NODE_FIELD</tt>, DOFs are supported by nodes [Default]
 *    <li><tt>ELEMENT_FIELD</tt>, DOFs are supported by elements
 *    <li><tt>SIDE_FIELD</tt>, DOFs are supported by sides
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
                             int       dof_type=NODE_FIELD,
                             int       flag=0);

/** \brief Localize Side Vector.
 *  @param [in] b Global vector to be localized
 *  <ul>
 *     <li><tt>NODE_FIELD</tt>, DOFs are supported by nodes [ default ]
 *     <li><tt>ELEMENT_FIELD</tt>, DOFs are supported by elements
 *     <li><tt>SIDE_FIELD</tt>, DOFs are supported by sides
 *  </ul>
 *  The resulting local vector can be accessed by attribute <tt>ePrev</tt>.
 *  @remark This member function is to be used if a constructor with Side was invoked.
 *  It uses the Side pointer <tt>_theSide</tt>
 */
    void SideVector(const Vect<T_>& b);

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
    void ElementAssembly(Matrix<T_>* A) { A->Assembly(TheElement,eMat.get()); }

#if defined (USE_PETSC)
/** \brief Assemble element matrix into global one
 *  @param A Reference to global matrix
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(PETScMatrix<T_>& A) { ElementAssembly(TheElement,A); }

/** \brief Assemble side matrix into global one
 *  @param A Reference to global matrix
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(PETScMatrix<T_>& A) { SideAssembly(TheSide,A); }

/** \brief Assemble element right-hand side vector into global one
 *  @param b Reference to global right-hand side vector
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(PETScVect<T_>& b) { ElementAssembly(TheElement,b); }

/** \brief Assemble side right-hand side vector into global one
 *  @param b Reference to global right-hand side vector
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(PETScVect<T_>& b) { SideAssembly(TheSide,b); }
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void ElementAssembly(Element& el, Matrix<T_>* A) { A->Assembly(el,eMat.get()); }
    void ElementAssembly(Element& el, SkMatrix<T_>& A) { A.Assembly(el,eMat.get()); }
    void ElementAssembly(Element& el, SkSMatrix<T_>& A) { A.Assembly(el,eMat.get()); }
    void ElementAssembly(Element& el, SpMatrix<T_>& A) { A.Assembly(el,eMat.get()); }
    void ElementAssembly(Element& el, BMatrix<T_>& A) { A.Assembly(el,eMat.get()); }
    void ElementAssembly(Element& el, TrMatrix<T_>& A) { A.Assembly(el,eMat.get()); }
    void SideAssembly(Side& sd, Matrix<T_>* A) { A->Assembly(sd,sMat.get()); }
    void SideAssembly(Side& sd, BMatrix<T_>& A) { A.Assembly(sd,sMat.get()); }
    void SideAssembly(Side& sd, SkMatrix<T_>& A) { A.Assembly(sd,sMat.get()); }
    void SideAssembly(Side& sd, SkSMatrix<T_>& A) { A.Assembly(sd,sMat.get()); }
    void SideAssembly(Side& sd, SpMatrix<T_>& A) { A.Assembly(sd,sMat.get()); }
    void ElementAssembly(Element& el, Vect<T_>& v) { v.Assembly(el,eRHS.get()); }
    void SideAssembly(Side& sd, Vect<T_>& v) { v.Assembly(sd,sRHS.get()); }
#if defined(USE_PETSC)
    void updateBC(const Element& el, const PETScVect<T_>& bc);
    void updateBC(const PETScVect<T_>& bc) { updateBC(TheElement,bc); }
    void ElementAssembly(Element& el, PETScMatrix<T_>& A) { A.Assembly(el,eMat.get()); }
    void SideAssembly(Side& sd, PETScMatrix<T_>& A) { A.Assembly(sd,sMat.get()); }
    void ElementAssembly(Element& el, PETScVect<T_>& v) { v.Assembly(el,eRHS.get()); }
    void SideAssembly(Side& sd, PETScVect<T_>& v) { v.Assembly(sd,sRHS.get()); }
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Assemble element matrix into global one
 *  @param A Global matrix stored as a BMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(BMatrix<T_>& A) { ElementAssembly(TheElement,A); }

/** \brief Assemble element matrix into global one
 *  @param A Global matrix stored as an SkSMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(SkSMatrix<T_>& A) { ElementAssembly(TheElement,A); }

/** \brief Assemble element matrix into global one
 *  @param [in] A Global matrix stored as an SkMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(SkMatrix<T_>& A) { ElementAssembly(TheElement,A); }

/** \brief Assemble element matrix into global one
 *  @param [in] A Global matrix stored as an SpMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(SpMatrix<T_>& A) { ElementAssembly(TheElement,A); }

/** \brief Assemble element matrix into global one
 *  @param [in] A Global matrix stored as an TrMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(TrMatrix<T_>& A) { ElementAssembly(TheElement,A); }

/** \brief Assemble element matrix into global one for the Discontinuous Galerkin approximation
 *  @param A Pointer to global matrix (abstract class: can be any of classes SkSMatrix, SkMatrix, SpMatrix)
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void DGElementAssembly(Matrix<T_>* A) { A->DGAssembly(TheElement,eMat.get()); }

/** \brief Assemble element matrix into global one for the Discontinuous Galerkin approximation
 *  @param A Global matrix stored as an SkSMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void DGElementAssembly(SkSMatrix<T_>& A) { A.DGAssembly(TheElement,eMat.get()); }

/** \brief Assemble element matrix into global one for the Discontinuous Galerkin approximation
 *  @param [in] A Global matrix stored as an SkMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void DGElementAssembly(SkMatrix<T_>& A) { A.DGAssembly(TheElement,eMat.get()); }

/** \brief Assemble element matrix into global one for the Discontinuous Galerkin approximation
 *  @param [in] A Global matrix stored as an SpMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void DGElementAssembly(SpMatrix<T_>& A) { A.DGAssembly(TheElement,eMat.get()); }

/** \brief Assemble element matrix into global one for the Discontinuous Galerkin approximation
 *  @param [in] A Global matrix stored as an TrMatrix instance
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void DGElementAssembly(TrMatrix<T_>& A) { A.DGAssembly(TheElement,eMat.get()); }

/** \brief Assemble side (edge or face) matrix into global one
 *  @param A Pointer to global matrix (abstract class: can be any of classes SkSMatrix, SkMatrix, SpMatrix)
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(Matrix<T_>* A) { SideAssembly(TheSide,A); }

/** \brief Assemble side (edge or face) matrix into global one
 *  @param [in] A Global matrix stored as an SkSMatrix instance
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(SkSMatrix<T_>& A) { SideAssembly(TheSide,A); }

/** \brief Assemble side (edge or face) matrix into global one
 *  @param [in] A Global matrix stored as an SkMatrix instance
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(SkMatrix<T_>& A) { SideAssembly(TheSide,A); }

/** \brief Assemble side (edge or face) matrix into global one
 *  @param [in] A Global matrix stored as an SpMatrix instance
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(SpMatrix<T_>& A) { SideAssembly(TheSide,A); }

/** \brief Assemble element vector into global one
 *  @param [in] v Global vector (Vect instance)
 *  @warning The element pointer is given by the global variable <tt>theElement</tt>
 */
    void ElementAssembly(Vect<T_>& v) { ElementAssembly(TheElement,v); }

/** \brief Assemble side (edge or face) vector into global one
 *  @param [in] v Global vector (Vect instance)
 *  @warning The side pointer is given by the global variable <tt>theSide</tt>
 */
    void SideAssembly(Vect<T_>& v) { SideAssembly(TheSide,v); }

/** \brief Assemble product of element matrix by element vector into global vector
 *  @param [in] el Reference to Element instance
 *  @param [in] x Global vector to multiply by (Vect instance)
 *  @param [out] b Global vector to add (Vect instance)
 */
    void AxbAssembly(const Element&  el,
                     const Vect<T_>& x,
                           Vect<T_>& b);

/** \brief Assemble product of side matrix by side vector into global vector
 *  @param [in] sd Reference to Side instance
 *  @param [in] x Global vector to multiply by (Vect instance)
 *  @param [out] b Global vector (Vect instance)
 */
    void AxbAssembly(const Side&     sd,
                     const Vect<T_>& x,
                           Vect<T_>& b);

/// \brief Return number of element nodes.
    size_t getNbNodes() const;

/// \brief Return number of element equations
    size_t getNbEq() const;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Set problem data
    void setBoundaryCondition(Vect<T_>& bc);

/// \brief Set Body force vector
    void setBodyForce(Vect<T_>& bf);

/// \brief Set surface force vector
    void setSurfaceForce(Vect<T_>& sf);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Set initial solution (previous time step)
    void setInitialSolution(const Vect<T_>& u);

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

/// \brief LocalMatrix instance containing local matrix associated to current element
    LocalMatrix<T_,NEE_,NEE_> eMat, eA0, eA1, eA2;

/// \brief LocalMatrix instance containing local matrix associated to current side
    LocalMatrix<T_,NSE_,NSE_> sMat;

/// \brief LocalVect instance containing local vector associated to current element
/// \details This vector has been stored as the one at previous iteration or previous
/// time step
    LocalVect<T_,NEE_> ePrev;

/// \brief LocalVect instance containing local right-hand side vector associated to current element
    LocalVect<T_,NEE_> eRHS;

/// \brief LocalVect instance containing local residual vector associated to current element
    LocalVect<T_,NEE_> eRes;

/// \brief LocalVect instance containing local right-hand side vector associated to current side
    LocalVect<T_,NSE_> sRHS;

 protected:

/// Set element arrays to zero
    void Init(const Element* el);

/// Set side arrays to zero
    void Init(const Side* sd);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    size_t                        _nb_dof, _label;
    const Element*                _theElement;
    const Side*                   _theSide;
    real_t                        _volume, _area, _length, _det;
    Point<real_t>                 _center;
    LocalVect<Point<real_t>,NEN_> _dSh;
    Vect<T_>*                     _BodyForce;
    LocalVect<Point<real_t>,NEN_> _x;
    T_                            _body_source, _bound_source;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};


///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
Equation<T_,NEN_,NEE_,NSN_,NSE_>::Equation(Mesh& mesh)
{
   _sol_given = _bc_given = _init_given = _bf_given = _sf_given = false; 
   _time_step = 0.1;
   initEquation(mesh,_time_step);
   _time = 0.;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
Equation<T_,NEN_,NEE_,NSN_,NSE_>::Equation(Mesh&     mesh,
                                           Vect<T_>& b,
                                           real_t&   t,
                                           real_t&   ts)
{
   _sol_given = _bc_given = _init_given = _bf_given = _sf_given = false; 
   initEquation(mesh,ts);
   _time = t;
   this->setInput(INITIAL_FIELD,b);
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
Equation<T_,NEN_,NEE_,NSN_,NSE_>::Equation(const Element* el)
{
   _label = el->n();
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
Equation<T_,NEN_,NEE_,NSN_,NSE_>::Equation(const Side* sd)
{
   _label = sd->n();
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
Equation<T_,NEN_,NEE_,NSN_,NSE_>::Equation(const Element*  el,
                                           const Vect<T_>& u,
                                                 real_t    time)
{
   _label = el->n();
   _time = time;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
Equation<T_,NEN_,NEE_,NSN_,NSE_>::Equation(const Side*     sd,
                                           const Vect<T_>& u,
                                                 real_t    time)
{
   _label = sd->n();
   _time = time;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::initEquation(Mesh&  mesh,
                                                    real_t ts)
{
   _theMesh = &mesh;
   _theMesh->removeImposedDOF();
   _time_step = ts;
   _b = new Vect<T_>(_theMesh->getNbEq());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::updateBC(const Element&  el,
                                                const Vect<T_>& bc)
{
   _nb_dof = NEE_/NEN_;
   size_t in=1;
   for (size_t i=1; i<=NEN_; ++i) {
      for (size_t k=1; k<=_nb_dof; ++k, ++in) {
         size_t jn=1;
         for (size_t j=1; j<=NEN_; ++j) {
            for (size_t l=1; l<=_nb_dof; ++l, ++jn) {
               if (el(j)->getCode(l) > 0)
                  eRHS(in) -= eMat(in,jn) * bc(el(j)->getFirstDOF()+l-1);
            }
         }
      }
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::DiagBC(int dof_type,
                                              int dof)
{
   size_t in, jn;
   _nb_dof = NEE_/NEN_;
   if (dof_type==NODE_FIELD) {
      if (dof) {
         for (size_t i=1; i<=NEN_; i++) {
            Node *nd1 = (*_theElement)(i);
            if (nd1->getCode(dof)) {
               for (size_t j=1; j<=NEN_; j++) {
                  Node *nd2 = (*_theElement)(j);
                  in = nd1->n(), jn = nd2->n();
                  if (in != jn)
                     eMat(in,jn) = T_(0.);
               }
            }
         }
      }
      else {
         size_t in=1;
         for (size_t i=1; i<=NEN_; i++) {
            for (size_t k=1; k<=_nb_dof; k++) {
               size_t jn = 1;
               for (size_t j=1; j<=NEN_; j++) {
                  for (size_t l=1; l<=_nb_dof; l++) {
                     if (in != jn && (*_theElement)(i)->getCode(k))
                        eMat(in,jn) = T_(0.);
                     jn++;
                  }
               }
               in++;
            }
         }
      }
   }

   else if (dof_type==SIDE_FIELD) {
      if (dof) {
         for (size_t i=1; i<=NEN_; i++) {
            Side *sd1 = _theElement->getPtrSide(i);
            if (sd1->getCode(dof))
               for (size_t j=1; j<=NEN_; j++)
                  if (i != j)
                     eMat(i,j) = T_(0.);
         }
      }
      else {
         for (size_t i=1; i<=NEN_; i++) {
            for (size_t k=1; k<=_nb_dof; k++) {
               if (_theElement->getPtrSide(i)->getCode(k)) {
                  size_t in = (*_theElement)(i)->getDOF(k);
                  for (size_t j=1; j<=NEN_; j++) {
                     for (size_t l=1; l<=_nb_dof; l++) {
                        size_t jn = _theElement->getPtrSide(j)->getDOF(l);
                        if (in != jn)
                           eMat(in,jn) = T_(0.);
                     }
                  }
               }
            }
         }
      }
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::setBodyForce(const Vect<T_>& f)
{
   _BodyForce = &f;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::setResidue()
{
   eMat.Mult(ePrev,eRes);
   eRes -= eRHS;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::LocalNodeVector(Vect<T_>& b)
{
   size_t k = 0;
   Node *nd;
   for (size_t n=1; n<=NEN_; ++n) {
      nd = (*_theElement)(n);
      for (size_t j=1; j<=nd->getNbDOF(); j++)
         ePrev[k++] = b(nd->n(),j);
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementNodeVector(const Vect<T_>&           b,
                                                               LocalVect<T_,NEE_>& be)
{
   size_t k = 0;
   Node *nd;
   for (size_t n=1; n<=NEN_; ++n) {
      nd = (*_theElement)(n);
      for (size_t j=1; j<=nd->getNbDOF(); j++)
         be[k++] = b(nd->getDOF(j));
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementNodeVector(const Vect<T_>&           b,
                                                               LocalVect<T_,NEN_>& be,
                                                               int                 dof)
{
   size_t k = 0;
   Node *nd;
   for (size_t n=1; n<=NEN_; ++n) {
      nd = (*_theElement)(n);
      for (size_t j=1; j<=nd->getNbDOF(); j++)
         if (j==dof)
            be[k++] = b(nd->getDOF(j));
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementNodeVectorSingleDOF(const Vect<T_>&           b,
                                                                        LocalVect<T_,NEN_>& be)
{
   for (size_t n=1; n<=NEN_; ++n)
      be(n) = b((*_theElement)(n)->n());
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementSideVector(const Vect<T_>&           b,
                                                               LocalVect<T_,NSE_>& be)
{
   size_t k = 0;
   Side *sd;
   for (size_t n=1; n<=NSE_; ++n) {
      sd = _theElement->getPtrSide(n);
      for (size_t j=1; j<=sd->getNbDOF(); j++)
         if (sd->getDOF(j))
            be[k++] = b(sd->getDOF(j));
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementVector(const Vect<T_>& b,
                                                           int       dof_type,
                                                           int       flag)
{
   size_t k=0;
   switch (dof_type) {

      case NODE_FIELD:
         if (flag==0) {
            for (size_t n=1; n<=NEN_; ++n) {
               Node *nd = (*_theElement)(n);
               for (size_t j=1; j<=nd->getNbDOF(); j++)
                  if (nd->getDOF(j))
                     ePrev[k++] = b(nd->getDOF(j));
            }
         }
         else {
            for (size_t n=1; n<=NEN_; ++n)
               ePrev[k++] = b(_theElement->getNodeLabel(n));
         }
         break;

      case SIDE_FIELD:
         if (flag==0) {
            for (size_t n=1; n<=NEN_; n++) {
               Side *sd=_theElement->getPtrSide(n);
               for (size_t j=1; j<=sd->getNbDOF(); j++)
                  if (sd->getDOF(j))
                     ePrev[k++] = b(sd->getDOF(j));
            }
         }
         else {
            for (size_t n=1; n<=NEN_; n++)
               ePrev[k++] = b(_theElement->getSideLabel(n));
         }
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideVector(const Vect<T_>& b)
{
   Node *nd;
   size_t k=0;
   for (size_t n=1; n<=NSN_; n++) {
      nd = (*_theSide)(n);
      for (size_t j=1; j<=nd->getNbDOF(); j++)
         ePrev[k++] = b(nd->getDOF(j));
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::ElementNodeCoordinates()
{
   for (size_t n=1; n<=NEN_; n++)
      _x[n-1] = (*_theElement)(n)->getCoord();
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::SideNodeCoordinates()
{
   for (size_t n=1; n<=NSN_; n++)
      _x[n-1] = (*_theSide)(n)->getCoord();
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
size_t Equation<T_,NEN_,NEE_,NSN_,NSE_>::getNbNodes() const
{
   if (_theElement)
      return NEN_;
   else
      return NSN_;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
size_t Equation<T_,NEN_,NEE_,NSN_,NSE_>::getNbEq() const
{
   if (_theElement)
      return NEE_;
   else
      return NSE_;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::AxbAssembly(const Element&  el,
                                                   const Vect<T_>& x,
                                                         Vect<T_>& b)
{
   size_t ii=0, jj=0, nb_dof=NEE_/NEN_, ik, jl;
   for (size_t i=1; i<=NEN_; ++i) {
      Node *ndi = el(i);
      for (size_t k=1; k<=nb_dof; ++k) {
         ii++, ik = ndi->getDOF(k);
         for (size_t j=1; j<=NEN_; ++j) {
            Node *ndj = el(j);
            for (size_t l=1; l<=nb_dof; ++l) {
               jj++, jl = ndj->getDOF(l);
               if (ik)
                  b(ik) += eMat(ii,jj)*x(jl);
            }
         }
      }
   }
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::AxbAssembly(const Side&     sd,
                                                   const Vect<T_>& x,
                                                         Vect<T_>& b)
{
   size_t ii=0, jj=0, nb_dof=NSE_/NSN_, ik, jl;
   for (size_t i=1; i<=NSN_; ++i) {
      Node *ndi = sd(i);
      for (size_t k=1; k<=nb_dof; ++k) {
         ii++, ik = ndi->getDOF(k);
         for (size_t j=1; j<=NSN_; ++j) {
            Node *ndj = sd(j);
            for (size_t l=1; l<=nb_dof; ++l) {
               jj++, jl = ndj->getDOF(l);
               if (ik)
                  b(ik) += sMat(ii,jj)*x(jl);
            }
         }
      }
   }
}


#if defined (USE_PETSC)
template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::updateBC(const Element&       el,
                                                const PETScVect<T_>& bc)
{
   size_t in = 0;
   for (size_t i=1; i<=NEN_; ++i) {
      for (size_t k=1; k<=_nb_dof; ++k) {
         in++;
         size_t nn = el(i)->getFirstDOF() + k - 1;
         if (el(i)->getDOF(k) == 0) {
            size_t jn = 0;
            for (size_t j=1; j<=NEN_; ++j) {
               for (size_t l=1; l<=_nb_dof; ++l) {
                  jn++;
                  if (el(j)->getDOF(l) > 0)
                     eRHS(jn) -= eMat(jn,in) * bc(nn);
               }
            }
         }
      }
   }
}
#endif


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
real_t Equation<T_,NEN_,NEE_,NSN_,NSE_>::setMaterialProperty(const string& exp,
                                                             const string& prop)
{
   int err;
   real_t d[4], r;
   theParser.Parse(exp,"x,y,z,t");
   d[0] = _center.x;
   d[1] = _center.y;
   d[2] = _center.z;
   d[3] = _time;
   r = theParser.Eval(d);
   try {
      if ((err=theParser.EvalError()))
         THROW_RT("setMaterialProperty(string,string): Error in evaluation in expression " + prop);
   }
   CATCH("Equation<>");
   return r;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::Init(const Element* el)
{
   _theElement = el;
   _theSide = NULL;
}


template<class T_, size_t NEN_, size_t NEE_, size_t NSN_, size_t NSE_>
void Equation<T_,NEN_,NEE_,NSN_,NSE_>::Init(const Side* sd)
{
   _theSide = sd;
   _theElement = NULL;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
