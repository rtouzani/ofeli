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

                            Definition of class DC2DT3
   for 2-D diffusion-convection scalar equation using 3-node triangular element

  ==============================================================================*/


#ifndef __DC2DT3_H
#define __DC2DT3_H

#include "equations/therm/Equa_Therm.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file DC2DT3.h
 *  \brief Definition file for class DC2DT3.
 */

/*! \class DC2DT3
 *  \ingroup Therm
 *  \brief Builds finite element arrays for thermal diffusion and convection in 2-D domains 
 *  using 3-Node triangles.
 *
 *  \details Note that members calculating element arrays have as an argument a real
 *  <tt>coef</tt> that will be multiplied by the contribution of the current element. 
 *  This makes possible testing different algorithms.
 *
 */

class DC2DT3 : public Equa_Therm<real_t,3,3,2,2>
{

 public:

    using Equa_Therm<real_t,3,3,2,2>::build;
    using Equa_Therm<real_t,3,3,2,2>::run;

/// \brief Default Constructor.
/// Constructs an empty equation.
    DC2DT3();

/// \brief Constructor using Mesh data
/// @param [in] ms Mesh instance
    DC2DT3(Mesh& ms);

/** \brief Constructor using Mesh and initial condition
 *  @param [in] ms Mesh instance
 *  @param [in] u Vect instance containing initial solution
 */
    DC2DT3(Mesh& ms, Vect<real_t>& u);

/// \brief Constructor for an element
/// @param el Pointer to Element instance
    DC2DT3(const Element* el);

/// \brief Constructor for a boundary side
/// @param [in] sd Pointer to Side instance
    DC2DT3(const Side* sd);

/** \brief Constructor for an element (transient case).
 *  @param [in] el Pointer to element
 *  @param [in] u Vect instance that contains solution at previous time step
 *  @param [in] time Current time value [Default: <tt>0</tt>]
 */
    DC2DT3(const Element*      el,
           const Vect<real_t>& u,
                 real_t        time=0.);

/** \brief Constructor for an element (transient case) with specification of time integration scheme.
 *  @param [in] el Pointer to element.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value [Default: <tt>0</tt>].
 *  @param [in] deltat Value of time step.
 *  @param [in] scheme Time Integration Scheme:
 *  <ul>
 *    <li><tt>FORWARD_EULER</tt> for Forward Euler scheme
 *    <li><tt>BACKWARD_EULER</tt> for Backward Euler scheme
 *    <li><tt>CRANK_NICOLSON</tt> for Crank-Nicolson Euler scheme
 *  </ul>
 */
    DC2DT3(const Element*      el,
           const Vect<real_t>& u,
                 real_t        time,
                 real_t        deltat,
                 int           scheme);

/** \brief Constructor for a boundary side (transient case).
 *  @param [in] sd Pointer to side.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value [Default: <tt>0</tt>]
 */
    DC2DT3(const Side*         sd,
           const Vect<real_t>& u,
                 real_t        time=0.);

/** \brief Constructor for a side (transient case) with specification of time integration scheme.
 *  @param [in] sd Pointer to side.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value [Default: <tt>0</tt>].
 *  @param [in] deltat Value of time step.
 *  @param [in] scheme Time Integration Scheme:
 *  <ul>
 *    <li><tt>FORWARD_EULER</tt> for Forward Euler scheme
 *    <li><tt>BACKWARD_EULER</tt> for Backward Euler scheme
 *    <li><tt>CRANK_NICOLSON</tt> for Crank-Nicolson Euler scheme
 *  </ul>
 */
    DC2DT3(const Side*         sd,
           const Vect<real_t>& u,
                 real_t        time,
                 real_t        deltat,
                 int           scheme);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Constructors for a mesh
    DC2DT3(Mesh&         mesh,
           Vect<real_t>& b,
           real_t&       t,
           real_t&       ts)
         : Equation<real_t,3,3,2,2>(mesh,b,t,ts) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Destructor
    ~DC2DT3();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
//    void withConvection(int fl) { fl=0; _with_convection=true; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Add lumped capacity matrix to left-hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void LCapacityToLHS(real_t coef=1);

/** \brief Add lumped capacity contribution to right-hand side after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void LCapacityToRHS(real_t coef=1);

/** \brief Add lumped capacity contribution to left and right-hand sides after multiplying 
 *  it by coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void LCapacity(real_t coef) { LCapacityToLHS(coef); LCapacityToRHS(coef); }

/** \brief Add Consistent capacity matrix to left-hand side after multiplying it 
 *  by coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void CapacityToLHS(real_t coef=1);

/** \brief Add Consistent capacity contribution to right-hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void CapacityToRHS(real_t coef=1);

/** \brief Add Consistent capacity contribution to left and right-hand sides after 
 *  multiplying it by coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Capacity(real_t coef=1) { CapacityToLHS(coef); CapacityToRHS(coef); }

/** \brief Add diffusion matrix to left-hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void Diffusion(real_t coef=1);

/** \brief Add diffusion matrix to left-hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  \details Case where the diffusivity matrix is given as an argument.
 *  @param [in] diff Diffusion matrix (class LocalMatrix).
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void Diffusion(const LocalMatrix<real_t,2,2>& diff,
                         real_t                   coef=1);

/** \brief Add diffusion contribution to right-hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  To be used for explicit diffusion term
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void DiffusionToRHS(real_t coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] v Constant velocity vector
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void Convection(const Point<real_t>& v,
                          real_t         coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient
 *  <tt>coef</tt>
 *  \details Case where velocity field is given by a vector <tt>v</tt>
 *  @param [in] v Velocity vector
 *  @param [in] coef Coefficient to multiply by added term (Default: <tt>1</tt>]
 */
    void Convection(const Vect<real_t>& v,
                          real_t        coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient
 *  <tt>coef</tt>
 *  \details Case where velocity field has been previouly defined
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void Convection(real_t coef=1);

/** \brief Add convection contribution to right-hand side after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  \details To be used for explicit convection term.
 *  @param [in] v Velocity vector
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void ConvectionToRHS(const Point<real_t>& v,
                               real_t         coef=1);

/** \brief Add convection contribution to right-hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  \details Case where velocity field has been previouly defined 
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void ConvectionToRHS(real_t coef=1);

/** \brief Add an edge linear exchange term to left and right-hand sides
 *  @param [in] coef Coefficient of exchange
 *  @param [in] T External value for exchange
 *  @remark This assumes a constant value of <tt>T</tt>
 */
    void LinearExchange(real_t coef,
                        real_t T);

/** \brief Add body right-hand side term to right hand side after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  @param [in] ud Instance of UserData or of a derived class. Contains
 *  a member function that provides body source.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void BodyRHS(UserData<real_t>& ud,
                 real_t            coef=1);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] bf Vector containing source at element nodes.
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size <tt>3</tt> or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BodyRHS(const Vect<real_t>& bf,
                       int           opt=GLOBAL_ARRAY);

/** \brief Add body right-hand side term to right hand side.
 *  \details Case where the body right-hand side is piecewise constant.
 *  @param [in] bf Value of thermal source (Constant in element).
 */
    void BodyRHS(real_t bf);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] ud Instance of UserData or of an inherited class. Contains
 *  a member function that provides body source.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void BoundaryRHS(UserData<real_t>& ud,
                     real_t            coef=1);

/// \brief Add boundary right-hand side flux to right hand side.
/// @param [in] flux Vector containing source at side nodes.
    void BoundaryRHS(real_t flux);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] b Vector containing source at side nodes
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size 2 or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BoundaryRHS(const Vect<real_t>& b,
                           int           opt=GLOBAL_ARRAY);

/** \brief Add contribution of periodic boundary condition (by a penalty technique).
 *  \details Boundary nodes where periodic boundary conditions are to be imposed must
 *  have codes equal to <tt>PERIODIC_A</tt> on one side and <tt>PERIODIC_B</tt> on the
 *  opposite side.
 *  @param [in] coef Value of penalty parameter [Default: <tt>1.e20</tt>]
 */
    void Periodic(real_t coef=1.e20);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// Return element contribution to energy.
    real_t Energy(const Vect<real_t>&     u,
                        UserData<real_t>& ud);
    real_t Energy(const Vect<real_t>& u);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// Return energy gradient in element.
    void EnerGrad(const Vect<real_t>& u,
                  UserData<real_t>&   ud);
    void EnerGrad(const Vect<real_t>& u,
                  const Vect<real_t>& f);
    void EnerGrad(const Vect<real_t>& u);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return (constant) heat flux in element.
    Point<real_t>& Flux() const;

/// \brief Return gradient of a vector in element.
/// @param [in] u Vector for which gradient is computed.
    Point<real_t>& Grad(const LocalVect<real_t,3>& u) const;

/** \brief Return gradient of a vector in element.
 *  @param [in] u Global vector for which gradient is computed.
 *  Vector <tt>u</tt> has as size the total number of nodes
 */
    Point<real_t>& Grad(const Vect<real_t>& u) const;

/** \brief Set equation input data
 *  @param [in] opt Parameter to select type of input (enumerated values)
 *  <ul>
 *     <li><tt>INITIAL_FIELD</tt>: Initial temperature
 *     <li><tt>BOUNDARY_CONDITION_DATA</tt>: Boundary condition (Dirichlet)
 *     <li><tt>SOURCE_DATA</tt>: Heat source
 *     <li><tt>FLUX_DATA</tt>: Heat flux (Neumann boundary condition)
 *     <li><tt>VELOCITY_FIELD</tt>: Velocity vector (for the convection term)
 *  </ul>
 *  @param [in] u Vector containing input data
 */
    void setInput(EqDataType    opt,
                  Vect<real_t>& u);

/** \brief Set Joule heating term as source
 *  @param [in] sigma Vect instance containing electric conductivity (elementwise)
 *  @param [in] psi Vect instance containing electric potential (elementwise)
 */
    void JouleHeating(const Vect<real_t>& sigma,
                      const Vect<real_t>& psi);

/** \brief Build the linear system of equations
 *  \details Before using this function, one must have properly selected 
 *  appropriate options for:
 *  <ul>
 *     <li>The choice of a steady state or transient analysis. By default, the analysis is stationary
 *     <li>In the case of transient analysis, the choice of a time integration scheme
 *         and a lumped or consistent capacity matrix. If transient analysis is chosen, the lumped
 *         capacity matrix option is chosen by default, and the implicit Euler scheme is used
 *         by default for time integration.
 *  </ul>
 */
    void build() { Equa_Therm<real_t,3,3,2,2>::build(); }
    void build(TimeStepping& s) { Equa_Therm<real_t,3,3,2,2>::build(s); }

/** \brief Run the equation
 *  \details If the analysis (see function setAnalysis) is <tt>STEADY_STATE</tt>, then
 *  the function solves the stationary equation.\n
 *  If the analysis is <tt>TRANSIENT</tt>, then the function performs time stepping
 *  until the final time is reached.
 */
//    int run() { return Equa_Therm<real_t,3,3,2,2>::run(); }


 protected:
    void set(const Element* el);
    void set(const Side* sd);

 private:
    mutable Point<real_t> _f;
    real_t _h;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
