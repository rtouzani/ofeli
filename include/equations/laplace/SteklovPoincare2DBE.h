/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

                       Definition of class SteklovPoincare2DBE
              for boundary elements for the 2-D Steklov Poincare solver
                          using P0 collocation elements

  ==============================================================================*/


#ifndef __STEKLOV_POINCARE_2DBE_H
#define __STEKLOV_POINCARE_2DBE_H

#include "mesh/Mesh.h"
#include "linear_algebra/DMatrix.h"
#include "solvers/Prec.h"
#include "solvers/GMRes.h"
#include "linear_algebra/Vect.h"
#include "shape_functions/Triang3.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file SteklovPoincare2DBE.h
 *  \brief Definition file for class SteklovPoincare2DBE.
 */

/*! \class SteklovPoincare2DBE
 *  \ingroup Laplace
 *  \brief Solver of the Steklov Poincare problem in 2-D geometries using piecewie constant boundary elemen
 *
 *  \details SteklovPoincare2DBE solves the Steklov Poincare problem in 2-D: Given the trace of a harmonic function
 *  on the boundary of a given (inner or outer) domain, this class computes the normal derivative of the function.
 *  The normal is considered as oriented out of the bounded (inner) domain in both inner and outer configurations.
 *  The numerical approximation uses piecewise constant (<tt>P<sub>0</sub></tt>) approximation on edges of the boundary.
 *  Solution is obtained from the GMRES iterative solver without preconditioning.
 *  The given data is the vector (instance of class Vect) of piecewise constant values of the harmonic function
 *  on the boundary and the returned solution is piecewise constant value of the normal derivative considered either
 *  as a Vect instance.
 *
 *  @note Although the mesh of the inner domain is not necessary to solve the problem, this one must be provided in 
 *  order to calculate the outward normal.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class SteklovPoincare2DBE {

 public:

/// \brief Default Constructor
/// @param [in] ext Boolean variable to say if the domain is external (true) or internal (false: Default value).
    SteklovPoincare2DBE(bool ext=false)
   { 
      _ext = -1;
      if (ext)
         _ext = 1; 
   }

/** \brief Constructor using mesh data.
 *  \details This constructor calls member function setMesh.
 *  @param [in] mesh Reference to mesh instance.
 *  @param [in] ext Boolean variable to say if the domain is external (true) or internal (false: Default value).
 */
    SteklovPoincare2DBE(const Mesh& mesh,
                              bool  ext=false);

/** \brief Constructor that solves the Steklov Poincare problem.
 *  \details This constructor calls member function setMesh and Solve.
 *  @param [in] mesh Reference to mesh instance.
 *  @param [in] g Vect instance that contains imposed solution on the boundary
 *  @param [in] b Vect instance that contains the left hand side in input and the solution in output
 *  @param [in] ext Boolean variable to say if the domain is external (true) or internal (false: Default value).
 */
    SteklovPoincare2DBE(const Mesh&         mesh,
                        const Vect<real_t>& g,
                              Vect<real_t>& b,
                              bool          ext=false);

/// \brief Destructor
    ~SteklovPoincare2DBE() { }

/// \brief set Mesh instance
/// @param [in] mesh Mesh instance
/// @param [in] ext Boolean variable to say if the domain is external (true) or internal (false: Default value).
    void setMesh(const Mesh& mesh,
                       bool  ext=false);

/** \brief Build equation left and right-hand sides for <tt>P<sub>0</sub></tt> (piecewise
 *  constant) approximation.
 *  \details This member function is to be used if the constructor using <tt>mesh</tt>,
 *  <tt>b</tt> and <tt>g</tt> has been used.
 */
    void Solve();

/** \brief Build equation left and right-hand sides for <tt>P<sub>0</sub></tt>
 *  (piecewise constant) approximation
 *  \details This member function is to be used if the constructor using \c mesh has been used.
 *  It concerns cases where the imposed boundary condition is given by sides
 *  @param [in] g Vector that contains imposed solution on the boundary
 *  @param [in] b Vector that contains the left hand side in input and the solution in output
 */
    int Solve(      Vect<real_t>& b,
              const Vect<real_t>& g);

 private:

   int                        _ext;
   const Mesh                 *_theMesh;
   size_t                     _nb_eq;
   Vect<Point<real_t> >       _nn, _ttg, _center;
   Vect<real_t>               _length;
   real_t                     _h;
   SpMatrix<real_t>           _A;

// Return integral of the Green function over the side sd at point x 
   inline real_t single_layer(      size_t         j,
                               const Point<real_t>& z) const
   {
      real_t t = z.NNorm();
      if (t < OFELI_EPSMCH)
         return -0.25/OFELI_PI*_h*I1(0.5*_h);
      else
         return -0.125/OFELI_PI*_h*I2(0.25*_h*_h,z*_ttg(j),t);
   }

// Return integral of the normal derivative of the Green function over the side sd at point x
   inline real_t double_layer(      size_t         j,
                               const Point<real_t>& z) const
   {
      real_t t = _nn(j)*z;
      if (fabs(t) < OFELI_EPSMCH)
         return 0;
      else
         return -0.25/OFELI_PI*t*I3(0.25*_h*_h,z*_ttg(j),z.NNorm()); 
   }

// Return integral of log|a*t| over (-1,1)
   inline real_t I1(real_t a) const { return 2*log(a)-2; }

// Return integral of log(a*t^2+b*t+c) over (-1,1) when 4*a*c-b^2>0.
   inline real_t I2(real_t a,
                    real_t b,
                    real_t c) const
   {
      real_t d = sqrt(4*a*c-b*b);
      return 0.5*(log(a-b+c)*(2*a-b)+log(a+b+c)*(2*a+b)+2*(atan((2*a+b)/d)+atan((2*a-b)/d))*sqrt(4*a*c-b*b)-8*a)/a;
   }

// Return integral of 1/(a*t^2+b*t+c) over (-1,1) when 4*a*c-b^2>0.
   inline real_t I3(real_t a,
                    real_t b,
                    real_t c) const
   {
      real_t d = sqrt(4*a*c-b*b);
      return 2*(atan((2*a-b)/d)+atan((2*a+b)/d))/d;
   }

// Calculate the normal to sides and their lengths
   void _util();
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
