/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

               Definition and Implementation of class 'SpaceTime'

  ==============================================================================*/

#ifndef __SPACE_TIME_H
#define __SPACE_TIME_H

#include "OFELI_Config.h"

#include <iostream>
using std::ostream;

#include <iomanip>
using std::setw;
using std::ios;

#include <math.h>

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file SpaceTime.h
 *  \brief Definition file and implementation for class SpaceTime.
 */

/*! \class SpaceTime
 * \ingroup Util
 * \brief Defines a space-time point.
 * \details Operators <tt>=</tt> and <tt>()</tt> are overloaded.
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

struct SpaceTime {

/// \brief Default constructor
    SpaceTime() { x = y = z = t = 0.0; }

/// \brief Constructor that assigns <tt>a</tt>, <tt>b</tt>, <tt>c</tt> and <tt>d</tt> to first, second, third coordinates and time respectively.
    SpaceTime(real_t a, real_t b, real_t c, real_t d) { x = a; y = b; z = c; t = d; }

/// \brief Copy constructor
    SpaceTime(const SpaceTime& p) { x = p.x; y = p.y; z = p.z; t = p.t; }

/// \brief Operator <tt>()</tt>: Non constant version.
/// \details Values <tt>i = 1, 2, 3, 4</tt> correspond to <tt>x</tt>, <tt>y</tt>, <tt>z</tt> 
/// and <tt>t</tt> respectively
    real_t &operator()(size_t i) 
    {
       switch (i) {
          case 1:  return x;
          case 2:  return y;
          case 3:  return z;
          case 4:  return t;
          default: return x;
       }
       return x;
    }

/// \brief Operator <tt>()</tt>: Constant version.
/// \details Values <tt>i = 1, 2, 3, 4</tt> correspond to <tt>x</tt>, <tt>y</tt>, <tt>z</tt> 
/// and <tt>t</tt> respectively
    const real_t &operator()(size_t i) const
    {
       switch (i) {
          case 1:  return x;
          case 2:  return y;
          case 3:  return z;
          case 4:  return t;
          default: return x;
       }
       return x;
    }

/// \brief Operator <tt>[]</tt>: Non constant version.
/// \details Values <tt>i = 0, 1, 2, 3</tt> correspond to <tt>x</tt>, <tt>y</tt>, <tt>z</tt>
/// and <tt>t</tt> respectively
    real_t &operator[](size_t i)
    {
       switch (i) {
          case 0:  return x;
          case 1:  return y;
          case 2:  return z;
          case 3:  return t;
          default: return x;
       }
       return x;
    }

/// \brief Operator <tt>[]</tt>: Constant version.
/// \details Values <tt>i = 0, 1, 2</tt> correspond to <tt>x</tt>, <tt>y</tt>, <tt>z</tt> 
/// and <tt>t</tt> respectively
    const real_t &operator[](size_t i) const
    {
       switch (i) {
          case 0:  return x;
          case 1:  return y;
          case 2:  return z;
          case 3:  return t;
          default: return x;
       }
       return x;
    }

/// \brief Operator <tt>+=</tt>
/// \details Add point <tt>p</tt> to current instance
    SpaceTime & operator+=(const SpaceTime& p) { x+=p.x; y+=p.y; z+=p.z; t+=p.t; return *this; }

/// \brief Operator <tt>-=</tt>
/// \details Subtract point <tt>p</tt> from current instance
    SpaceTime & operator-=(const SpaceTime& p) { x-=p.x; y-=p.y; z-=p.z; t-=p.t; return *this; }

/// \brief Operator <tt>+=</tt>
/// \details Add constant <tt>a</tt> to current instance coordinates
    SpaceTime & operator+=(const real_t& a) { x+=a; y+=a; z+=a; t+= a; return *this; }

/// \brief Operator <tt>-=</tt>
/// \details Subtract constant <tt>a</tt> from current instance coordinates
    SpaceTime & operator-=(const real_t& a) { x-=a; y-=a; z-=a; t-=a; return *this; }

/// \brief Operator <tt>*=</tt>
/// \details Multiply constant <tt>a</tt> by current instance coordinates
    SpaceTime & operator*=(const real_t& a) { x*=a; y*=a; z*=a; t*=a; return *this; }

/// \brief Operator <tt>/=</tt>
/// \details Divide current instance coordinates by <tt>a</tt>
    SpaceTime & operator/=(const real_t& a) { x/=a; y/=a; z/=a; t/=a; return *this; }

/// \brief Operator <tt>==</tt>
/// \details Return <tt>true</tt> if current instance is equal to <tt>p</tt>,
/// <tt>false</tt> otherwise.
    bool operator==(const SpaceTime& p)
    {
       if (p.x==x && p.y==y && p.z==z && p.t==t)
          return true;
       return false;
    }

/// \brief Operator <tt>!=</tt>
/// \details Return <tt>false</tt> if current instance is equal to <tt>p</tt>,
/// <tt>true</tt> otherwise.
    bool operator!=(const SpaceTime& p)
    {
       return !(*this==p);
    }

/// \brief Return squared euclidean norm of vector
    real_t NNorm() const { return (x*x + y*y + z*z); }

/// \brief Return norm (length) of vector
    real_t Norm() const { return sqrt(NNorm()); }

/// \brief Normalize vector.
/// \details Divide vector components by its 2-norm
    void Normalize()
    {
       real_t n=Norm();
       if (n==0)
          return;
       x /= n;
       y /= n;
       z /= n;
    }

/// \brief Return Director (Normalized vector)
    SpaceTime Director(const SpaceTime& p) const
    {
       real_t n=p.Norm();
       return SpaceTime(p.x/n,p.y/n,p.z/n,0.);
    }

/// \brief Return <tt>true</tt> if current point is close to instance <tt>a</tt> 
/// (up to tolerance <tt>toler</tt>)
/// \details Default value for <tt>toler</tt> is the <tt>OFELI_TOLERANCE</tt> constant.
    bool isCloseTo(const SpaceTime& a,
                   real_t           toler=OFELI_TOLERANCE) const
    {
       return (((x-a.x)*(x-a.x)+(y-a.y)*(y-a.y)+(z-a.z)*(z-a.z)+(t-a.t)*(t-a.t)) < toler*toler);
    }

/** \brief Return Dot (scalar) product of two vectors
 *  \details A typical use of this operator is <tt>double a = (p,q)</tt>
 *  where <tt>p</tt> and <tt>q</tt> are 2 instances of <tt>SpaceTime</tt>
 *  @param [in] p Point instance by which the current instance is multiplied
 */
    real_t operator,(const SpaceTime& p) const
    {
       return (x*p.x + y*p.y + z*p.z);
    }

/// \brief Coordinates
    real_t x, y, z, t;

};


///////////////////////////////////////////////////////////////////////////////
//                             ASSOCIATED  FUNCTIONS                         //
///////////////////////////////////////////////////////////////////////////////

/** \fn bool operator== (const SpaceTime &a, const SpaceTime &b)
 *  \brief Operator <tt>==</tt>
 *  \ingroup Util
 *  \details Return <tt>true</tt> if <tt>a=b</tt>, <tt>false</tt> if not.
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline bool operator== (const SpaceTime& a,
                        const SpaceTime& b)
{ 
   return (a.x==b.x && a.y==b.y && a.z==b.z && a.t==b.t); 
}

/** \fn SpaceTime operator+ (const SpaceTime &a, const SpaceTime &b)
 *  \brief Operator <tt>+</tt>
 *  \ingroup Util
 *  \details Return sum of two points <tt>a</tt> and <tt>b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline SpaceTime operator+ (const SpaceTime& a,
                            const SpaceTime& b)
{
   return SpaceTime(a.x+b.x,a.y+b.y,a.z+b.z,a.t+b.t);
}

/** \fn SpaceTime operator+ (const SpaceTime &a, const real_t &x)
 *  \brief Operator <tt>+</tt>
 *  \ingroup Util
 *  \details Translate <tt>a</tt> by <tt>x</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline SpaceTime operator+ (const SpaceTime& a,
                            const real_t&    x)
{
   return SpaceTime(a.x+x,a.y+x,a.z+x,a.t+x);
}


/** \fn SpaceTime operator- (const SpaceTime &a)
 * \brief Unary Operator <tt>-</tt>
 * \ingroup Util
 * \details Return minus <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline SpaceTime operator- (const SpaceTime& a)
{
   return SpaceTime(-a.x,-a.y,-a.z,-a.t);
}


/** \fn SpaceTime operator- (const SpaceTime &a, const SpaceTime &b)
 *  \brief Operator <tt>-</tt>
 *  \ingroup Util
 *  \details Return point <tt>a</tt> minus point <tt>b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline SpaceTime operator- (const SpaceTime& a,
                            const SpaceTime& b)
{
   return SpaceTime(a.x-b.x,a.y-b.y,a.z-b.z,a.t-b.t);
}


/** \fn SpaceTime operator- (const SpaceTime &a, const real_t &x)
 *  \brief Operator <tt>-</tt>
 *  \ingroup Util
 *  \details Translate <tt>a</tt> by <tt>-x</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline SpaceTime operator- (const SpaceTime& a,
                            const real_t&    x)
{
   return SpaceTime(a.x-x,a.y-x,a.z-x,a.t-x);
}


/** \fn SpaceTime operator* (const real_t &a, const SpaceTime &b)
 *  \brief Operator <tt>*</tt>
 *  \ingroup Util
 *  \details Return point <tt>b</tt> premultiplied by constant <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline SpaceTime operator* (const real_t&    a,
                            const SpaceTime& b)
{
   return SpaceTime(a*b.x,a*b.y,a*b.z,a*b.t);
}


/** \fn SpaceTime operator* (const int &a, const SpaceTime &b)
 *  \brief Operator *.
 *  \ingroup Util
 *  \details Return point <tt>b</tt> divided by integer constant <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline SpaceTime operator* (const int&       a,
                            const SpaceTime& b)
{
   return SpaceTime(a*b.x,a*b.y,a*b.z,a*b.t);
}


/** \fn SpaceTime operator* (const SpaceTime &b, const real_t &a)
 *  \brief Operator <tt>/</tt>
 *  \ingroup Util
 *  \details Return point <tt>b</tt> multiplied by constant <tt>a</tt>
 */
inline SpaceTime operator* (const SpaceTime& b,
                            const real_t&    a)
{
   return SpaceTime(a*b.x,a*b.y,a*b.z,a*b.t);
}


/** \fn SpaceTime operator* (const SpaceTime &b, const int &a)
 *  \brief Operator <tt>*</tt>
 *  \ingroup Util
 *  \details Return point <tt>b</tt> postmultiplied by constant <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline SpaceTime operator* (const SpaceTime& b,
                            const int&       a)
{
   return SpaceTime(a*b.x,a*b.y,a*b.z,a*b.t);
}

/** \fn real_t operator* (const SpaceTime &b, const SpaceTime &a)
 * \brief Operator <tt>*</tt>
 * \ingroup Util
 * \details Return inner (scalar) product of points <tt>a</tt> and <tt>b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline real_t operator* (const SpaceTime& a,
                         const SpaceTime& b)
{
   return (a.x*b.x+a.y*b.y+a.z*b.z+b.t*b.t);
}

/** \fn SpaceTime operator/ (const SpaceTime &b, const real_t &a)
 * \brief Operator <tt>/</tt>
 * \ingroup Util
 * \details Return point <tt>b</tt> divided by constant <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline SpaceTime operator/ (const SpaceTime& b,
                            const real_t&    a) 
{
   return SpaceTime(b.x/a,b.y/a,b.z/a,b.t/a);
}

/// \fn SpaceTime CrossProduct(const SpaceTime &lp, const SpaceTime &rp)
/// \brief Return Cross product of two vectors <tt>lp</tt> and <tt>rp</tt>
inline SpaceTime CrossProduct(const SpaceTime& lp,
                              const SpaceTime& rp)
{
   return SpaceTime(lp.y*rp.z-lp.z*rp.y,lp.z*rp.x-lp.x*rp.z,lp.x*rp.y-lp.y*rp.x,0.);
}

/** \fn bool areClose(const SpaceTime &a, const SpaceTime &b, double toler=OFELI_EPSMCH)
 *  \brief Return <tt>true</tt> if both instances of class Point<double> are distant with less then <tt>toler</tt>
 *  \ingroup Util
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline bool areClose(const SpaceTime& a,
                     const SpaceTime& b,
                     real_t           toler=OFELI_TOLERANCE)
{
   if (((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z)+(a.t-b.t)*(a.t-b.t)) < toler*toler)
      return true;
   else
      return false;
}

/** \fn double SqrDistance(const SpaceTime& a, const SpaceTime& b)
 *  \brief Return squared euclidean distance between points <tt>a</tt> and <tt>b</tt>
 *  \ingroup Util
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline real_t SqrDistance(const SpaceTime& a,
                          const SpaceTime& b)
{
   return ((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z) + (a.t-b.t)*(a.t-b.t));
}

/** \fn double Distance(const SpaceTime &a, const SpaceTime &b)
 *  \brief Return euclidean distance between points <tt>a</tt> and <tt>b</tt>
 *  \ingroup Util
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline real_t Distance(const SpaceTime& a,
                       const SpaceTime& b)
{
   return sqrt(SqrDistance(a,b));
}


/** \fn bool operator< (const SpaceTime &a, const SpaceTime &b)
 * \brief Comparison operator. Returns true if all components of first vector are lower than
 * those of second one
 * \ingroup Util
 * \details Return minus <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline bool operator< (const SpaceTime& a,
                       const SpaceTime& b)
{
   if (a.x==b.x && a.y==b.y && a.z<b.z)
      return true;
   return false;
}


/** \fn ostream & operator<<(std::ostream &s, const SpaceTime &a)
 *  \brief Output space-time point coordinates.
 *  \ingroup Util
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline std::ostream & operator<<(std::ostream&    s,
                                 const SpaceTime& a)
{
   s.setf(ios::scientific);
   s << "( " << std::setprecision(8) << setw(12) << a.x << " , ";
   s << setw(12) << a.y << " , " << setw(12) << a.z  << " , " << setw(12) << a.t << " )";
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
