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

            Definition and Implementation of template class 'Point'

  ==============================================================================*/

#ifndef __POINT_H
#define __POINT_H

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

/*! \file Point.h
 *  \brief Definition file and implementation for class Point.
 */

/*! \class Point
 * \ingroup Util
 * \brief Defines a point with arbitrary type coordinates.
 * \details Operators <tt>=</tt> and <tt>()</tt> are overloaded.
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

template<class T_>
struct Point {

/// \brief Default constructor
    Point() { x = y = z = T_(0); }

/// \brief Constructor that assigns <tt>a</tt> to first coordinate.
    Point(T_ a) { x = a; y = z = T_(0); }

/// \brief Constructor that assigns <tt>a</tt> and <tt>b</tt> to first and second coordinates respectively.
    Point(T_ a, T_ b) { x = a; y = b; z = T_(0); }

/// \brief Constructor that assigns <tt>a</tt>, <tt>b</tt> and <tt>c</tt> to first, second and third coordinates respectively.
    Point(T_ a, T_ b, T_ c) { x = a; y = b; z = c; }

/// \brief Copy constructor
    Point(const Point<T_>& p) { x = p.x; y = p.y; z = p.z; }

/// \brief Operator <tt>()</tt>: Non constant version.
/// \details Values <tt>i = 1, 2, 3</tt> correspond to <tt>x</tt>, <tt>y</tt> and <tt>z</tt> respectively
    T_ &operator()(size_t i) 
    {
       switch (i) {
         case 1: return x;
         case 2: return y;
         case 3: return z;
       }
       return x;
    }

/// \brief Operator <tt>()</tt>: Constant version.
/// \details Values <tt>i = 1, 2, 3</tt> correspond to <tt>x</tt>, <tt>y</tt> and <tt>z</tt> 
/// respectively
    const T_ &operator()(size_t i) const
    {
       switch (i) {
         case 1: return x;
         case 2: return y;
         case 3: return z;
       }
       return x;
    }

/// \brief Operator <tt>[]</tt>: Non constant version.
/// \details Values <tt>i = 0, 1, 2</tt> correspond to <tt>x</tt>, <tt>y</tt> and <tt>z</tt> respectively
    T_ &operator[](size_t i)
    {
       switch (i) {
         case 0: return x;
         case 1: return y;
         case 2: return z;
       }
       return x;
    }

/// \brief Operator <tt>[]</tt>: Constant version.
/// \details Values <tt>i = 0, 1, 2</tt> correspond to <tt>x</tt>, <tt>y</tt> and <tt>z</tt> respectively
    const T_ &operator[](size_t i) const
    {
       switch (i) {
         case 0: return x;
         case 1: return y;
         case 2: return z;
       }
       return x;
    }

/// \brief Operator <tt>+=</tt>
/// \details Add point <tt>p</tt> to current instance
    Point<T_> & operator+=(const Point<T_>& p) { x+=p.x; y+=p.y; z+=p.z; return *this; }

/// \brief Operator <tt>-=</tt>
/// \details Subtract point <tt>p</tt> from current instance
    Point<T_> & operator-=(const Point<T_>& p) { x-=p.x; y-=p.y; z-=p.z; return *this; }

/// \brief Operator <tt>=</tt>
/// \details Assign constant <tt>a</tt> to current instance coordinates
    Point<T_> & operator=(const T_& a)  { x = y = z = a; return *this; }

/// \brief Operator <tt>+=</tt>
/// \details Add constant <tt>a</tt> to current instance coordinates
    Point<T_> & operator+=(const T_& a) { x+=a; y+=a; z+=a; return *this; }

/// \brief Operator <tt>-=</tt>
/// \details Subtract constant <tt>a</tt> from current instance coordinates
    Point<T_> & operator-=(const T_& a) { x-=a; y-=a; z-=a; return *this; }

/// \brief Operator <tt>*=</tt>
/// \details Multiply constant <tt>a</tt> by current instance coordinates
    Point<T_> & operator*=(const T_& a) { x*=a; y*=a; z*=a; return *this; }

/// \brief Operator <tt>/=</tt>
/// \details Divide current instance coordinates by <tt>a</tt>
    Point<T_> & operator/=(const T_& a) { x/=a; y/=a; z/=a; return *this; }

/// \brief Operator <tt>==</tt>
/// \details Return <tt>true</tt> if current instance is equal to <tt>p</tt>,
/// <tt>false</tt> otherwise.
    bool operator==(const Point<T_>& p)
    {
       if (p.x==x && p.y==y && p.z==z)
          return true;
       return false;
    }

/// \brief Operator <tt>!=</tt>
/// \details Return <tt>false</tt> if current instance is equal to <tt>p</tt>,
/// <tt>true</tt> otherwise.
    bool operator!=(const Point<T_>& p)
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
       x /= n, y /= n, z /= n;
    }

/// \brief Return Director (Normalized vector)
    Point<real_t> Director(const Point<real_t>& p) const
    {
       real_t n=p.Norm();
       return Point<real_t>(p.x/n,p.y/n,p.z/n);
    }

/// \brief Return <tt>true</tt> if current point is close to instance <tt>a</tt> 
/// (up to tolerance <tt>toler</tt>)
/// \details Default value for <tt>toler</tt> is the <tt>OFELI_TOLERANCE</tt> constant.
    bool isCloseTo(const Point<real_t>& a,
                   real_t               toler=OFELI_TOLERANCE) const
    {
       return areClose(*this,a,toler);
    }

/** \brief Return Dot (scalar) product of two vectors
 *  \details A typical use of this operator is <tt>double a = (p,q)</tt>
 *  where <tt>p</tt> and <tt>q</tt> are 2 instances of <tt>Point<double></tt>
 *  @param [in] p Point instance by which the current instance is multiplied
 */
    T_ operator,(const Point<T_>& p) const
    {
       return (x*p.x + y*p.y + z*p.z);
    }

/// \brief Point coordinates
    T_ x, y, z;

};


///////////////////////////////////////////////////////////////////////////////
//                             ASSOCIATED  FUNCTIONS                         //
///////////////////////////////////////////////////////////////////////////////

/** \fn bool operator== (const Point<T_> &a, const Point<T_> &b)
 *  \brief Operator <tt>==</tt>
 *  \ingroup Util
 *  \details Return <tt>true</tt> if <tt>a=b</tt>, <tt>false</tt> if not.
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
bool operator== (const Point<T_>& a,
                 const Point<T_>& b)
{ 
   return (a.x==b.x && a.y==b.y && a.z==b.z); 
}

/** \fn Point<T_> operator+ (const Point<T_> &a, const Point<T_> &b)
 *  \brief Operator <tt>+</tt>
 *  \ingroup Util
 *  \details Return sum of two points <tt>a</tt> and <tt>b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
Point<T_> operator+ (const Point<T_>& a,
                     const Point<T_>& b)
{
   return Point<T_>(a.x+b.x,a.y+b.y,a.z+b.z);
}

/** \fn Point<T_> operator+ (const Point<T_> &a, const T_ &x)
 *  \brief Operator <tt>+</tt>
 *  \ingroup Util
 *  \details Translate <tt>a</tt> by <tt>x</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
Point<T_> operator+ (const Point<T_>& a,
                     const T_&        x)
    { return Point<T_>(a.x+x,a.y+x,a.z+x); }


/** \fn Point<T_> operator- (const Point<T_> &a)
 * \brief Unary Operator <tt>-</tt>
 * \ingroup Util
 * \details Return minus <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
Point<T_> operator- (const Point<T_>& a) { return Point<T_>(-a.x,-a.y,-a.z); }


/** \fn Point<T_> operator- (const Point<T_> &a, const Point<T_> &b)
 *  \brief Operator <tt>-</tt>
 *  \ingroup Util
 *  \details Return point <tt>a</tt> minus point <tt>b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
Point<T_> operator- (const Point<T_>& a,
                     const Point<T_>& b)
   { return Point<T_>(a.x-b.x,a.y-b.y,a.z-b.z); }


/** \fn Point<T_> operator- (const Point<T_> &a, const T_ &x)
 *  \brief Operator <tt>-</tt>
 *  \ingroup Util
 *  \details Translate <tt>a</tt> by <tt>-x</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
Point<T_> operator- (const Point<T_>& a,
                     const T_&        x)
   { return Point<T_>(a.x-x,a.y-x,a.z-x); }

/** \fn Point<T_> operator* (const T_ &a, const Point<T_> &b)
 *  \brief Operator <tt>*</tt>
 *  \ingroup Util
 *  \details Return point <tt>b</tt> premultiplied by constant <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
Point<T_> operator* (const T_&        a,
                     const Point<T_>& b)
    { return Point<T_>(a*b.x,a*b.y,a*b.z); }


/** \fn Point<T_> operator* (const int &a, const Point<T_> &b)
 *  \brief Operator *.
 *  \ingroup Util
 *  \details Return point <tt>b</tt> divided by integer constant <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
Point<T_> operator* (const int&       a,
                     const Point<T_>& b)
     { return Point<T_>(a*b.x,a*b.y,a*b.z); }


/** \fn Point<T_> operator* (const Point<T_> &b, const T_ &a)
 *  \brief Operator <tt>/</tt>
 *  \ingroup Util
 *  \details Return point <tt>b</tt> multiplied by constant <tt>a</tt>
 */
template <class T_>
Point<T_> operator* (const Point<T_>& b,
                     const T_&        a)
      { return Point<T_>(a*b.x,a*b.y,a*b.z);}


/** \fn Point<T_> operator* (const Point<T_> &b, const int &a)
 *  \brief Operator <tt>*</tt>
 *  \ingroup Util
 *  \details Return point <tt>b</tt> postmultiplied by constant <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
Point<T_> operator* (const Point<T_>& b,
                     const int&       a)
          { return Point<T_>(a*b.x,a*b.y,a*b.z);}

/** \fn T_ operator* (const Point<T_> &b, const Point<T_> &a)
 * \brief Operator <tt>*</tt>
 * \ingroup Util
 * \details Return inner (scalar) product of points <tt>a</tt> and <tt>b</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
T_ operator* (const Point<T_>& a,
              const Point<T_>& b)
   { return (a.x*b.x+a.y*b.y+a.z*b.z); }

/** \fn Point<T_> operator/ (const Point<T_> &b, const T_ &a)
 * \brief Operator <tt>/</tt>
 * \ingroup Util
 * \details Return point <tt>b</tt> divided by constant <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
Point<T_> operator/ (const Point<T_>& b,
                     const T_&        a) 
          { return Point<T_>(b.x/a,b.y/a,b.z/a);}

/// \fn Point<double> CrossProduct(const Point<double> &lp, const Point<double> &rp)
/// \brief Return Cross product of two vectors <tt>lp</tt> and <tt>rp</tt>
inline Point<real_t> CrossProduct(const Point<real_t>& lp,
                                  const Point<real_t>& rp)
{
   return Point<real_t>(lp.y*rp.z-lp.z*rp.y,lp.z*rp.x-lp.x*rp.z,lp.x*rp.y-lp.y*rp.x);
}

/** \fn bool areClose(const Point<double> &a, const Point<double> &b, double toler=OFELI_EPSMCH)
 *  \brief Return <tt>true</tt> if both instances of class Point<double> are distant with less then <tt>toler</tt>
 *  \ingroup Util
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline bool areClose(const Point<real_t>& a,
                     const Point<real_t>& b,
                     real_t               toler=OFELI_TOLERANCE)
{
   if (((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z)) < toler*toler)
      return true;
   else
      return false;
}

/** \fn double SqrDistance(const Point<double>& a, const Point<double>& b)
 *  \brief Return squared euclidean distance between points <tt>a</tt> and <tt>b</tt>
 *  \ingroup Util
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline real_t SqrDistance(const Point<real_t>& a,
                          const Point<real_t>& b)
{
   return ((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
}

/** \fn double Distance(const Point<double> &a, const Point<double> &b)
 *  \brief Return euclidean distance between points <tt>a</tt> and <tt>b</tt>
 *  \ingroup Util
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline real_t Distance(const Point<real_t>& a,
                       const Point<real_t>& b)
{
   return sqrt(SqrDistance(a,b));
}


/** \fn bool operator< (const Point<size_t> &a, const Point<size_t> &b)
 * \brief Comparison operator. Returns true if all components of first vector are lower than
 * those of second one
 * \ingroup Util
 * \details Return minus <tt>a</tt>
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
inline bool operator< (const Point<size_t>& a,
                       const Point<size_t>& b)
{
   if (a.x < b.x)
      return true;
   if (a.x==b.x && a.y<b.y)
      return true;
   if (a.x==b.x && a.y==b.y && a.z<b.z)
      return true;
   return false;
}


/** \fn ostream & operator<<(std::ostream &s, const Point<T_> &a)
 *  \brief Output point coordinates.
 *  \ingroup Util
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
template <class T_>
std::ostream & operator<<(std::ostream&    s,
                          const Point<T_>& a)
{
   s.setf(ios::scientific);
   s << "( " << std::setprecision(8) << setw(12) << a.x << " , ";
   s << setw(12) << a.y << " , " << setw(12) << a.z << " )";
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
