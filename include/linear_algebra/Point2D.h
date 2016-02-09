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

                Definition and Implementation of template class 'Point2D'

  ==============================================================================*/

#ifndef __POINT_2D_H
#define __POINT_2D_H

#include "OFELI_Config.h"
#include "Point.h"

namespace OFELI {

/*! \file Point2D.h
 *  \brief Definition file for class Point2D.
 */

/*! \class Point2D
 * \ingroup Util
 *  \brief Defines a 2-D point with arbitrary type coordinates.
 *
 * \details Operators <tt>=</tt> and <tt>()</tt> are overloaded. The actual 
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 */

template<class T_>
struct Point2D {

/// \brief Default constructor
    Point2D() { x = y = T_(0); }

/// \brief Constructor that assigns <tt>a</tt>, <tt>b</tt> to <tt>x</tt>-, <tt>y</tt>- 
/// and <tt>y</tt>-coordinates respectively.
/// \details Default value for \a b is 0
    Point2D(T_ a, T_ b=T_(0)) { x = a; y = b; }

/// \brief Initialize point coordinates with C-array \a a
    Point2D(T_ *a) { x = a[0]; y = a[1]; }

/// \brief Copy constructor
    Point2D(const Point2D<T_>& pt) { x = pt.x; y = pt.y; }

/// \brief Copy constructor from class Point
    Point2D(const Point<T_>& pt) { x = pt.x; y = pt.y; }

/// \brief Operator() : Non constant version.
/// \details Values <tt>i = 1,2</tt> correspond to <tt>x</tt> and <tt>y</tt> respectively
    T_ &operator()(size_t i) 
    {
       if (i==2)
          return y;
       else
          return x;
    }

/// \brief Operator() : Constant version.
/// \details Values <tt>i=1,2</tt> correspond to <tt>x</tt> and <tt>y</tt> respectively
    const T_ &operator()(size_t i) const
    {
       if (i==2)
          return y;
       else
          return x;
    }

/// \brief Operator<tt>[]</tt>: Non constant version.
/// \details Values <tt>i=0,1</tt> correspond to <tt>x</tt> and <tt>y</tt> respectively
    T_ &operator[](size_t i)
    {
       if (i==1)
          return y;
       else
          return x;
    }

/// \brief Operator<tt>[]</tt> Constant version
/// \details Values <tt>i=0,1</tt> correspond to <tt>x</tt> and <tt>y</tt> respectively
    const T_ &operator[](size_t i) const
    {
       if (i==1)
          return y;
       else
          return x;
    }

/// \brief Operator <tt>=</tt>
/// \details Assign point <tt>p</tt> to current instance
    Point2D<T_> & operator=(const Point2D<T_>& p) { x = p.x; y = p.y; return *this; }

/// \brief Operator <tt>+=</tt>
/// \details Add point <tt>p</tt> to current instance
    Point2D<T_> & operator+=(const Point2D<T_>& p) { x+=p.x; y+=p.y; return *this; }

/// \brief Operator <tt>-=</tt>
/// \details Subtract point <tt>p</tt> from current instance
    Point2D<T_> & operator-=(const Point2D<T_>& p) { x-=p.x; y-=p.y; return *this; }

/// \brief Operator <tt>=</tt>
/// \details Assign constant <tt>a</tt> to current instance coordinates
    Point2D<T_> & operator=(const T_& a)  { x = y = a; return *this; }

/// \brief Operator <tt>+=</tt>
/// \details Add constant <tt>a</tt> to current instance coordinates
    Point2D<T_> & operator+=(const T_& a) { x+=a; y+=a; return *this; }

/// \brief Operator <tt>-=</tt>
/// \details Subtract constant <tt>a</tt> from current instance coordinates
    Point2D<T_> & operator-=(const T_& a) { x-=a; y-=a; return *this; }

/// \brief Operator <tt>*=</tt>
/// \details Multiply constant <tt>a</tt> by current instance coordinates
    Point2D<T_> & operator*=(const T_& a) { x*=a; y*=a; return *this; }

/// \brief Operator <tt>/=</tt>
/// \details Divide current instance coordinates by <tt>a</tt>
    Point2D<T_> & operator/=(const T_& a) { x/=a; y/=a; return *this; }

/// \brief Operator <tt>==</tt>
/// \details Return <tt>true</tt> if current instance is equal to <tt>p</tt>,
/// <tt>false</tt> otherwise.
    bool operator==(const Point2D<T_>& p)
    {
       if (p.x!=x && p.y!=y)
          return false;
       return true;
    }

/// \brief Operator <tt>!=</tt>
/// \details Return <tt>false</tt> if current instance is equal to <tt>p</tt>,
/// <tt>true</tt> otherwise.
    bool operator!=(const Point2D<T_>& p)
    {
       if (p.x!=x && p.y!=y)
          return true;
       return false;
    }

/// \brief Return Cross product of two vectors <tt>lp</tt> and <tt>rp</tt>
    real_t CrossProduct(const Point2D<real_t>& lp,
                        const Point2D<real_t>& rp)
    {
       return (lp.x*rp.y - lp.y*rp.x);
    }

/// \brief Return squared norm (length) of vector
    real_t NNorm() const { return (x*x + y*y); }

/// \brief Return norm (length) of vector
    real_t Norm() const { return sqrt(NNorm()); }

/// \brief Return Director (Normalized vector)
    Point2D<real_t> Director(const Point2D<real_t>& p) const
    {
       real_t n = p.Norm();
       return Point2D<real_t>(p.x/n,p.y/n);
    }

/// \brief Return <tt>true</tt> if current point is close to instance <tt>a</tt> 
/// (up to tolerance <tt>toler</tt>)
    bool isCloseTo(const Point2D<real_t>& a,
                         real_t           toler=OFELI_TOLERANCE) const
    {
       return (*this,a,toler);
    }

/// \brief First coordinate of point
    T_ x;

/// \brief Second coordinate of point
    T_ y;

};


//-----------------------------------------------------------------------------
// Associated functions
//-----------------------------------------------------------------------------

/// \fn bool operator== (const Point2D<T_>& a, const Point2D<T_>& b)
/// \brief Operator ==.
/// \ingroup Util
/// \details Return <tt>true</tt> if <tt>a=b</tt>, false if not.
template <class T_>
bool operator== (const Point2D<T_>& a,
                 const Point2D<T_>& b)
{ 
   return (a.x==b.x && a.y==b.y);
}

/// \fn Point2D<T_> operator+ (const Point2D<T_>& a, const Point2D<T_>& b)
/// \brief Operator +.
/// \ingroup Util
/// \details Return sum of two points <tt>a</tt> and <tt>b</tt>
template <class T_>
Point2D<T_> operator+ (const Point2D<T_>& a,
                       const Point2D<T_>& b)
   { return Point2D<T_>(a.x+b.x,a.y+b.y); }

/// \fn Point2D<T_> operator+ (const Point2D<T_> &a, const T_ &x)
/// \brief Operator +.
/// \ingroup Util
/// \details Translate <tt>a</tt> by <tt>x</tt>
template <class T_>
Point2D<T_> operator+ (const Point2D<T_>& a,
                       const T_&          x)
    { return Point2D<T_>(a.x+x,a.y+x); }


/// \fn Point2D<T_> operator- (const Point2D<T_>& a)
/// \brief Unary Operator <tt>-</tt>
/// \ingroup Util
/// \details Return minus <tt>a</tt>
template <class T_>
Point2D<T_> operator- (const Point2D<T_>& a) { return Point2D<T_>(-a.x,-a.y); }


/// \fn Point2D<T_> operator- (const Point2D<T_>& a, const Point2D<T_>& b)
/// \brief Operator <tt>-</tt>
/// \ingroup Util
/// \details Return point <tt>a</tt> minus point <tt>b</tt>
template <class T_>
Point2D<T_> operator- (const Point2D<T_>& a,
                       const Point2D<T_>& b)
    { return Point2D<T_>(a.x-b.x,a.y-b.y); }


/// \fn Point2D<T_> operator- (const Point2D<T_>& a, const T_& x)
/// \brief Operator <tt>-</tt>
/// \ingroup Util
/// \details Translate <tt>a</tt> by <tt>-x</tt>
template <class T_>
Point2D<T_> operator- (const Point2D<T_>& a,
                       const T_&          x)
    { return Point2D<T_>(a.x-x,a.y-x); }


/// \fn Point2D<T_> operator* (const T_ &a, const Point2D<T_> &b)
/// \brief Operator *.
/// \ingroup Util
/// \details Return point <tt>b</tt> premultiplied by constant <tt>a</tt>
template <class T_>
Point2D<T_> operator* (const T_&          a,
                       const Point2D<T_>& b)
   { return Point2D<T_>(a*b.x,a*b.y); }


/// \fn Point2D<T_> operator* (const int& a, const Point2D<T_>& b)
/// Operator *.
/// \ingroup Util
/// Return point <tt>b</tt> divided by integer constant <tt>a</tt>
template <class T_>
Point2D<T_> operator* (const int&         a,
                       const Point2D<T_>& b)
   { return Point2D<T_>(a*b.x,a*b.y); }


/// \fn Point2D<T_> operator* (const Point2D<T_>& b, const T_& a)
/// \brief Operator <tt>/</tt>
/// \ingroup Util
/// \details Return point <tt>b</tt> postmultiplied by constant <tt>a</tt>
template <class T_>
Point2D<T_> operator* (const Point2D<T_>& b,
                       const T_&          a)
   { return Point2D<T_>(a*b.x,a*b.y);}


/// \fn Point2D<T_> operator* (const Point2D<T_> &b, const int &a)
/// \brief Operator <tt>*</tt>
/// \ingroup Util
/// \details Return point <tt>b</tt> postmultiplied by constant <tt>a</tt>
template <class T_>
Point2D<T_> operator* (const Point2D<T_>& b,
                       const int&         a)
   { return Point2D<T_>(a*b.x,a*b.y);}

/// \fn T_ operator* (const Point2D<T_>& b, const Point2D<T_>& a)
/// \brief Operator *.
/// \ingroup Util
/// \details Return point \a b postmultiplied by integer constant \a a.
template <class T_>
T_ operator* (const Point2D<T_>& b,
              const Point2D<T_>& a)
   { return (a.x*b.x+a.y*b.y); }

/// \fn Point2D<T_> operator/ (const Point2D<T_>& b, const T_& a)
/// \brief Operator <tt>/</tt>
/// \ingroup Util
/// \details Return point <tt>b</tt> divided by constant <tt>a</tt>
template <class T_>
Point2D<T_> operator/ (const Point2D<T_>& b,
                       const T_&          a)
   { return Point2D<T_>(b.x/a,b.y/a);}

/** \fn bool areClose(const Point2D<real_t> &a, const Point2D<real_t> &b, real_t toler=OFELI_EPSMCH)
 *  \brief Return <tt>true</tt> if both instances of class \b Point2D<real_t> are distant with
 *  less then \a toler [Default: <tt>OFELI_EPSMCH</tt>].
 * \ingroup Util
 */
inline bool areClose(const Point2D<real_t>& a,
                     const Point2D<real_t>& b,
                           real_t           toler=OFELI_TOLERANCE)
{
   if (((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)) < toler*toler)
      return true;
   else
      return false;
}


/// \fn real_t SqrDistance(const Point2D<real_t> &a, const Point2D<real_t> &b)
/// \brief Return squared euclidean distance between points <tt>a</tt> and <tt>b</tt>
/// \ingroup Util
inline real_t SqrDistance(const Point2D<real_t>& a,
                          const Point2D<real_t>& b)
{
   return ((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

/// \fn real_t Distance(const Point2D<real_t>& a, const Point2D<real_t>& b)
/// \brief Return euclidean distance between points <tt>a</tt> and <tt>b</tt>
/// \ingroup Util
inline real_t Distance(const Point2D<real_t>& a,
                       const Point2D<real_t>& b)
{
   return sqrt(SqrDistance(a,b));
}

/// \fn ostream & operator<<(std::ostream& s, const Point2D<T_>& a)
/// \brief Output point coordinates.
/// \ingroup Util
template <class T_>
std::ostream & operator<<(      std::ostream& s,
                          const Point2D<T_>&  a)
{
   s.setf(ios::scientific);
   s << "( " << std::setprecision(8) << setw(12) << a.x << " , ";
   s << setw(12) << a.y << " , " << " )";
   return s;
}

} /* namespace OFELI */

#endif
