/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani
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

                              Some utility functions

  ==============================================================================*/

#ifndef __UTIL_H
#define __UTIL_H

#include <complex>
#include <string>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Point.h"
using std::string;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file util.h
 *  \brief File that contains various utility functions.
 */


/** \fn real_t Sgn(real_t a)
 * \ingroup Util
 * \brief Return sign of <tt>a</tt>: <tt>-1</tt> or <tt>1</tt>. 
 */
inline int Sgn(real_t a) { return (((a) > 0) ? (1) : -(1)); }


/** \fn real_t Sgn(real_t a, real_t b)
 * \ingroup Util
 * \brief Return \c |a| if \c b>0, \c -|a| if not
 */
inline real_t Sgn(real_t a, real_t b) { return ((b) >= 0.0 ? std::abs(a) : -std::abs(a)); }


/** \fn real_t Abs2(complex_t a)
 *  \ingroup Util
 *  \brief Return square of modulus of complex number <tt>a</tt>
 */
inline real_t Abs2(complex_t a) { return (a.real()*a.real()+a.imag()*a.imag()); }


/** \fn real_t Abs2(real_t a)
 *  \ingroup Util
 *  \brief Return square of real number <tt>a</tt>
 */
inline real_t Abs2(real_t a) { return (a*a); }


/** \fn real_t Abs(real_t a)
 * \ingroup Util
 * \brief Return absolute value of <tt>a</tt>
 */
inline real_t Abs(real_t a) { return (((a) > 0) ? (a) : -(a)); }


/** \fn real_t Abs(complex_t a)
 *  \ingroup Util
 *  \brief Return modulus of complex number <tt>a</tt>
 */
inline real_t Abs(complex_t a) { return sqrt(Abs2(a)); }


/** \fn real_t Abs(const Point<real_t>& p)
 *  \ingroup Util
 *  \brief Return norm of vector <tt>a</tt>
 */
inline real_t Abs(const Point<real_t>& p) { return p.Norm(); }


/** \fn real_t Conjg(real_t x)
 *  \ingroup Util
 *  \brief Return complex conjugate of real number <tt>a</tt>
 */
inline real_t Conjg(real_t a) { return a; }


/** \fn complex_t Conjg(complex_t x)
 *  \ingroup Util
 *  \brief Return complex conjugate of complex number <tt>a</tt>
 */
inline complex_t Conjg(complex_t a) { return complex_t(a.real(),-a.imag()); }


/** \fn string complex_string(real_t x, real_t y)
 *  \ingroup Util
 *  \brief Return string to conviently display a complex number
 *  @param [in] z Complex number
 */
inline string out_complex(complex_t z)
{
   std::ostringstream out;  
   if (z.imag()==0.0) {
      out << z.real();
      return out.str();
   }
   string ss = " + ";
   if (z.imag()<0.0)
      ss = " - ";
   out << z.real() << ss << std::abs(z.imag()) << " I";
   return out.str();
}


/** \fn string complex_string(real_t x, real_t y)
 *  \ingroup Util
 *  \brief Return string to conviently display a complex number
 *  @param [in] x Real part
 *  @param [in] y Imaginary part
 */
inline string out_complex(real_t x, real_t y)
{
   std::ostringstream out;  
   if (y==0.0) {
      out << x;
      return out.str();
   }
   string ss = " + ";
   if (y<0.0)
      ss = " - ";
   out << x << ss << std::abs(y) << " I";
   return out.str();
}


/** \fn real_t Max(real_t a, real_t b, real_t c)
 *  \ingroup Util
 *  \brief Return maximum value of real numbers <tt>a</tt>, <tt>b</tt> and <tt>c</tt>
 */
inline real_t Max(real_t a,
                  real_t b,
                  real_t c) { return std::max(std::max(a,b),c);}


/** \fn size_t Kronecker(int i, int j)
 * \ingroup Util
 * \brief Return Kronecker delta of <tt>i</tt> and <tt>j</tt>. 
 */
inline int Kronecker(int i, int j) { return ((i==j) ? 1 : 0); }


/** \fn const T_& Max(const T_& a, const T_& b)
 *  \ingroup Util
 *  \brief Return maximum value of elements <tt>a</tt> and <tt>b</tt>
 */
template<class T_>
const T_& Max(const T_& a,
              const T_& b)
{
   return (a < b) ? b : a;
}


/** \fn const T_& Min(const T_& a, const T_& b)
 *  \ingroup Util
 *  \brief Return minimum value of elements <tt>a</tt> and <tt>b</tt>
 */
template<class T_>
const T_& Min(const T_& a,
              const T_& b)
{
   return (a > b) ? b : a;
}


/** \fn int Max(int a, int b, int c)
 *  \ingroup Util
 *  \brief Return maximum value of integer numbers <tt>a</tt>, <tt>b</tt> and <tt>c</tt>
 */
inline int Max(int a,
               int b,
               int c) { return std::max(std::max(a,b),c); }


/** \fn real_t Min(real_t a, real_t b, real_t c)
 *  \ingroup Util
 *  \brief Return minimum value of real numbers <tt>a</tt>, <tt>b</tt> and <tt>c</tt>
 */
inline real_t Min(real_t a,
                  real_t b,
                  real_t c) { return std::min(std::min(a,b),c); }

/** \fn Min(int a, int b, int c)
 *  \ingroup Util
 *  \brief Return minimum value of integer numbers <tt>a</tt>, <tt>b</tt> and <tt>c</tt>
 */
inline int Min(int a,
               int b,
               int c) { return std::min(std::min(a,b),c); }


/** \fn real_t Max(real_t a, real_t b, real_t c, real_t d)
 *  \ingroup Util
 *  \brief Return maximum value of integer numbers <tt>a</tt>, <tt>b</tt>, <tt>c</tt>
 *  and <tt>d</tt>
 */
inline real_t Max(real_t a,
                  real_t b,
                  real_t c,
                  real_t d) { return std::max(Max(a,b,c),d); }


/** \fn int Max(int a, int b, int  c, int d)
 *  \ingroup Util
 * \brief Return maximum value of integer numbers <tt>a</tt>, <tt>b</tt>, <tt>c</tt>
 * and <tt>d</tt>
 */
inline int Max(int a,
               int b,
               int c,
               int d) { return std::max(Max(a,b,c),d); }


/** \fn real_t Min(real_t a, real_t b, real_t c, real_t d)
 *  \ingroup Util
 *  \brief Return minimum value of real numbers <tt>a</tt>, <tt>b</tt>, <tt>c</tt>
 *  and <tt>d</tt>
 */
inline real_t Min(real_t a,
                  real_t b,
                  real_t c,
                  real_t d) { return std::min(Min(a,b,c),d); }


/// \fn int Min(int a, int b, int c, int d)
/// \ingroup Util
/// \brief Return minimum value of integer numbers <tt>a</tt>, <tt>b</tt>, <tt>c</tt> and <tt>d</tt>
inline int Min(int a,
               int b,
               int c,
               int d) { return std::min(Min(a,b,c),d); }


/// \fn real_t Arg(complex_t x)
/// \ingroup Util
/// \brief Return argument of complex number <tt>x</tt>
inline real_t Arg(complex_t x)
{
  real_t a = 0.;
  if (x.real() == 0.)
    a = 0.5*OFELI_PI;
  else
    a = atan(x.imag()/x.real());
  return a;
}


/** \fn complex_t Log(complex_t x)
 *  \ingroup Util
 *  \brief Return principal determination of logarithm of complex number <tt>x</tt>
 */
inline complex_t Log(complex_t x) { return complex_t(log(Abs(x)),Arg(x)); }


/// \fn T_ Sqr(T_ x)
/// \ingroup Util
/// \brief Return square of value <tt>x</tt>
template<class T_>
inline T_ Sqr(T_ x) { return x*x; }


/** \fn void Scale(T_ a, const vector<T_>& x, vector<T_>& y)
 *  \ingroup Util
 *  \brief Mutiply vector <tt>x</tt> by <tt>a</tt> and save result in vector <tt>y</tt>
 *  \details <tt>x</tt> and <tt>y</tt> are instances of class vector<T_>
 */
template<class T_>
inline void Scale(T_                a,
                  const vector<T_>& x,
                  vector<T_>&       y)
{
   for (size_t i=0; i<x.size(); ++i)
      y[i] = a * x[i];
}


/** \fn void Scale(T_ a, const Vect<T_>& x, Vect<T_>& y)
 *  \ingroup Util
 *  \brief Mutiply vector <tt>x</tt> by <tt>a</tt> and save result in vector <tt>y</tt>
 *  \details <tt>x</tt> and <tt>y</tt> are instances of class Vect<T_>
 */
template<class T_>
inline void Scale(T_              a,
                  const Vect<T_>& x,
                  Vect<T_>&       y)
{
   for (size_t i=0; i<x.size(); ++i)
      y[i] = a * x[i];
}


/** \fn void Scale(T_ a, vector<T_>& x)
 *  \ingroup Util
 *  \brief Mutiply vector <tt>x</tt> by <tt>a</tt>
 *  \details <tt>x</tt> is an instance of class vector<T_>
 */
template<class T_>
inline void Scale(T_          a,
                  vector<T_>& x)
{
   for (size_t i=0; i<x.size(); ++i)
      x[i] *= a;
}


/** \fn void Xpy(size_t n, T_* x, T_* y)
 *  \ingroup Util
 * \brief Add array <tt>x</tt> to <tt>y</tt>
 */
template<class T_>
inline void Xpy(size_t n,
                T_*    x,
                T_*    y)
{
   for (size_t i=0; i<n; ++i)
      y[i] += x[i];
}


/** \fn void Xpy(const vector<T_>& x, vector<T_>& y)
 *  \ingroup Util
 *  \brief Add vector <tt>x</tt> to <tt>y</tt>
 *  \details <tt>x</tt> and <tt>y</tt> are instances of class vector<T_>
 */
template<class T_>
inline void Xpy(const vector<T_>& x,
                vector<T_>&       y)
{
   for (size_t i=0; i<x.size(); i++)
      y[i] += x[i];
}


/** \fn void Axpy(size_t n, T_ a, T_* x, T_* y)
 *  \ingroup Util
 *  \brief Multiply array <tt>x</tt> by <tt>a</tt> and add result to <tt>y</tt>
 *  \details <tt>n</tt> is the arrays size.
 */
template<class T_>
void Axpy(size_t n,
          T_     a,
          T_*    x,
          T_*    y)
{
   for (size_t i=0; i<n; ++i)
      y[i] += a*x[i];
}


/** \fn void Axpy(T_ a, const vector<T_>& x, vector<T_>& y)
 *  \ingroup Util
 *  \brief Multiply vector <tt>x</tt> by <tt>a</tt> and add result to <tt>y</tt>
 *  \details <tt>x</tt> and <tt>y</tt> are instances of class vector<T_>
 */
template<class T_>
void Axpy(T_                a,
          const vector<T_>& x,
          vector<T_>&       y)
{
   for (size_t i=0; i<x.size(); i++)
      y[i] += a*x[i];
}


/** \fn void Copy(size_t n, T_* x, T_* y)
 *  \ingroup Util
 *  \brief Copy array <tt>x</tt> to <tt>y</tt>
 *  <tt>n</tt> is the arrays size.
 */
template<class T_>
inline void Copy(size_t n,
                 T_*    x,
                 T_*    y)
{
   for (size_t i=0; i<n; ++i)
      y[i] = x[i];
}


/** \fn real_t Error2(const vector<real_t>& x, const vector<real_t>& y)
 *  \ingroup Util
 *  \brief Return absolute L2 error between vectors <tt>x</tt> and <tt>y</tt>
 */
inline real_t Error2(const vector<real_t> &x,
                     const vector<real_t>& y)
{
   real_t s=0;
   for (size_t i=0; i<x.size(); i++)
      s += (x[i]-y[i])*(x[i]-y[i]);
   return fabs(s/x.size());
}


/** \fn real_t RError2(const vector<real_t>& x, const vector<real_t>& y)
 *  \ingroup Util
 *  \brief Return absolute L<sup>2</sup> error between vectors <tt>x</tt> and <tt>y</tt>
 */
inline real_t RError2(const vector<real_t>& x,
                      const vector<real_t>& y)
{
   real_t s=0, t=0;
   for (size_t i=0; i<x.size(); i++) {
      s += (x[i]-y[i])*(x[i]-y[i]);
      t += x[i]*x[i];
   }
   return fabs(s/t);
}


/** \fn real_t ErrorMax(const vector<real_t>& x, const vector<real_t>& y)
 *  \ingroup Util
 *  \brief Return absolute Max. error between vectors <tt>x</tt> and <tt>y</tt>
 */
inline real_t ErrorMax(const vector<real_t>& x,
                       const vector<real_t>& y)
{
   real_t s=0;
   for (size_t i=0; i<x.size(); i++)
      s = std::max(s,fabs(x[i]-y[i]));
   return s;
}


/** \fn real_t RErrorMax(const vector<real_t>& x, const vector<real_t>& y)
 *  \ingroup Util
 *  \brief Return relative Max. error between vectors <tt>x</tt> and <tt>y</tt>
 */
inline real_t RErrorMax(const vector<real_t>& x,
                        const vector<real_t>& y)
{
   real_t s=0, t=0;
   for (size_t i=0; i<x.size(); i++) {
      s = std::max(s,fabs(x[i]-y[i]));
      t = std::max(t,fabs(x[i]));
   }
   return s/t;
}


/** \fn T_ Dot(size_t n, T_* x, T_* y)
 *  \ingroup Util
 *  \brief Return dot product of arrays <tt>x</tt> and <tt>y</tt>
 *  \details <tt>n</tt> is the arrays size.
 */
template<class T_>
inline T_ Dot(size_t n,
              T_*    x,
              T_*    y)
{
   T_ d = T_(0);
   for (size_t i=0; i<n; ++i)
      d += x[i]*Conjg(y[i]);
   return d;
}


/** \fn double Dot(const vector<double>& x, const vector<double>& y)
 *  \ingroup Util
 *  \brief Return dot product of vectors <tt>x</tt> and <tt>y</tt>.
 *  \details <tt>x</tt> and <tt>y</tt> are instances of class vector<double>
 */
inline real_t Dot(const vector<real_t>& x,
                  const vector<real_t>& y)
{
   real_t d = 0;
   for (size_t i=0; i<x.size(); i++)
      d += x[i]*y[i];
   return d;
}


/** \fn real_t operator*(const vector<real_t>& x, const vector<real_t>& y)
 *  \brief Operator * (Dot product of 2 vector instances)
 *  \ingroup VectMat
 *  \return <tt>x.y</tt>
 */
inline real_t operator*(const vector<real_t>& x,
                        const vector<real_t>& y)
{
   return Dot(x,y);
}


/** \fn T_ Dot(const Point<T_>& x, const Point<T_>& y)
 *  \ingroup Util
 *  \brief Return dot product of <tt>x</tt> and <tt>y</tt>
 */
template<class T_>
inline T_ Dot(const Point<T_>& x,
              const Point<T_>& y)
{
   return (x.x*y.x + x.y*y.y + x.z*y.z);
}
    

/** \fn double exprep(double x)
 *  \ingroup Util
 *  \brief Compute the exponential function with avoiding over and underflows
 */
inline real_t exprep(real_t x)
{
   real_t e;
   if (x > 174.)
      e = 3.69e75;
   else if (x < -180.)
      e = 0.;
   else
      e = exp(x);
   return e;
}


/** \fn void Assign(vector<T_>& v, const T_ &a)
 *  \ingroup Util
 *  \brief Assign the value <tt>a</tt> to all entries of a vector <tt>v</tt>
 */
template<class T_>
inline void Assign(vector<T_>& v,
                   const T_&   a)
{
   for (size_t i=0; i<v.size(); i++)
      v[i] = a;
}


/** \fn void clear(vector<T_>& v)
 *  \ingroup Util
 *  \brief Assign <tt>0</tt> to all entries of a vector
 *  @param [in] v Vector to clear
 */
template<class T_>
inline void clear(vector<T_>& v)
{ 
   for (size_t i=0; i<v.size(); i++)
      v[i] = 0;
}


/** \fn void clear(Vect<T_>& v)
 *  \ingroup Util
 *  \brief Assign <tt>0</tt> to all entries of a vector
 *  @param [in] v Vector to clear
 */
template<class T_>
inline void clear(Vect<T_>& v)
{
   for (size_t i=1; i<=v.size(); i++)
      v.set(i,0);
}


/** \fn real_t Nrm2(size_t n, real_t* x)
 *  \ingroup Util
 *  \brief Return 2-norm of array <tt>x</tt>
 *  @param [in] n is Array length
 *  @param [in] x Array to treat
 */
inline real_t Nrm2(size_t  n,
                   real_t* x) { return sqrt(Dot(n,x,x)); }


/** \fn real_t Nrm2(const vector<real_t>& x)
 *  \ingroup Util
 *  \brief Return 2-norm of vector <tt>x</tt>
 */
inline real_t Nrm2(const vector<real_t>& x) { return sqrt(Dot(x,x)); }


/** \fn real_t Nrm2(const Point<T_>& a)
 *  \ingroup Util
 *  \brief Return 2-norm of <tt>a</tt>
 */
template<class T_>
inline real_t Nrm2(const Point<T_>& a) { return sqrt(Dot(a,a)); }


/** \fn bool Equal(real_t x, real_t y, real_t toler=OFELI_EPSMCH)
 *  \ingroup Util
 *  \brief Function to return true if numbers <tt>x</tt> and <tt>y</tt> are close
 *  up to a given tolerance <tt>toler</tt>
 *  \details Default value of tolerance is the constant <tt>OFELI_EPSMCH</tt> 
 */
inline bool Equal(real_t x,
                  real_t y,
                  real_t toler=OFELI_EPSMCH)
{
   if (fabs(x-y) < toler)
      return true;
   else
      return false;
}


/** \fn char itoc(int i)
 *  \ingroup Util
 *  \brief Function to convert an integer to a character.
 */
inline char itoc(int i)
{
   static char buf[10];
   sprintf(buf,"%d",i);
   return buf[0];
}


/** \fn T_ stringTo(const std::string& s)
 *  \ingroup Util
 *  \brief Function to convert a string to a template type parameter.
 */
template <class T_>
T_ stringTo(const std::string& s)
{
   T_ val;
   std::stringstream sstr(s);
   sstr >> val;
   return val;
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \fn void RTrim(char* s)
 *  \ingroup Util
 *  \brief Function to remove blanks at the end of a string.
 */
inline void RTrim(char* s)
{
   size_t i=0, l=string(s).length()-1;
   if (l > 0) {
      i = l;
      while (s[i] == ' ')
         i--;
      if (i < l)
         s[i+1] = 0;
   }
}


/** \fn void LTrim(char* s)
 *  \ingroup Util
 *  \brief Function to remove blanks at the beginning of a string.
 */
inline void LTrim(char* s)
{
   size_t i=0, j=0, l=string(s).length()-1;
   if (l > 0) {
      while (s[i] == ' ')
         i++;
      if (i > 0) {
         for (j=0; j<=l-i; j++)
            s[j] = s[i+j];
         s[j] = 0;
      }
   }
}


/** \fn void Trim(char* s)
 *  \ingroup Util
 *  \brief Function to remove blanks at the beginning and end of a string.
 */
inline void Trim(char* s) { RTrim(s); LTrim(s); }

/** \fn void Swap(T_& a, T_& b)
 *  \ingroup Util
 *  \brief Swap elements \a a and \a b.
 */
template<class T_>
inline void Swap(T_& a,
                 T_& b)
{
   T_ t = a;
   a = b;
   b = t;
}


inline std::string zeros(size_t m,
                         size_t n=3)
{
   int k=1;
   std::string s;
   while (m >= pow(10.,k))
      n--, k++;
   for (size_t l=0; l<n; l++)
      s += "0";
   return s+std::to_string(int(m));
}


static inline std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
   str.erase(0, str.find_first_not_of(chars));
   return str;
}


static inline std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
   str.erase(str.find_last_not_of(chars) + 1);
   return str;
}


static inline std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
   return ltrim(rtrim(str, chars), chars);
}

inline int MaxQuad(const real_t& a, const real_t& b, const real_t& c, real_t& x)
{
   real_t d = b*b - a*c;
   if (d<0.)
      return -1;
   d = sqrt(d);
   x = fmax((-b-d)/a,(-b+d)/a);
   return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
