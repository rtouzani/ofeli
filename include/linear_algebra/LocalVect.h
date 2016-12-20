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

          Definition of Template Class LocalVect for small size vectors

  ==============================================================================*/


#ifndef __LOCAL_VECT_H
#define __LOCAL_VECT_H

#include "OFELI_Config.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "mesh/Node.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file LocalVect.h
 *  \brief Definition file for class LocalVect.
 */

/*! \class LocalVect
 *  \ingroup VectMat
 *  \brief Handles small size vectors like element vectors.
 *
 *  \details The template class LocalVect treats small size vectors. Typically, this class is
 *  recommended to store element and side arrays.
 *  Operators \b =, \b [] and \b () are overloaded so that one can write for instance:
 *
 * \verbatim
        LocalVect<double,10> u, v;
        v = -1.0;
        u = v;
        u(3) = -2.0;
   \endverbatim

 * to set vector \b v entries to \b -1, copy vector \b v into vector \b u and
 * assign third entry of \b v to \b -2. Notice that entries of \b v are here
 * \b v(1), \b v(2), ..., \b v(10), \e i.e. vector entries start at index \b 1.\n
 * Internally, no dynamic storage is used.
 *
 * \tparam T_ Data type (double, float, complex<double>, ...)
 * \tparam N_ Vector size
 */

template<class T_> class Vect;
class Element;
class Side;


template<class T_,size_t N_>
class LocalVect
{

public:

/// \brief Default constructor.
    LocalVect();

/// \brief Constructor using a C-array
    LocalVect(const T_ *a);

/// \brief Constructor using Element pointer
    LocalVect(const Element *el);

/// \brief Constructor using Side pointer
    LocalVect(const Side *sd);

/// \brief Copy constructor
    LocalVect(const LocalVect<T_,N_>& v);

/** \brief Constructor of an element vector from a global Vect instance.
 *  \details The constructed vector has local numbering of nodes
 *  @param [in] el Pointer to Element to localize
 *  @param [in] v Global vector to localize
 *  @param [in] opt Option for DOF treatment
 *     - = 0, Normal case [Default]
 *     - Any other value : only one DOF is handled (Local vector has as dimension number of degrees of freedom)
 */
    LocalVect(const Element*  el,
              const Vect<T_>& v,
                    int       opt=0);

/** \brief Constructor of a side vector from a global Vect instance.
 *  \details The constructed vector has local numbering of nodes
 *  @param [in] sd Pointer to Side to localize
 *  @param [in] v Global vector to localize
 *  @param [in] opt Option for DOF treatment
 *     - = 0, Normal case [Default]
 *     - Any other value : only one DOF is handled (Local vector has as dimension number of degrees of freedom)
 */
    LocalVect(const Side*     sd,
              const Vect<T_>& v,
                    int       opt=0);

/// \brief Destructor
    ~LocalVect() { }

/** \brief Localize an element vector from a global Vect instance.
 *  \details The constructed vector has local numbering of nodes
 *  This function is called by the constructor:
 *  <tt>LocalVect(const Element *el, const Vect<T_> &v)</tt>
 *  @param [in] el Pointer to Element to localize
 *  @param [in] v Global vector to localize
 *  @param [in] type Type of element. This is to be chosen
 *  among enumerated values:
 *  <tt>LINE2</tt>, <tt>TRIANG3</tt>, <tt>QUAD4</tt>, <tt>TETRA4</tt>, <tt>HEXA8</tt>,
 *  <tt>PENTA6</tt>
 */
    void getLocal(const Element&  el,
                  const Vect<T_>& v,
                        int       type);

/** \brief Localize an element vector from a global Vect instance.
 *  \details The constructed vector has local numbering of nodes
 *  This function is called by the constructor:
 *  \b LocalVect(const Element *el, const Vect<T_> &v)
 *  @param [in] el Pointer to Side to localize
 *  @param [in] v Global vector to localize
 *  @param [in] k Degree of freedom to localize [Default: All degrees of freedom are stored]
 */
    void Localize(const Element*  el,
                  const Vect<T_>& v,
                        size_t    k=0);

/** \brief Localize a side vector from a global Vect instance.
 *  \details The constructed vector has local numbering of nodes
 *  This function is called by the constructor:
 *  \b LocalVect(const Side *sd, const Vect<T_> &v)
 *  @param [in] sd Pointer to Side to localize
 *  @param [in] v Global vector to localize
 *  @param [in] k Degree of freedom to localize [Default: All degrees of freedom are stored]
 */
    void Localize(const Side*     sd,
                  const Vect<T_>& v,
                        size_t    k=0);

/// \brief Operator <tt>[]</tt> (Non constant version).
/// \details <tt>v[i]</tt> starts at <tt>v[0]</tt> to <tt>v[size()-1]</tt>
    T_ &operator[](size_t i) { return _v[i]; }

/// \brief Operator <tt>[]</tt> (Constant version).
/// \details <tt>v[i]</tt> starts at <tt>v[0]</tt> to <tt>v[size()-1]</tt>
    T_ operator[](size_t i) const { return _v[i]; }

/** \brief Operator <tt>()</tt> (Non constant version).
 *  \details <tt>v(i)</tt> starts at <tt>v(1)</tt> to <tt>v(size())</tt>.
 *  <tt>v(i)</tt> is the same element as <tt>v[i-1]</tt>
 */
    T_ &operator()(size_t i) { return _v[i-1]; }

/** \brief Operator <tt>()</tt> (Constant version).
 *  \details <tt>v(i)</tt> starts at <tt>v(1)</tt> to <tt>v(size())</tt>
 *  <tt>v(i)</tt> is the same element as <tt>v[i-1]</tt>
 */
    T_ operator()(size_t i) const { return _v[i-1]; }

/// \brief Return pointer to Element if vector was constructed
/// using an element and <tt>NULL</tt> otherwise
    Element *El() { return _el; }

/// \brief Return pointer to Side if vector was constructed
/// using a side and <tt>NULL</tt> otherwise
    Side *Sd() { return _sd; }

/// \brief Operator <tt>=</tt>
/// \details Copy a LocalVect instance to the current one
    LocalVect<T_,N_> & operator=(const LocalVect<T_,N_>& v);

/// \brief Operator <tt>=</tt>
/// \details Assign value <tt>x</tt> to all vector entries
    LocalVect<T_,N_> & operator=(const T_& x);

/// \brief Operator <tt>+=</tt>
/// \details Add vector <tt>v</tt> to this instance
    LocalVect<T_,N_> & operator+=(const LocalVect<T_,N_>& v);

/// \brief Operator <tt>+=</tt>
/// \details Add constant <tt>a</tt> to vector entries
    LocalVect<T_,N_> & operator+=(const T_& a);

/// \brief Operator <tt>-=</tt>
/// \details Subtract vector <tt>v</tt> from this instance
    LocalVect<T_,N_> & operator-=(const LocalVect<T_,N_>& v);

/// \brief Operator <tt>-=</tt>
/// \details Subtract constant <tt>a</tt> from vector entries
    LocalVect<T_,N_> & operator-=(const T_& a);

/// \brief Operator <tt>*=</tt>
/// \details Multiply vector by constant <tt>a</tt>
    LocalVect<T_,N_> & operator*=(const T_& a);

/// \brief Operator <tt>/=</tt>
/// \details Divide vector by constant <tt>a</tt>
    LocalVect<T_,N_> & operator/=(const T_& a);

/// \brief Return pointer to vector as a C-Array
    T_ *get() { return _v; }
    
/** \brief Return Dot (scalar) product of two vectors
 *  \details A typical use of this operator is <tt>double a = (v,w)</tt>
 *  where <tt>v</tt> and <tt>w</tt> are 2 instances of <tt>LocalVect<double,n></tt>
 *  @param [in] v LocalVect instance by which the current instance is multiplied
 */
    T_ operator,(const LocalVect<T_,N_>& v) const;

 private:

   T_      _v[N_];
   Element *_el;
   Side    *_sd;
};

///////////////////////////////////////////////////////////////////////////////
//                         I M P L E M E N T A T I O N                       //
///////////////////////////////////////////////////////////////////////////////


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect()
{
   _el = NULL;
   _sd = NULL;
   for (size_t i=0; i<N_; ++i)
      _v[i] = 0;
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const T_* a)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] = a[i];
   _el = NULL;
   _sd = NULL;
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const Element* el)
{
   _el = el;
   _sd = NULL;
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const Side* sd)
{
   _el = NULL;
   _sd = sd;
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const LocalVect<T_,N_>& v)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] = v[i];
   _el = v._el;
   _sd = v._sd;
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const Element*  el,
                            const Vect<T_>& v,
                                  int       opt)
{
   if (opt==0)
      Localize(el,v);
   else {
      size_t i = 0;
      for (size_t n=1; n<=el->getNbNodes(); ++n)
         _v[i++] = v((*el)(n)->n());
   }
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const Side*     sd,
                            const Vect<T_>& v,
                                  int       opt)
{
   if (opt==0)
      Localize(sd,v);
   else {
      size_t i = 0;
      for (size_t n=1; n<=sd->getNbNodes(); ++n)
         _v[i++] = v(sd->getNodeLabel(n));
   }
}


template<class T_,size_t N_>
void LocalVect<T_,N_>::getLocal(const Element&  el,
                                const Vect<T_>& v,
                                      int       type)
{
   size_t i=0;
   if (type==LINE2) {
      for (size_t n=1; n<=2; ++n) {
         Node *nd = el(n);
         size_t k = nd->getFirstDOF();
         for (size_t j=1; j<=nd->getNbDOF(); ++j)
            _v[i++] = v(k++);
      }
   }
   else if (type==TRIANG3) {
      for (size_t n=1; n<=3; ++n) {
         Node *nd = el(n);
         size_t k = nd->getFirstDOF();
         for (size_t j=1; j<=nd->getNbDOF(); ++j)
            _v[i++] = v(k++);
      }
   }
   else if (type==QUAD4) {
      for (size_t n=1; n<=4; ++n) {
         Node *nd = el(n);
         size_t k = nd->getFirstDOF();
         for (size_t j=1; j<=nd->getNbDOF(); ++j)
            _v[i++] = v(k++);
      }
   }
   else if (type==TETRA4) {
      for (size_t n=1; n<=4; ++n) {
         Node *nd = el(n);
         size_t k = nd->getFirstDOF();
         for (size_t j=1; j<=nd->getNbDOF(); ++j)
            _v[i++] = v(k++);
      }
   }
   else if (type==HEXA8) {
      for (size_t n=1; n<=8; ++n) {
         Node *nd = el(n);
         size_t k = nd->getFirstDOF();
         for (size_t j=1; j<=nd->getNbDOF(); ++j)
            _v[i++] = v(k++);
      }
   }
   else if (type==PENTA6) {
      for (size_t n=1; n<=6; ++n) {
         Node *nd = el(n);
         size_t k = nd->getFirstDOF();
         for (size_t j=1; j<=nd->getNbDOF(); ++j)
            _v[i++] = v(k++);
      }
   }
}


template<class T_,size_t N_>
void LocalVect<T_,N_>::Localize(const Element*  el,
                                const Vect<T_>& v,
                                      size_t    k)
{
   size_t i = 0;
   if (k==0) {
      for (size_t n=1; n<=el->getNbNodes(); ++n) {
         Node *nd = (*el)(n);
         for (size_t j=1; j<=nd->getNbDOF(); ++j)
            _v[i++] = v(nd->getDOF(j));
      }
   }
   else {
      for (size_t n=1; n<=el->getNbNodes(); ++n)
         _v[i++] = v((*el)(n)->getDOF(k));
   }
}


template<class T_,size_t N_>
void LocalVect<T_,N_>::Localize(const Side*     sd,
                                const Vect<T_>& v,
                                      size_t    k)
{
   size_t i = 0;
   if (k==0) {
      for (size_t n=1; n<=sd->getNbNodes(); ++n) {
         for (size_t j=1; j<=sd->getNbDOF(); ++j)
            _v[i++] = v(sd->getDOF(j));
      }
   }
   else {
      for (size_t n=1; n<=sd->getNbNodes(); ++n)
         _v[i++] = v((*sd)(n)->getDOF(k));
   }
}


template<class T_,size_t N_>
LocalVect<T_,N_> & LocalVect<T_,N_>::operator=(const LocalVect<T_,N_>& v)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] = v._v[i];
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_> & LocalVect<T_,N_>::operator=(const T_& x)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] = x;
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_> & LocalVect<T_,N_>::operator+=(const LocalVect<T_,N_>& v)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] += v[i];
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_> & LocalVect<T_,N_>::operator+=(const T_& a)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] += a;
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_> & LocalVect<T_,N_>::operator-=(const LocalVect<T_,N_>& v)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] -= v[i];
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_> & LocalVect<T_,N_>::operator-=(const T_& a)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] -= a;
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_> & LocalVect<T_,N_>::operator*=(const T_& a)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] *= a;
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_> & LocalVect<T_,N_>::operator/=(const T_& a)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] /= a;
   return *this;
}
    
    
template<class T_,size_t N_>
T_ LocalVect<T_,N_>::operator,(const LocalVect<T_,N_>& v) const
{
   T_ p = 0;
   for (size_t i=0; i<N_; i++)
      p += _v[i] * v[i];
   return p;
}

	
///////////////////////////////////////////////////////////////////////////////
//                A S S O C I A T E D    F U N C T I O N S                   //
///////////////////////////////////////////////////////////////////////////////

/** \fn LocalVect<T_,N_> operator+(const LocalVect<T_,N_>& x, const LocalVect<T_,N_>& y)
 *  \brief Operator + (Add two vectors)
 *  \ingroup VectMat
 *  \return    x+y
 */
template<class T_,size_t N_>
LocalVect<T_,N_> operator+(const LocalVect<T_,N_>& x,
                           const LocalVect<T_,N_>& y)
{
   LocalVect<T_,N_> z(x);
   for (size_t i=0; i<N_; ++i)
      z[i] += y[i];
   return z;
}

/** \fn LocalVect<T_,N_> operator-(const LocalVect<T_,N_> &x, const LocalVect<T_,N_> &y)
 *  \brief Operator - (Subtract two vectors)
 *  \ingroup VectMat
 *  \return    x-y
 */
template<class T_,size_t N_>
LocalVect<T_,N_> operator-(const LocalVect<T_,N_>& x,
                           const LocalVect<T_,N_>& y)
{
   LocalVect<T_,N_> z(x);
   for (size_t i=0; i<N_; ++i)
      z[i] -= y[i];
   return z;
}


/** \fn LocalVect<T_,N_> operator*(T_ a, const LocalVect<T_,N_> &x)
 *  \brief Operator * (Premultiplication of vector by constant)
 *  \ingroup VectMat
 *  \return    a*x
 */
template<class T_,size_t N_>
LocalVect<T_,N_> operator*(      T_                a,
                           const LocalVect<T_,N_>& x)
{
   LocalVect<T_,N_> v(x);
   for (size_t i=0; i<N_; ++i)
      v[i] *= a;
   return v;
}


/** \fn LocalVect<T_,N_> operator/(T_ a, const LocalVect<T_,N_> &x)
 *  \brief Operator / (Division of vector by constant)
 *  \ingroup VectMat
 *  \return    x/a
 */
template<class T_,size_t N_>
LocalVect<T_,N_> operator/(      T_                a,
                           const LocalVect<T_,N_>& x)
{
   LocalVect<T_,N_> v(x);
   for (size_t i=0; i<N_; ++i)
      v[i] /= a;
   return v;
}


///////////////////////////////////////////////////////////////////////////////
//                       T H E   M E T A P R O G R A M S                     //
///////////////////////////////////////////////////////////////////////////////

// Dot Product

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<size_t I_>
struct _meta_dot
{
  template<class T_, size_t N_>
  static T_ _dot(const LocalVect<T_,N_> &a, const LocalVect<T_,N_> &b)
  { return a[I_]*b[I_] + _meta_dot<I_-1>::_dot(a,b); }
};


template<>
struct _meta_dot<0>
{
  template<class T_, size_t N_>
  static T_ _dot(const LocalVect<T_,N_>& a,
                 const LocalVect<T_,N_>& b)
  { return a[0]*b[0]; }
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \fn double Dot(const LocalVect<T_,N_> &a, const LocalVect<T_,N_> &b)
 *  \brief Calculate dot product of 2 vectors (instances of class LocalVect)
 *  \ingroup VectMat
 *  \return  Dot product
 */
template<class T_, size_t N_>
inline real_t Dot(const LocalVect<T_,N_>& a,
                  const LocalVect<T_,N_>& b)
{
   return _meta_dot<N_-1>::_dot(a,b);
}


// scale

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<size_t I_>
struct _meta_scale
{
   template<class T_, size_t N_>
   static void _scale(      T_                a,
                      const LocalVect<T_,N_>& x,
                            LocalVect<T_,N_>& y)
   {
      y[I_] = a*x[I_];
      _meta_scale<I_-1>::_scale(a,x,y);
   }
};


template<>
struct _meta_scale<0>
{
   template<class T_, size_t N_>
   static void _scale(      T_                a,
                      const LocalVect<T_,N_>& x,
                            LocalVect<T_,N_>& y)
   {
      y[0] = a*x[0];
   }
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \fn void Scale(T_ a, const LocalVect<T_,N_> &x, LocalVect<T_,N_> &y)
 *  \ingroup VectMat
 *  \brief Multiply vector <tt>x</tt> by constant <tt>a</tt> and store result in <tt>y</tt>.
 */
template<class T_, size_t N_>
inline void Scale(      T_                a,
                  const LocalVect<T_,N_>& x,
                        LocalVect<T_,N_>& y)
{
   _meta_scale<N_-1>::_scale(a,x,y);
}


// scale (in place)

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<size_t I_>
struct _meta_scale1
{
  template<class T_, size_t N_>
  static void _scale(T_                a,
                     LocalVect<T_,N_>& x)
  { x[I_] *= a; _meta_scale1<I_-1>::_scale(a,x); }
};

template<>
struct _meta_scale1<0>
{
  template<class T_, size_t N_>
  static void scale_(T_                a,
                     LocalVect<T_,N_>& x)
  { x[0] *= a; }
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \fn void Scale(T_ a, LocalVect<T_,N_> &x)
/// \ingroup VectMat
/// \brief Multiply vector <tt>x</tt> by constant <tt>a</tt> and store result in <tt>x</tt>.
template<class T_, size_t N_>
inline void Scale(T_                a,
                  LocalVect<T_,N_>& x)
{
   _meta_scale1<N_-1>::_scale(a,x);
}

// axpy

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<size_t I_>
struct _meta_axpy
{
  template<class T_, size_t N_>
  static void _axpy(      T_                a,
                    const LocalVect<T_,N_>& x,
                          LocalVect<T_,N_>& y)
  {
     y[I_] += a*x[I_];
     _meta_axpy<I_-1>::_axpy(a,x,y);
  }
};


template<>
struct _meta_axpy<0>
{
  template<class T_, size_t N_>
  static void _axpy(      T_                a,
                    const LocalVect<T_,N_>& x,
                          LocalVect<T_,N_>& y)
  {
     y[0] += a*x[0];
  }
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \fn Axpy(T_ a, const LocalVect<T_,N_> &x, LocalVect<T_,N_> &y)
 *  \ingroup VectMat
 * \brief Add <tt>a*x</tt> to vector <tt>y</tt>.
 */
template<class T_, size_t N_>
inline void Axpy(      T_                a,
                 const LocalVect<T_,N_>& x,
                       LocalVect<T_,N_>& y)
{
   _meta_axpy<N_-1>::_axpy(a,x,y);
}


// xpy

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<size_t I_>
struct _meta_xpy
{
  template<class T_, size_t N_>
  void xpy_(const LocalVect<T_,N_>& x,
                  LocalVect<T_,N_>& y)
  {
     y[I_] += x[I_];
     _meta_axpy<I_-1>::_xpy(x,y);
  }
};


template<>
struct _meta_xpy<0>
{
  template<class T_, size_t N_>
  void _xpy(const LocalVect<T_,N_>& x,
                  LocalVect<T_,N_>& y)
  {
     y[0] += x[0];
  }
};


/** \fn Xpy(const LocalVect<T_,N_> &x, LocalVect<T_,N_> &y)
 *  \ingroup VectMat
 *  Add <tt>x</tt> to vector <tt>y</tt>.
 */
template<class T_, size_t N_>
inline void Xpy(const LocalVect<T_,N_>& x,
                      LocalVect<T_,N_>& y)
{
   _meta_xpy<N_-1>::_xpy(x,y);
}

// copy

template<size_t I_>
struct _meta_copy
{
  template<class T_, size_t N_>
  void _copy(const LocalVect<T_,N_>& x,
                   LocalVect<T_,N_>& y)
  {
     y[I_] = x[I_];
     _meta_axpy<I_-1>::_xpy(x,y);
  }
};


template<>
struct _meta_copy<0>
{
  template<class T_, size_t N_>
  void _copy(const LocalVect<T_,N_>& x,
                   LocalVect<T_,N_>& y)
  { y[0] = x[0]; }
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \fn void Copy(const LocalVect<T_,N_> &x, LocalVect<T_,N_> &y)
 *  \ingroup VectMat
 *  \brief Copy vector \a x into vector \a y.
 */
template<class T_, size_t N_>
inline void Copy(const LocalVect<T_,N_>& x,
                       LocalVect<T_,N_>& y)
{
   _meta_copy<N_-1>::_copy(x,y);
}


/** \fn ostream& operator<<(ostream& s, const LocalVect<T_,N_> & v)
 *  \brief Output vector in output stream
 *  \ingroup VectMat
 */
template<class T_,size_t N_>
ostream& operator<<(      ostream&          s,
                    const LocalVect<T_,N_>& v)
{
   s.setf(ios::scientific);
   for (size_t i=1; i<=N_; i++)
      s << setw(6) << i << "  " << std::setprecision(8) << std::setw(18) << v(i) << std::endl;
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
