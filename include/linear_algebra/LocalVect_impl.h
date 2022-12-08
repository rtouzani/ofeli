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

                        Implementation of class LocalVect

  ==============================================================================*/


#ifndef __LOCAL_VECT_IMPL_H
#define __LOCAL_VECT_IMPL_H

#include "linear_algebra/LocalVect.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "mesh/Node.h"
#include "OFELIException.h"

namespace OFELI {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect()
                 : _el(nullptr), _sd(nullptr)
{
   _v.resize(N_);
   for (size_t i=0; i<N_; ++i)
      _v[i] = 0;
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const T_ *a)
                 : _el(nullptr), _sd(nullptr)
{
   _v.resize(N_);
   for (size_t i=0; i<N_; ++i)
      _v[i] = a[i];
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const Element *el)
                 : _el(el), _sd(nullptr)
{
   _v.resize(N_);

}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const Side *sd)
                 : _el(nullptr), _sd(sd)
{
   _v.resize(N_);
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const LocalVect<T_,N_>& v)
{
   _v.resize(N_);
   for (size_t i=0; i<N_; ++i)
      _v[i] = v[i];
   _el = v._el;
   _sd = v._sd;
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const Element*  el,
                            const Vect<T_>& v,
                            int             opt)
{
   _v.resize(N_);
   if (opt==0)
      Localize(el,v);
   else {
      size_t i = 0;
      for (size_t n=1; n<=el->getNbNodes(); ++n)
         _v[i++] = v((*el)(n)->n());
   }
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const Element&  el,
                            const Vect<T_>& v,
                            int             opt)
{
   LocalVect<T_,N_>(&el,v,opt);
}


template<class T_,size_t N_>
LocalVect<T_,N_>::LocalVect(const Side*     sd,
                            const Vect<T_>& v,
                            int             opt)
{
   _v.resize(N_);
   if (opt==0)
      Localize(sd,v);
   else {
      size_t i = 0;
      for (size_t n=1; n<=sd->getNbNodes(); ++n)
         _v[i++] = v(sd->getNodeLabel(n));
   }
}


template<class T_,size_t N_>
LocalVect<T_,N_>::~LocalVect() { }


template<class T_,size_t N_>
void LocalVect<T_,N_>::getLocal(const Element&  el,
                                const Vect<T_>& v,
                                int             type)
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
                                size_t          k)
{
   size_t i=0;
   if (k==0) {
      for (size_t n=1; n<=el->getNbNodes(); ++n) {
         Node *nd = (*el)(n);
         for (size_t j=1; j<=nd->getNbDOF(); ++j)
            _v[i++] = v(nd->getFirstDOF()+j-1);
      }
   }
   else
      for (size_t n=1; n<=el->getNbNodes(); ++n)
         _v[i++] = v((*el)(n)->getFirstDOF()+k-1);
}


template<class T_,size_t N_>
void LocalVect<T_,N_>::Localize(const Side*     sd,
                                const Vect<T_>& v,
                                size_t          k)
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
T_& LocalVect<T_,N_>::operator[](size_t i) { return _v[i]; }


template<class T_,size_t N_>
T_ LocalVect<T_,N_>::operator[](size_t i) const { return _v[i]; }


template<class T_,size_t N_>
T_& LocalVect<T_,N_>::operator()(size_t i) { return _v[i-1]; }


template<class T_,size_t N_>
T_ LocalVect<T_,N_>::operator()(size_t i) const { return _v[i-1]; }


template<class T_,size_t N_>
Element* LocalVect<T_,N_>::El() { return _el; }


template<class T_,size_t N_>
Side* LocalVect<T_,N_>::Sd() { return _sd; }


template<class T_,size_t N_>
LocalVect<T_,N_>& LocalVect<T_,N_>::operator=(const LocalVect<T_,N_>& v)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] = v._v[i];
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_>& LocalVect<T_,N_>::operator=(const T_& x)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] = x;
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_>& LocalVect<T_,N_>::operator+=(const LocalVect<T_,N_>& v)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] += v[i];
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_>& LocalVect<T_,N_>::operator+=(const T_& a)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] += a;
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_>& LocalVect<T_,N_>::operator-=(const LocalVect<T_,N_>& v)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] -= v[i];
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_>& LocalVect<T_,N_>::operator-=(const T_& a)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] -= a;
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_>& LocalVect<T_,N_>::operator*=(const T_& a)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] *= a;
   return *this;
}


template<class T_,size_t N_>
LocalVect<T_,N_>& LocalVect<T_,N_>::operator/=(const T_& a)
{
   for (size_t i=0; i<N_; ++i)
      _v[i] /= a;
   return *this;
}


template<class T_,size_t N_>
T_* LocalVect<T_,N_>::get() { return &_v[0]; }


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


template<class T_,size_t N_>
LocalVect<T_,N_> operator+(const LocalVect<T_,N_>& x,
                           const LocalVect<T_,N_>& y)
{
   LocalVect<T_,N_> z(x);
   for (size_t i=0; i<N_; ++i)
      z[i] += y[i];
   return z;
}


template<class T_,size_t N_>
LocalVect<T_,N_> operator-(const LocalVect<T_,N_>& x,
                           const LocalVect<T_,N_>& y)
{
   LocalVect<T_,N_> z(x);
   for (size_t i=0; i<N_; ++i)
      z[i] -= y[i];
   return z;
}


template<class T_,size_t N_>
LocalVect<T_,N_> operator*(T_                      a,
                           const LocalVect<T_,N_>& x)
{
   LocalVect<T_,N_> v(x);
   for (size_t i=0; i<N_; ++i)
      v[i] *= a;
   return v;
}


template<class T_,size_t N_>
LocalVect<T_,N_> operator/(T_                      a,
                           const LocalVect<T_,N_>& x)
{
   LocalVect<T_,N_> v(x);
   for (size_t i=0; i<N_; ++i)
      v[i] /= a;
   return v;
}


template<class T_, size_t N_>
void Scale(T_                      a,
           const LocalVect<T_,N_>& x,
           LocalVect<T_,N_>&       y)
{
   _meta_scale<N_-1>::_scale(a,x,y);
}


template<class T_, size_t N_>
void Scale(T_                a,
           LocalVect<T_,N_>& x)
{
   _meta_scale1<N_-1>::_scale(a,x);
}


template<class T_, size_t N_>
void Axpy(T_                      a,
          const LocalVect<T_,N_>& x,
          LocalVect<T_,N_>&       y)
{
   _meta_axpy<N_-1>::_axpy(a,x,y);
}


template<class T_, size_t N_>
void Xpy(const LocalVect<T_,N_>& x,
         LocalVect<T_,N_>&       y)
{
   _meta_xpy<N_-1>::_xpy(x,y);
}


template<class T_, size_t N_>
void Copy(const LocalVect<T_,N_>& x,
          LocalVect<T_,N_>&       y)
{
   _meta_copy<N_-1>::_copy(x,y);
}


template<class T_,size_t N_>
ostream& operator<<(ostream&                s,
                    const LocalVect<T_,N_>& v)
{
   s.setf(ios::scientific);
   for (size_t i=1; i<=N_; i++)
      s << setw(6) << i << "  " << std::setprecision(8) << std::setw(18) << v(i) << std::endl;
   return s;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
