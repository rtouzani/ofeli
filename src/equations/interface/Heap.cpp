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

                         Implementation of class 'Heap'

  ==============================================================================*/

#include "equations/interface/Heap.h"

namespace OFELI {


Heap::Heap()
     : _max_size(0), _size(0)
{
}


Heap::Heap(size_t size)
     : _heap(size), _max_size(0), _size(0)
{
}


Heap::~Heap()
{
}


Heap::Heap(Vect<IPoint>& v)
     : _max_size(v.size()*v.size()), _size(0)
{
   _heap.resize(_max_size);
   for (size_t i=0; i<v.size(); i++, _size++)
      Add(v[i]);
}


Heap::Heap(Heap& h)
     : _max_size(h._max_size), _size(h._size)
{
   _heap.resize(_max_size);
   for (size_t i=0; i<_size; ++i)
      _heap[i] = h._heap[i];
}


void Heap::set(size_t size)
{
   _size = 0;
   _max_size = size;
   _heap.resize(size);
}


Heap& Heap::operator=(const Heap& h)
{
   _size = h._size;
   _heap.resize(_size);
   for (size_t i=0; i<_size; ++i)
      _heap[i] = h._heap[i];
   return (*this);
}


IPoint Heap::operator[](size_t i) const
{
   if (i>_size-1)
      throw OFELIException("Heap::operator[]: Out of bounds");
   return _heap[i];
}


IPoint& Heap::operator[](size_t i)
{
   if (i>_size-1)
      throw OFELIException("Heap::operator[]: Out of bounds");
   return _heap[i];
}


IPoint Heap::Current()
{
   IPoint el;
   if (_size != 0) {
      el = _heap[0];
      _heap[0] = _heap[_size-1];
      _size--;
      Down();
   }
   return el;
}


void Heap::Add(IPoint el)
{
   _heap[_size++] = el;
   Up(_size);
}


void Heap::Down(size_t rank)
{
    while (rank <= _size/2) {
      size_t rgFils = 2*rank;  // son on the left rank
      if (rgFils+1<=_size) {   // son on the right exist
         if (_heap[rgFils-1]>_heap[rgFils])
            rgFils++;
      }
      if (_heap[rank-1]>_heap[rgFils-1]) {
         IPoint aux = _heap[rank-1];
         _heap[rank-1] = _heap[rgFils-1];
         _heap[rgFils-1] = aux;
      }
      rank = rgFils;
   }
}


void Heap::Up(size_t rank)
{
   while (rank>1 && _heap[rank/2-1]>_heap[rank-1]) {
      IPoint aux = _heap[rank/2-1];
      _heap[rank/2-1] = _heap[rank-1];
      _heap[rank-1] = aux;
      rank /= 2;
   }
}


bool Heap::Find(const IPoint& pt,
                size_t&       i)
{
   for (size_t j=0; j<_size; ++j) {
      if (_heap[j] == pt) {
         i = j;
         return true;
      }
   }
   return false;
}


void Heap::Update(const real_t& val,
                  size_t        rg)
{
   _heap[rg].val = val;
   Up(rg);
}


ostream & operator<<(ostream&    s,
                     const Heap& h)
{
   s << "[ \n";
   for (size_t i=0; i<h._size; ++i) {
      s << "\t(" << h[i].i << "," << h[i].j;
      if (h[i].k != 0)
         s << "," << h[i].k;
      s << "," << h[i].val << ") \n";
   }
   s << " ]" << endl;
   return s;
}

} /* namespace OFELI */
