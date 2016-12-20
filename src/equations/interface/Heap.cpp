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

                           Implementation of class 'Heap'

  ==============================================================================*/

#include "equations/interface/Heap.h"

namespace OFELI {

Heap::Heap(size_t size):_heap(size)
{
   _size = 0;
   _max_size = _size;
}


Heap::~Heap()
{
}


Heap::Heap(Vect<IPoint>& vec)
{
   _size = 0;
   _max_size = vec.size()*vec.size();
   _heap.resize(_max_size);
   for (size_t i=0; i<vec.size(); ++i) {
      _size++;
      Add(vec[i]);
   }
}


Heap::Heap(Heap& H)
{
   _max_size = H._max_size;
   _size = H._size;
   _heap.resize(_max_size);
   for (size_t i=0; i<_size; ++i)
      _heap[i] = H._heap[i];
}


Heap & Heap::operator=(const Heap& H)
{
   _size = H._size;
   _heap.resize(_size);
   for (size_t i=0; i<_size; ++i)
      _heap[i] = H._heap[i];
   return (*this);
}


IPoint Heap::operator[](size_t ind) const
{
   if (ind>_size) {
      cerr << "ERROR: OutOfBound\n";
      exit(-1);
   }
   return _heap[ind];
}


IPoint & Heap::operator[](size_t ind)
{
   if (ind>_size) {
      cerr << "ERROR: OutOfBound\n";
      exit(-1);
   }
   return _heap[ind];
}


IPoint Heap::Service()
{
   IPoint elem;
   if (_size != 0) {
      elem = _heap(1);
      _heap(1) = _heap(_size);
      _size--;
      Down_Heap();
   }
   return elem;
}


void Heap::Add(IPoint el)
{
   _size++;
   _heap(_size) = el;
   Up_Heap(_size);
}


void Heap::Down_Heap(size_t rank)
{
   while (rank <= _size/2) {
      size_t rgFils = 2*rank;  // son's on the left rank
      if (rgFils+1 <= _size) { // son on the right exist
         if (_heap(rgFils) > _heap(rgFils + 1) )
            rgFils++;
      }
      if (_heap(rank) > _heap(rgFils)) {
         IPoint aux = _heap(rank);
         _heap(rank) = _heap(rgFils);
         _heap(rgFils) = aux;
      }
      rank = rgFils;
   }
}


void Heap::Up_Heap(size_t rank)
{
   while (rank > 1 && (_heap(rank/2) > _heap(rank))) {
      IPoint aux = _heap(rank/2);
      _heap(rank/2) = _heap(rank);
      _heap(rank)  = aux;
      rank = rank/2;
   }
}


bool Heap::Find(const IPoint& pt,
                      size_t& ind)
{
   for (size_t i=0; i<_size; ++i) {
      if (_heap[i] == pt) {
         ind = i;
         return true;
      }
   }
   return false;
}


void Heap::Update(const real_t& val,
                        size_t  rg)
{
   _heap[rg].setValue(val);
   Up_Heap(rg);
}


ostream & operator<<(      ostream& s,
                     const Heap&    H)
{
   s << " The heap: [ \n";
   for (size_t i=0; i<H._size; ++i) {
      s << "\t(" << H[i].getX() << "," << H[i].getY();
      if (H[i].getZ() != 0)
         s << "," << H[i].getZ();
      s << "," << H[i].getValue() << ") \n";
   }
   s << " ] \n";
   return s;
}

} /* namespace OFELI */
