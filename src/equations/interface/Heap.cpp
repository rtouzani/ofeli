/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

Heap::Heap(size_t size)
     : _heap(size), _max_size(0), _size(0)
{
}


Heap::~Heap()
{
}


Heap::Heap(Vect<IPoint>& vec)
     : _max_size(vec.size()*vec.size()), _size(0)
{
   _heap.resize(_max_size);
   for (size_t i=0; i<vec.size(); i++, _size++)
      Add(vec[i]);
}


Heap::Heap(Heap& H)
     : _max_size(H._max_size), _size(H._size)
{
   _heap.resize(_max_size);
   for (size_t i=0; i<_size; ++i)
      _heap[i] = H._heap[i];
}


Heap& Heap::operator=(const Heap& H)
{
   _size = H._size;
   _heap.resize(_size);
   for (size_t i=0; i<_size; ++i)
      _heap[i] = H._heap[i];
   return (*this);
}


IPoint Heap::operator[](size_t ind) const
{
   if (ind > _size)
      throw OFELIException("Heap::operator[]: Out of bounds");
   return _heap[ind];
}


IPoint& Heap::operator[](size_t ind)
{
   if (ind > _size)
      throw OFELIException("Heap::operator[]: Out of bounds");
   return _heap[ind];
}


IPoint Heap::Current()
{
   IPoint el;
   if (_size != 0) {
      el = _heap(1);
      _heap(1) = _heap(_size);
      _size--;
      Down_Heap();
   }
   return el;
}


void Heap::Add(IPoint el)
{
   _heap(++_size) = el;
std::cout<<"             Added: "<<el.getX()<<"  "<<el.getY()<<endl;
   Up_Heap(_size);
}


void Heap::Down_Heap(size_t rank)
{
    while (rank <= _size/2) {
      size_t rgFils = 2*rank;  // son on the left rank
      if (rgFils+1 <= _size) { // son on the right exist
         if (_heap(rgFils) > _heap(rgFils+1))
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
   while (rank>1 && _heap(rank/2)>_heap(rank)) {
      IPoint aux = _heap(rank/2);
      _heap(rank/2) = _heap(rank);
      _heap(rank)  = aux;
      rank /= 2;
   }
}


bool Heap::Find(const IPoint& pt,
                size_t&       ind)
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
                  size_t        rg)
{
   _heap[rg].setValue(val);
   Up_Heap(rg);
}


ostream & operator<<(ostream&    s,
                     const Heap& H)
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
