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

      Definition and implementation of class 'Heap' to manage binary heaps

  ==============================================================================*/

#ifndef __HEAP_H
#define __HEAP_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <iostream>
#include <vector>
#include <algorithm>

/*!
 * \file Heap.h
 * \brief Management of binary heaps
 */

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*!
 * \class Heap
 * \brief class for the management of binary heaps
 *
 * \details This class manages the binary heaps. It is used in particular
 * for fast marching algorithm, but can be used for any other purposes.
 * @remark The operator > must be available for the template parameter T_
 */

template<class T_>
class Heap
{
 public:

// Default constructor
   Heap() { }

// Destructor
   ~Heap() { }

// Insert Element into a Heap
   void insert(T_ el)
   {
      _heap.push_back(el);
      up(_heap.size()-1);
   }

// Remove Minimum Element
   int remove_min()
   {
      if (_heap.size()==0)
         return -1;
      _heap.front();
      _heap[0] = _heap.at(_heap.size()-1);
      _heap.pop_back();
      down(0);
      return 0;
   }

// Extract Minimum Element
   int get_min(T_ el) const 
   { 
      if (_heap.size()==0)
         return -1;
      else 
         el = _heap.front();
      return 0;
   }

// Return heap size
   int size() const { return _heap.size(); }

// Return i-th entry of the heap
   T_ &operator[](int i) { return _heap[i]; }

// Find an element in heap
   int find(T_ el)
   {
      auto it = std::find(_heap.begin(),_heap.end(),el);
      int n = it - _heap.begin();
      if (n<_heap.size())
         return n;
      return -1;
   }

 private:
   std::vector<T_> _heap;

// Return Left Child
   int l(int parent)
   {
      int i = (parent<<1) + 1; // 2 * parent + 1
      return (i<_heap.size()) ? i : -1;
   }

// Return Right Child
   int r(int parent)
   {
      int i = (parent<<1) + 2; // 2 * parent + 2
      return (i<_heap.size()) ? i : -1;
   }

// Return Parent
   int parent(int child)
   {
      if (child==0)
         return -1;
      else
         return (child-1) >> 1;
   }

// Heapify- Maintain Heap Structure bottom up
   void up(int i)
   {
      if (i>=0 && parent(i)>=0 && _heap[parent(i)]>_heap[i]) {
         T_ tmp = _heap[i];
         _heap[i] = _heap[parent(i)];
         _heap[parent(i)] = tmp;
         up(parent(i));
      }
   }

// Heapify- Maintain Heap Structure top down
   void down(int i)
   {
      int child = l(i);
      if (child>0 && r(i)>0 && _heap[child]>_heap[r(i)])
         child = r(i);
      if (child>0) {
         T_ tmp = _heap[i];
         _heap[i] = _heap[child];
         _heap[child] = tmp;
         down(child);
      }
   }
};


template<class T_> 
std::ostream& operator<<(std::ostream &s, Heap<T_> &bh)
{
   for (size_t i=0; i<bh.size(); ++i)
      s << bh[i] << std::endl;
   return s;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} /* namespace OFELI */

#endif
