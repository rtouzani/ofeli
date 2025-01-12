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

           Definition and implementation of class 'FMHeap' to manage 
               binary heaps for the Fast Marching algorithm

  ==============================================================================*/

#ifndef __FM_HEAP_H
#define __FM_HEAP_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include "OFELI_Config.h"
#include <iostream>
#include <vector>
#include <algorithm>

/*!
 * \file FMHeap.h
 * \brief Management of binary heaps for Fast marching
 */

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*!
 * \struct Pt
 * \brief class for the management of points in grid
 *
 */

struct Pt {

   int i, j, k, sgn, state;
   real_t x, y, z, v;

   enum { FROZEN, ALIVE, FAR } st;
   Pt(): i(0), j(0), k(0) { }
   Pt(int ii, int jj): i(ii), j(jj), k(0), state(FAR) { } 
   Pt(int ii, int jj, real_t vv): i(ii), j(jj), k(0), state(FAR), v(vv) { } 
   Pt(int ii, int jj, int kk): i(ii), j(jj), k(kk), state(FAR) { } 
   Pt(int ii, int jj, int kk, real_t vv): i(ii), j(jj), k(kk), state(FAR), v(vv) { } 

};


inline std::ostream& operator<<(std::ostream &s, const Pt &p)
{
   if (p.k==0)
      s << "(" << p.i << "," << p.j << "):  value: " << p.sgn*p.v;
   else
      s << "(" << p.i << "," << p.j << "," << p.k << "):  value: " << p.sgn*p.v;
   if (p.state==Pt::ALIVE)
      s << ", status: alive";
   else if (p.state==Pt::FAR)
      s << ", status: far";
   else if (p.state==Pt::FROZEN)
      s << ", status: frozen";
   return s;
}

/*!
 * \class FMHeap
 * \brief class for the management of binary heaps for the Fast Marching algorithm
 *
 * \details This class manages the binary heaps. It is used in particular
 * for fast marching algorithm.
 */

class FMHeap
{
 public:

/// Default constructor
    FMHeap() { }

/// Destructor
    ~FMHeap() { }

/// Insert Element into a Heap
    void insert(Pt *p)
    {
       _heap.push_back(p);
       up(_heap.size()-1);
    }

/// Remove Minimum Element
    int remove()
    {
       if (_heap.size()==0)
          return -1;
//       _heap.front();
       _heap[0] = _heap.at(_heap.size()-1);
       _heap.pop_back();
       down(0);
       return 0;
    }

/// Extract Minimum Element
    Pt *get() const
    { 
       if (_heap.size()==0)
          return nullptr;
       return _heap.front();
    }

/// Return heap size
    int size() const { return _heap.size(); }

/// Return i-th entry of the heap
    Pt * &operator[](int i) { return _heap[i]; }

/// Find an element in heap
/// @param [in] p Pointer to element to find
    int find(Pt *p)
    {
       auto it = std::find(_heap.begin(),_heap.end(),p);
       int n = it - _heap.begin();
       if (n<size())
          return n;
       return -1;
    }

 private:
   std::vector<Pt *> _heap;

// Return Left Child
   int l(int parent)
   {
      size_t i = (parent<<1) + 1; // 2 * parent + 1
      return (i<_heap.size()) ? i : -1;
   }

// Return Right Child
   int r(int parent)
   {
      size_t i = (parent<<1) + 2; // 2 * parent + 2
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
      if (i>=0 && parent(i)>=0 && _heap[parent(i)]->v>_heap[i]->v) {
         Pt *tmp = _heap[i];
         _heap[i] = _heap[parent(i)];
         _heap[parent(i)] = tmp;
         up(parent(i));
      }
   }

// Heapify- Maintain Heap Structure top down
   void down(int i)
   {
      int child = l(i), child1 = r(i);
      if (child>=0 && child1>=0 && _heap[child]->v>_heap[child1]->v)
         child = child1;
      if (child>0) {
         Pt *tmp = _heap[i];
         _heap[i] = _heap[child];
         _heap[child] = tmp;
         down(child);
      }
   }
};


/// Output heap contents
    inline std::ostream& operator<<(std::ostream &s, FMHeap &h)
    {
       for (int i=0; i<h.size(); ++i)
          s << *h[i] << std::endl;
       return s;
    }

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif
