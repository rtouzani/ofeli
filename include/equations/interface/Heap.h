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

                 Definition of class 'Heap' for heap management

  ==============================================================================*/

#ifndef __HEAP_H
#define __HEAP_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/*!
 * \file Heap.h
 * \brief The Heap
 * \author Boris Meden, Mohamed Sylla
 * \version 0.1
 */

#include "equations/interface/IPoint.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*!
 * \class Heap
 * \brief Heap for priority queue
 *
 * \details This class manages the priority queue
 */
class Heap
{
 private:
   Vect<IPoint> _heap;
   size_t _max_size, _size;

 public:
/*!
 * \brief Constructor
 * \param [in] size the heap size
 */
    Heap(size_t size = 0);

/*!
 * \brief Constructor
 * Constructor using a vector of point of type IPoint
 * \param [in] vec vector of point of type IPoint 
 */
   Heap(Vect<IPoint>& vec);

/*!
 * \brief Constructor
 * Constructor by copy
 * \param [in] tas Heap to copy from
 */
   Heap(Heap& tas);

/// \brief Destructor
   ~Heap();

/// \brief Operator <tt>=</tt>
/// \param [in] tas Heap to copy from
    Heap & operator=(const Heap& H);

/*!
 * \brief Operator <tt>[]</tt>
 * Read only
 * \param [in] index in the Heap
 * \return point at the index <tt>ind</tt>
 */
    IPoint operator[](size_t ind) const;

/*!
 * \brief Operator <tt>[]</tt>
 * \param [in] index in the Heap
 * \return point at the index <tt>ind</tt>
 */
    IPoint &operator[](size_t ind);

/// \brief Set Heap size
/// \param [in] size heap new size
    void setSize(size_t size) { _max_size = size; }

    void setHeap(Vect<IPoint> vec);
    void Pop() { _size--; }

/// \brief return Heap size
/// \return Heap size
    size_t getSize() const { return _size; }

/// \brief extract the heap's first element
/// \return element at the top of the Heap
    IPoint Current();

/// \brief Add element to the Heap
/// \param [in] el node to add to the Heap
    void Add(IPoint el);

/*! \brief Update the value of element at index <tt>rg</tt>
 * \param [in] val the new value of the node at index <tt>rg</tt>
 * \param [in] rg the index of the node to update
 */
    void Update(const real_t& val,
                      size_t  rg);

/*!
 * \brief Find a node in the Heap
 * \param [in] pt node to find
 * \param [out] ind the node index in the Heap if it's found
 * \return true if the node is in the Heap, false otherwise
 */
    bool Find(const IPoint& pt,
                    size_t& ind);

    void Down_Heap(size_t rank = 1);
    void Up_Heap(size_t rank);

    friend ostream & operator<< (ostream&    s,
                                 const Heap& H);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif
