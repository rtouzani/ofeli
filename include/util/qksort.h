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

                           Functions for quick sorting

  ==============================================================================*/

#ifndef __QKSORT_H
#define __QKSORT_H

#include <vector>

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file qksort.h
 *  \brief File that contains template quick sorting function.
 */


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
bool _less(const T_ &x, const T_ &y);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \fn void QuickSort(std::vector<T_> &a, int begin, int end)
 * \ingroup Util
 * \brief Function to sort a vector.
 *
 * \details qksort uses the famous quick sorting algorithm.
 * @param [in,out] a Vector to sort.
 * @param [in] begin index of starting iterator
 * @param [in] end index of ending iterator
 *
 * The calling program must provide an overloading of the operator <
 * for the type \b T_
 */
template<class T_>
void QuickSort(std::vector<T_>& a,
               int              begin,
               int              end)
{
   if (begin >= end)
      return;
   size_t pivot = begin + (rand() % (end+1-begin));
   _swap(a[pivot],a[end]);

   T_ val = a[end];
   int left=begin-1, right=end;
   for (;;) {
      do left++;
      while (a[left] < val)
         ;
      do right--;
      while (right>=begin && val < a[right])
         ;
      if (left >= right)
         break;
      _swap(a[left],a[right]);
   }
   T_ tmp = a[left];
   a[left] = val;
   a[end] = tmp;
   qksort(a, begin, left-1);
   qksort(a, left+1, end);
}


/** \fn void qksort(std::vector<T_> &a, int begin, int end)
 * \ingroup Util
 * \brief Function to sort a vector.
 *
 * \details qksort uses the famous quick sorting algorithm.
 * @param [in,out] a Vector to sort.
 * @param [in] begin index of starting index (default value is 0)
 * @param [in] end index of ending index (default value is the vector size - 1)
 */
template<class T_>
void qksort(std::vector<T_>& a,
            int              begin,
            int              end)
{
   if (begin >= end)
      return;
   size_t pivot = begin + (rand() % (end+1-begin));
   _swap(a[pivot],a[end]);

   T_ val = a[end];
   int left=begin-1, right=end;
   for (;;) {
      do left++;
      while (a[left]<val)
         ;
      do right--;
      while (right>=begin && val<a[right])
         ;
      if (left >= right)
         break;
      _swap(a[left],a[right]);
   }
   T_ tmp = a[left];
   a[left] = val;
   a[end] = tmp;
   qksort(a, begin, left-1);
   qksort(a, left+1, end);
}


/** \fn void qksort(std::vector<T_> &a, int begin, int end, C_ compare)
 * \ingroup Util
 * \brief Function to sort a vector according to a key function.
 *
 * \details qksort uses the famous quick sorting algorithm.
 * @param [in,out] a Vector to sort.
 * @param [in] begin index of starting index (0 for the beginning of the vector)
 * @param [in] end index of ending index
 * @param [in] compare A function object that implements the ordering. The user
 * must provide this function that returns a boolean function that is true if
 * the first argument is less than the second and false if not.
 */
template<class T_, class C_>
void qksort(std::vector<T_>& a,
            int              begin,
            int              end,
            C_               compare)
{
   if (begin >= end)
      return;
   size_t pivot = begin + (rand() % (end+1-begin));
   _swap(a[pivot],a[end]);

   T_ val = a[end];
   int left=begin-1, right=end;
   for (;;) {
      do left++;
      while (compare(a[left], val))
         ;
      do right--;
      while (right>=begin && compare(val,a[right]))
            ;
      if (left >= right)
         break;
      _swap(a[left],a[right]);
   }
   T_ tmp = a[left];
   a[left] = val;
   a[end] = tmp;
   qksort(a, begin, left-1, compare);
   qksort(a, left+1, end, compare);
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T_>
inline void _swap(T_& a,
                  T_& b)
{
   T_ tmp = a;
   a = b;
   b = tmp;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
