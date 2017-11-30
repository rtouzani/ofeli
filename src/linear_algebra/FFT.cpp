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

                          Implementation of class 'FFT'

  ==============================================================================*/

#include "util/constants.h"
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;

#include "linear_algebra/FFT.h"
#include "OFELIException.h"

namespace OFELI {

FFT::FFT() 
    : Vect<complex_t>(), _dir(false), _called_already(false)
{ }


FFT::FFT(const Vect<real_t>& v) 
    : Vect<complex_t>(v.size()), _dir(false), _called_already(false)
{
   for (size_t i=0; i<size(); ++i)
      (*this)[i] = complex_t(v[i],0.0);
}


FFT::FFT(const FFT& t) 
    : Vect<complex_t>(t), _called_already(false)
{
   _dir = t._dir;
}


FFT & FFT::operator=(const FFT& t)
{
   if (this==&t)
      return *this;
   FFT::operator=(t);
   _dir = t._dir;
   return *this;
}


void FFT::run(bool dir) 
{
// Number of complex points in the FFT from the vector's size.
   size_t N = size();

/*  Compute log base 2 of N. Perturb slightly to avoid roundoff error.
 *  For example, int(4.99999999999) = 4 (wrong!), but
 *  (int)( 4.999999999 + .01 ) = 5 (correct!)
 */
   int logn = int(log(real_t(N))/(log(2.0))+0.01);

/*  Error checking.  Print a message to the terminal if any one of these
 *  conditions is true:
 *
 *       N < 2.
 *       log2(N) > MAX_FFT_SIZE + 1.
 *       N is not a power of 2.
 */
   if (N<2)
      throw OFELIException("In FFT::run(), FFT has less than two points");

   else if ((logn-1) > MAX_FFT_SIZE)
      throw OFELIException("In FFT::run(), FFT has too many points for its sin/cos table");

   else if (int(pow(2.0,logn)) != N)
      throw OFELIException("In FFT::run(), Number of points in the FFT is not a power of 2");

/*
 *  If _called_already is false, generate the sine and cosine table for
 *  the first time, then set _called_already to true.  Now later calls to 
 *  this function don't recreate the sine and cosine table.
 */
   if (!_called_already) {
      for (int i=0; i<=MAX_FFT_SIZE-1; ++i) { 
//       Current angle in sine and cosine table.
         real_t angle = -OFELI_PI / pow(real_t(2),real_t(i));

//       e^t = cos t + i sin t
         _starting_multiplier[i] = complex_t(cos(angle),sin(angle));
      }
      _called_already = true;
   }

/*  Compute the FFT according to its signal flow graph with 
 *  computations arranged into butterflies, groups of butterflies
 *  with the same multiplier, and columns.
 */
   int butterfly_height=N/2, number_of_groups=1;

/*
 *  A column refers to a column of the FFT signal flow graph.  A group is
 *  a collection of butterflies which have the same multiplier W^p.  A
 *  butterfly is the basic computational unit of the FFT.
 */
   for (int k=1; k<=logn; ++k) {
//    The multiplier W^p for a group of butterflies.
      complex_t multiplier = _starting_multiplier[k-1];

//    Modify starting multiplier W^p when computing the inverse transform.
      if (dir==false)
         multiplier = conj(multiplier);

//    Difference between current multiplier and next.
      complex_t increment = multiplier;

//    Compute the first group.
      for (int i=0; i<=butterfly_height-1; ++i) {
//       Lower branch of a butterfly.
         int offset = i + butterfly_height;

         complex_t temp = (*this)[i] + (*this)[offset];
         (*this)[offset] = (*this)[i] - (*this)[offset];
         (*this)[i] = temp;
      }

//    Now do the remaining groups.
      for (int group=2; group<=number_of_groups; ++group) {
//       Array index of first butterfly in group.
         int start_of_group = reverse_bits(group-1,logn);

//       Array index of last butterfly in group. 
         int end_of_group = start_of_group + butterfly_height - 1;

//       Compute all butterflies in a group.
         for (int i=start_of_group; i<=end_of_group; ++i) {
            int offset = i + butterfly_height;

//          Temporary storage in a butterfly calculation.
            complex_t product = (*this)[offset]*multiplier;
            complex_t temp    = (*this)[i] + product;
            (*this)[offset]   = (*this)[i] - product;
            (*this)[i]        = temp;
         }
         multiplier *= increment;
      }

      butterfly_height = butterfly_height/2;
      number_of_groups = number_of_groups*2;
   }

// Normalize the results by \/ N and permute them into order, two at a time.
   double normalizer = double(1.0)/sqrt(double(N));

   for (int i=0; i<=N-1; ++i) {
//    Bit-reversed array index i.
      int rev_i = reverse_bits(i,logn);

      if (rev_i>i) { // Exchange is not yet done.
         complex_t temp = (*this)[i];
         (*this)[i]     = (*this)[rev_i]*normalizer;
         (*this)[rev_i] = temp*normalizer;
      }
      else if (rev_i==i) // exchange with itself
         (*this)[i] = (*this)[i]*normalizer;
   }
   return;
}


int FFT::reverse_bits(int n,
                      int num_bits)
{
   int remainder = n,   // Bits remaining to be reversed.
       reversed_n = 0;  // Bit-reversed value of n.
   for (int bit_num=1;  bit_num<=num_bits; ++bit_num) {
      int bit    = remainder%2;           // Next least significant bit.
      reversed_n = bit + 2*reversed_n;    // Shift left and add to buffer.
      remainder  = int(remainder/2);      // Remaining bits.
   }
   return(reversed_n) ;
}

} /* namespace OFELI */
