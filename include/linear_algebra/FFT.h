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

        Definition of class 'FFT' for Fast Fourier Transform calculation

  ==============================================================================*/

#ifndef __FFT_H
#define __FFT_H

/*!
 * \file FFT.h
 * \brief Fast Fourier Transform
 */

#include "linear_algebra/Vect.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*!
 * \class FFT
 * \ingroup VectMat
 * \brief class for the Fast Fourier Transform calculation
 * \details This class enables computing an in-place one dimensional complex discrete
 * Fourier transform. It uses an FFT algorithm for a number of 
 * input points which is a power of two.
 *
 *  We must set the size of the sine and cosine table used in the FFT.
 *  The table is called starting_multiplier.
 *     
 *
 *  to be the size of the table.  This table can now be used by the FFT for 
 *  any number of points from 2 up to MAX.
 *
 *  For example, if MAX_TABLE_SIZE = 14, then we can transform anywhere from 
 *  2 to 2^15 = 32,768 points, using the same sine and cosine table.
 *
 *  @example 
 *  Consider the vector \c v given by
 *
 *                            ^
 *   The Fourier coefficients v[k] are defined for k = 0, . . ., N - 1
 *   by the formula
 *
 *                       N-1           2 Pi i
 *   ^           1      ____         - ------ j k
 *   v[k] =  --------   \                N
 *              ---     /     x(j) e
 *           \ / N      ----
 *            v          j=0
 *
 *   where e is 2.1718..., Pi = 3.14159... and i is the square root of -1.
 *
 *   This is the formula for the forward transform.  The inverse transform
 *   is the same, except that the minus sign in the exponent above is a
 *   plus sign instead.
 *                                 ---
 *    The normalizing factor 1 / \/ N  in front was picked so that the
 *    forward and inverse transforms look the same.
 *
 *    The particular type of FFT algorithm we use is a
 *
 *                              Complex,
 *                            Two to the N,
 *                       Decimation in  frequency,
 *                 Input in order, output bit reversed
 *
 *    type of FFT.  It is described in the book
 *
 *                  THE FAST FOURIER TRANSFORM,
 *                 F. Oran Brigham, Prentice Hall, 1974
 *
 * \author 
 * Sean O'Connor artificer!AT!seanerikoconnor!DOT!freeservers!DOT!com
 *
 * BUGS
 * 
 *   The bit-reversal routine could be made faster by using bit shifting operations.
 *   You might want to use long double precision complex numbers to get more accuracy.
 *   Check for a power of 2 probably ought to use the machine's floating point
 *   epsilon (computed automatically).
 */

class FFT : public Vect<complex_t>
{
   public:

/** \brief Default constructor
 *  \details This constructor just calls the base class constructor for
 *  the vector type, then initializes the derived class field.
 */
    FFT();

/// \brief Destructor
    ~FFT() { }

/** \brief Constructor using a double precision real vector
 *  \details Using this constructor, this class computes the FFT of the
 *  given vector as a class derived from Vect<complex<double> >
 *
 *  @param [in] v Vect instance containing a real valued vector for which 
 *  the Fourier transform is evaluated
 */
    FFT(const Vect<real_t>& v);

/// \brief Copy constructor
    FFT(const FFT& t);

/// \brief Assignment operator
    FFT& operator=(const FFT& t);

/** \brief Run the FFT algorithm
 *  \details This member function does an in-place one dimensional complex
 *  discrete Fourier transform. It uses an FFT algorithm for a number
 *  of input points which is a power of two. As a result, the instance of
 *  class FFT is derived from class Vect.
 */
    void run() { run(true); }

/** \brief Run the inverse FFT algorithm
 *  \details This member function does an in-place one dimensional complex
 *  discrete inverse Fourier transform. It uses an FFT algorithm for a number
 *  of input points which is a power of two. As a result, the instance of
 *  class FFT is derived from class Vect.
 */
    void run_inv() { run(false); }

 private:
    bool _dir;

//  Table of sines and cosines whose values are created the first 
//  time any object calls the fft() member function and which are
//  shared among all subsequent objects.
    complex_t _starting_multiplier[MAX_FFT_SIZE];
    bool _called_already;

/** This function reverses the order of the bits in a number.
 *  INPUT
 *    n              The number whose bits are to be reversed.
 *    num_bits       The number of bits in N, counting leading zeros.
 *
 *  OUTPUT
 *                   The number N with all its num_bits bits reversed.
 *
 *  EXAMPLE CALLING SEQUENCE
 *
 *     reverse_bits( (10)  , 2 ) = (01) = 1, but
 *                       2             2
 *    reverse_bits( (10)  , 3 ) = reverse_bits( (010)  , 3 ) = (010)  = (10)
 *                      2                            2              2       2
 *    = 10
 *
 *  METHOD
 *
 *     The algorithm works by a method similar to base conversion.
 *     As the low order bits are broken off, one by one, they are
 *     multiplied by 2 and accumulated to generate the bit-reversed
 *     number.
 *
 * BUGS
 *
 *    The bit-reversal routine could be made faster by using C shift operators.
 */
    static int reverse_bits(int n, int num_bits);
    void run(bool dir);
};


/*! @} End of Doxygen Groups */

} /* namespace OFELI */


#endif
