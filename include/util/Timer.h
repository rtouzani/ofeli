/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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

                            Definition of class 'Timer'

  ==============================================================================*/


#ifndef __TIMER_H
#define __TIMER_H

#include <time.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/timeb.h>

#include <iostream>

#include "OFELI_Config.h"
#include "OFELIException.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Timer.h
 *  \brief Definition file for class Timer.
 */

/*! \class Timer
 *  \ingroup Util
 *  \brief To handle elapsed time counting
 *  \details This class is to be used when testing program performances.
 *  A normal usage of the class is, once an instance is constructed,
 *  to use alternatively, Start, Stop and Resume. Elapsed time can be obtained
 *  once the member function Stop is called.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Timer
{

 public:

/// \brief Default constructor
    Timer()
    {
       _ct = 0;
       _st0 = _st1 = false;
    }

/// \brief Destructor
    ~Timer() { }

/// \brief Say if time counter has started
/// \details Return true if time has started, false if not 
    bool Started() const
    {
       return _st0;
    }
   
/** \brief Start (or resume) time counting
 *  \details This member function is to be used to start or resume time counting
 */
    void Start()
    { 
       _t0 = clock();
       _st0 = true;
    }

/** \brief Stop time counting
 *  \details This function interrupts time counting. This one can be resumed
 *  by the function Start
 */
    void Stop()
    {
       if (!_st0)
          throw OFELIException("Timer::Stop(): The function Start() must have been called before.");
       _t1 = clock();
       _ct += _t1-_t0;
       _st1 = true;
    }

/// \brief Clear time value (Set to zero)
    void Clear()
    { 
       _t0 = clock();
       _ct = 0;
    }

/// \brief Return elapsed time (in seconds)
    real_t get() const
    {
       if (!_st0 || !_st1)
         throw OFELIException("Timer::get(): The function Start() must have been called before.");
       return _ct/real_t(CLOCKS_PER_SEC);
   }
   
/// \brief Return elapsed time (in seconds)
/// \details Identical to get
    real_t getTime() const { return get(); }
   
 private:
    clock_t _t0, _t1, _ct;
    bool    _st0, _st1;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
