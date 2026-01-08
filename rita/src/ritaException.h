/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2026 Rachid Touzani

    This file is part of rita.

    rita is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rita is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

  ==============================================================================

      Definition of class 'ritaException' for handling exceptions in rita

  ==============================================================================*/

#pragma once

#include <stdexcept>

class ritaException : public std::runtime_error
{

 public:

/// \brief This form will be used most often in a throw.
    ritaException(const std::string& s) : runtime_error(s)
    { };

/// \brief Throw with no error message.
    ritaException() : std::runtime_error("Exception thrown in rita:\n") { }; 
};

#define CATCH_RITA_EXCEPTION catch(ritaException &e) {                                 \
                                std::cout << "rita error: " << e.what() << endl;       \
                                return 1;                                              \
                             }                                                         \
                             catch(runtime_error &e) {                                 \
                                std::cout << "Runtime error: " << e.what() << endl;    \
                                return 1;                                              \
                             }                                                         \
                             catch( ... ) {                                            \
                                std::cout << "Unexpected error: " << endl;             \
                                return 1;                                              \
                             }
