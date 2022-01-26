/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

                       Definition of for class 'Tabulation'

  ==============================================================================*/

#ifndef __TABULATION_H
#define __TABULATION_H

#include <stdlib.h>
#include <math.h>

#include <string>
using std::string;

#include <map>

#include "OFELI_Config.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Tabulation.h
 *  \brief Definition file for class Tabulation.
 */

/*! \class Tabulation
 * \ingroup Util
 * \brief To read and manipulate tabulated functions.
 * \details This class enables reading a tabulated function of one to three variables and 
 * calculating the value of the function using piecewise multilinear interpolation.\n
 * The file defining the function is an XML file where any function is introduced via the
 * tag "Function".
 * The abcissae are uniformly distributed.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct fct {
   size_t NbVar;
   string Label;
   vector<string>          VarLabel;
   vector<real_t>          Min, Max;
   vector<size_t>          Np;
   Vect<real_t>            Val;
   std::map<string,size_t> Id;
   fct() {
      Min.resize(5);
      Max.resize(5);
      VarLabel.resize(5);
      Np.resize(5);
   }
   fct& operator=(real_t a)
   {
      for (size_t i=0; i<NbVar; ++i)
         Val[i] = a;
      return *this;
   }
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

class Tabulation
{

 public:

/// \brief Default constructor
    Tabulation();

/// \brief Constructor using file name
    Tabulation(string file);

/// \brief Destructor
    ~Tabulation();

/// \brief Set file name.
/// \details This function is to be used when the default constructor is invoked.
    void setFile(string file);

/** \brief Return the calculated value of the function.
 *  \details Case of a function of one variable
 *  @param [in] funct Name of the function to be evaluated, as read from input file
 *  @param [in] x Value of the variable
 *  @return Computed value of the function
 */
    real_t getValue(string funct,
                    real_t x);

/** \brief Return the derivative of the function at a given point.
 *  \details Case of a function of one variable
 *  @param [in] funct Name of the function to be evaluated, as read from input file
 *  @param [in] x Value of the variable
 *  @return Derivative value
 */
    real_t getDerivative(string funct,
                         real_t x);

/** \brief Return the calculated value of the function.
 *  \details Case of a function of two variables
 *  @param [in] funct Name of the function to be evaluated, as read from input file
 *  @param [in] x Value of the first variable
 *  @param [in] y Value of the second variable
 *  @return Computed value of the function
 */
    real_t getValue(string funct,
                    real_t x,
                    real_t y);

/** \brief Return the calculated value of the function.
 *  \details Case of a function of three variables
 *  @param [in] funct Name of the funct to be evaluated, as read from input file
 *  @param [in] x Value of the first variable
 *  @param [in] y Value of the second variable
 *  @param [in] z Value of the third variable
 *  @return Computed value of the function
 */
    real_t getValue(string funct,
                    real_t x,
                    real_t y,
                    real_t z);

/** \brief Return the calculated value of the function.
 *  \details Case of a function of three variables
 *  @param [in] funct Name of the funct to be evaluated, as read from input file
 *  @param [in] x Value of the first variable
 *  @param [in] y Value of the second variable
 *  @param [in] z Value of the third variable
 *  @param [in] t Value of the fourth variable
 *  @return Computed value of the function
 */
    real_t getValue(string funct,
                    real_t x,
                    real_t y,
                    real_t z,
                    real_t t);

/// @brief Get the Number of read functions
/// @return size_t Number of functions
    size_t getNbFuncts() const { return _nb_funct; }

/** @brief Get number of variables of a given function
 *  @param [in] n index of function
 *  @return Number of variables
 */
    size_t getNbVar(size_t n) const { return Funct[n-1].NbVar; }

/** @brief Get the name of a read function
 *  @param [in] n index of function
 *  @return Name of function
 */
    string getFunctName(size_t n) const { return Funct[n-1].Label; }

/** @brief Get number of points defining tabulation
 *  @param [in] n index of function (Starting from 1)
 *  @param [in] i index of variable (Between 1 and 3)
 *  @return Size
 */
    size_t getSize(size_t n, size_t i) const { return Funct[n-1].Np[i-1]; }

/** @brief Get minimal value of a variable
 *  @param [in] n index of function
 *  @param [in] i index of variable (Between 1 and 3)
 *  @return Minimal value
 */
    real_t getMinVar(size_t n,
                     size_t i) const { return Funct[n-1].Min[i-1]; }

/** @brief Get maximal value of a variable
 *  @param [in] n index of function
 *  @param [in] i index of variable (between 1 and 3)
 *  @return 
 */
    real_t getMaxVar(size_t n,
                     size_t i) const { return Funct[n-1].Max[i-1]; }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    friend ostream& operator<<(ostream& s, const Tabulation &t);
    friend class XMLParser;
    vector<fct> Funct;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 private:

   bool                     _uniform;
   size_t                   _nb_funct;
   std::map<string,size_t>  _funct_id;
   void setFunction(string label);
   void setVariable(string label);
   void setSizes();
};

/// \fn ostream & operator<<(ostream& s, const Tabulation &t)
/// \brief Output Tabulated function data.
/// \ingroup Util
    ostream& operator<<(ostream&          s,
                        const Tabulation& t);

} /* namespace OFELI */

#endif
