/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

                            Definition of class 'Fct'
                    Class to define functions by an expression

  ==============================================================================*/

#ifndef __FCT_H
#define __FCT_H

#include <string>
#include <vector>
using std::string;
using std::vector;

#include "OFELI_Config.h"
#include "linear_algebra/Point.h"

namespace exprtk {
   template <typename T> class parser;
   template <typename T> class expression;
   template <typename T> class expression_helper;
   template <typename T> class symbol_table;
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

class Fct {

 public:

   Fct();
   Fct(const string& exp, const vector<string>& v);
   Fct(const string& exp);
   Fct(const string& exp, const string& v);
   Fct(const string& n, const string& exp, const vector<string>& v);
   ~Fct();
   string getErrorMessage();
   int set(const string& exp, const vector<string>& v, int opt=0);
   int set(const string& exp, int opt=0);
   int set(const string& exp, const string &v, int opt=0);
   int set(const string& n, const string &exp, const vector<string>& v, int opt=0);
   real_t operator()(real_t x);
   real_t operator()(real_t x, real_t y);
   real_t operator()(real_t x, real_t y, real_t z);
   real_t operator()(real_t x, real_t y, real_t z, real_t t);
   real_t operator()(const Point<real_t>& x);
   real_t operator()(const Point<real_t>& x, real_t t);
   real_t operator()(const vector<real_t>& x);
   real_t D(real_t x);
   real_t D(const vector<real_t>& x, size_t i);
   int check();
   string name, expr;
   vector<string> var;
   size_t nb_var;
   friend std::ostream& operator<<(std::ostream& s, const Fct& f);

 private:

   exprtk::parser<real_t> *_p;
   exprtk::symbol_table<real_t> *_st;
   exprtk::expression<real_t> *_ex;
   void add_constants();
   bool exp_ok, var_ok;
   int err;
   string error_message;
   vector<real_t> xvar;
};

} /* namespace OFELI */
/*! @} End of Doxygen Groups */
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif
