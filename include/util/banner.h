/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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
                               Display ofeli's banner
  ==============================================================================*/

#ifndef __BANNER_H
#define __BANNER_H

#include "OFELI_Config.h"

#include <string>
using std::string;

#include <iostream>
using std::cout;
using std::endl;

/*! \file banner.h
 *  \brief A function to output a banner.
 */

namespace OFELI {

/** \fn void banner(const string &prog=" ")
 * \ingroup Util
 * \brief Outputs a banner as header of any developed program.
 * @param [in] prog Calling program name.
 * Enables writing a copyright notice accompanying the program.
 */
inline void banner(const string &prog=" ")
{
   if (prog[0] != ' ') {
      cout << endl;
      cout << "====================================================================" << endl;
      cout << prog << ", Copyright (c) 1998 - 2016  Rachid Touzani\n\n";
      cout << "This program is free software: you can redistribute it and/or modify\n";
      cout << "it under the terms of the GNU Lesser General Public License as published by\n";
      cout << "the Free Software Foundation, either version 3 of the License, or\n";
      cout << "(at your option) any later version.\n\n";
      cout << "This program is distributed in the hope that it will be useful,\n";
      cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
      cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
      cout << "GNU Lesser General Public License for more details.\n\n";
      cout << "You should have received a copy of the GNU Lesser General Public License\n";
      cout << "along with this program.  If not, see <http://www.gnu.org/licenses/>.\n";
      cout << "====================================================================" << endl;
      cout << prog << " uses the ofeli package" << endl;
   }
   else {
      cout << endl;
      cout << "====================================================================" << endl;
      cout << "Copyright (c) 1998 - 2016 by Rachid Touzani\n\n";
      cout << "This is free software, and your are allowed to redistribute it\n";
      cout << "under certain conditions. Details are distributed with the software." << endl << endl;
      cout << "====================================================================" << endl;
   }
   cout << "version : " << OFELI_VERSION << endl;
   cout << "Date of Release " << OFELI_RELEASE_DATE << endl;
   cout << "---------------------------------------------------" << endl;
   cout << "Date of latest library building : " << __DATE__ << endl;
   cout << "====================================================================" << endl;
   cout << "ofeli comes with ABSOLUTELY NO WARRANTY.\n";
   cout << "====================================================================" << endl << endl;
}

} /* namespace OFELI */

#endif
