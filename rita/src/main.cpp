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

                             Main calling program

  ==============================================================================*/


#define _USE_MATH_DEFINES 
#undef __STRICT_ANSI__

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <typeinfo>
#include <ctime>

//--- muparserx framework -------------------------------------------------------------------------
#include "../muparserx/mpParser.h"
#include "../muparserx/mpDefines.h"
#include "../muparserx/mpTest.h"

//--- rita includes -------------------------------------------------------------------------------
#include "rita.h"

using namespace std;
using namespace mup;

const string currentDateTime() {
   time_t     now = time(0);
   struct tm  tstruct;
   char       buf[80];
   tstruct = *localtime(&now);
   strftime(buf,sizeof(buf),"%Y-%m-%d  %X",&tstruct);
   return buf;
}

int main(int argc, char *argv[])
{
   std::time_t t = std::time(0);
   std::tm* now = std::localtime(&t);
   string release = RITA_VERSION;
   string year = RITA_DATE;
   cout << "\n     R I T A     " << release << endl;
   cout << "     Last update " << (now->tm_year+1900) << '-' << (now->tm_mon+1)
        << '-' <<  now->tm_mday << "\n";
   cout << "     Copyright (C) 2021-" << year << ", Rachid Touzani\n\n";
   cout << "Type \"help\" or \"license\" for more information." << endl;
   cout << "rita web site:   http://ofeli.org/rita\n" << endl;
   int ret = 0;

   try {
      RITA::rita r;
      r.setRelease(release,year);
      if (argc>1) {
         if (string(argv[1])=="-h") {
            cout << r.hh->H0;
            return 0;
         }
         else if (string(argv[1])=="--help") {
            cout << r.hh->H1;
            return 0;
         }
         else if (string(argv[1])=="-v" || string(argv[1])=="--version") {
            cout << "rita, Release " << release << ", " << year << endl;
            return 0;
         }
         else if (argv[1][0]=='-') {
            cout << "Input error: Unknown argument." << endl;
            return 1;
         }
         else
            r.setInput(string(argv[1]),0);
      }
      ret = r.run();
   } CATCH_RITA_EXCEPTION
   return ret;
}
