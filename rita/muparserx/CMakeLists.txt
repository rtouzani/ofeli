# ==============================================================================
#
#                                 r  i  t  a
#
#            An environment for Modelling and Numerical Simulation
#
# ==============================================================================
#
#   Copyright (C) 2021 - 2022  Rachid Touzani
#
#   This file is part of rita.
#
#   rita is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   rita is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
# ==============================================================================

target_sources (${PROJECT_NAME} PRIVATE
                mpError.cpp
                mpFuncCmplx.cpp			
                mpFuncCommon.cpp		
                mpFuncMatrix.cpp		
                mpFuncNonCmplx.cpp		
                mpFuncStr.cpp			
                mpICallback.cpp			
                mpIOprt.cpp			
                mpIOprtBinShortcut.cpp		
                mpIPackage.cpp			
                mpIToken.cpp			
                mpIValReader.cpp		
                mpIValue.cpp			
                mpIfThenElse.cpp		
                mpOprtBinAssign.cpp		
                mpOprtBinCommon.cpp		
                mpOprtBinShortcut.cpp		
                mpOprtCmplx.cpp			
                mpOprtIndex.cpp			
                mpOprtMatrix.cpp
                mpOprtNonCmplx.cpp
                mpOprtPostfixCommon.cpp
                mpPackageCmplx.cpp
                mpPackageCommon.cpp
                mpPackageMatrix.cpp
                mpPackageNonCmplx.cpp
                mpPackageStr.cpp
                mpPackageUnit.cpp
                mpParser.cpp
                mpParserBase.cpp
                mpParserMessageProvider.cpp
                mpRPN.cpp
                mpScriptTokens.cpp
                mpTest.cpp
                mpTokenReader.cpp
                mpValReader.cpp
                mpValue.cpp
                mpValueCache.cpp
                mpVariable.cpp
               )