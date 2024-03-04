#! /bin/sh
# install-sdk.sh
# A script to install the ofeli Software Development Kit
# the following line can be modified according to the installator usage
PREFIX="/usr/local"

echo "Execution of this script my require administrator privileges"
echo "For this, type: sudo install.sh"

mkdir -p ${PREFIX}/bin
cp bin/rita ${PREFIX}/bin/.
cp bin/gmsh ${PREFIX}/bin/.
mkdir -p ${PREFIX}/share/ofeli/material
cp material/* ${PREFIX}/share/ofeli/material/.
mkdir -p ${PREFIX}/lib
cp -a lib/libofeli.a ${PREFIX}/lib/.
cp -a lib/libgmsh.* ${PREFIX}/lib/
mkdir -p ${PREFIX}/include/ofeli
cd include/
cp -rf * ${PREFIX}/include/ofeli/
cd ..
echo "ofeli installation complete"
