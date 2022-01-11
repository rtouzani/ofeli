#! /bin/sh
# install-sdk.sh

PREFIX="/usr/local"

echo "Execution of this script my require administrator privileges"
echo "In this case, type: sudo install.sh"

mkdir -p ${PREFIX}/bin
cp bin/* ${PREFIX}/bin/.
mkdir -p ${PREFIX}/share/ofeli/material
cp material/* ${PREFIX}/share/ofeli/material/.
mkdir -p ${PREFIX}/lib
cp -a lib/libofeli.a ${PREFIX}/lib/.
mkdir -p ${PREFIX}/include/ofeli
cd include/
cp -rf * ${PREFIX}/include/ofeli/
cd ..
echo "ofeli installation complete"
