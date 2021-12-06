#! /bin/sh
# install.sh

PREFIX="/usr/local"

mkdir -p $PREFIX/bin
cp bin/* $PREFIX/bin/.
mkdir -p $PREFIX/share/ofeli/material
cp material/* $PREFIX/share/ofeli/material/.
mkdir -p $PREFIX/lib
cp -a lib/libofeli.a $PREFIX/lib/.
mkdir -p $PREFIX/include/ofeli
cp -a lib/libofeli.a $PREFIX/include/ofeli/.
cp -rf include/ $PREFIX/include/ofeli/

echo "ofeli installation complete"
