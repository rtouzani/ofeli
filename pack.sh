#! /bin/sh
# pack.sh 
# A script to build OFELI distribution

if [ $# -eq 0 ]; then
   echo 1>&2 "Usage: pack.sh <release> [MacOSX|Linux64|Win64]"
   exit 2
fi

RELEASE=$1
REL=$1
PREFIX="/usr/local"
SYSTEM=$2
RELEASE=ofeli-${RELEASE}

if [${SYSTEM} == ""]; then
   echo "Creating Source Package ..."
   tar --exclude='./build' --exclude='.DS_Store' --exclude='.git' --exclude='.vscode' -czf ofeli.tar.gz ofeli
   mv ofeli ${RELEASE}
   tar --exclude='./build' --exclude='.DS_Store' --exclude='.git' --exclude='.vscode' --exclude='./pack.sh' --exclude='install.sh' -czf ${RELEASE}-src.tar.gz ${RELEASE}
   mv ${RELEASE} ofeli

   echo "Creating Demo Package ..."
   DEMO=${RELEASE}-demos
   mkdir -p ${DEMO}
   cd ${DEMO}
   cp -rf ${PREFIX}/share/ofeli/demos/* .
   find . -name "CMakeLists.txt" -exec rm {} \;
   cd ..
   tar czf ${DEMO}.tar.gz ${DEMO}/
   rm -rf ${DEMO}

   echo "Creating Doc Package ..."
   DOC=${RELEASE}-doc
   mkdir -p ${DOC}
   cd ${DOC}
   cp -rf "$PREFIX"/share/ofeli/doc/ .
   find . -name "CMakeLists.txt" -exec rm {} \;
   cd ..
   tar czf ${DOC}.tar.gz ${DOC}/
   rm -rf ${DOC}
   exit
fi

if [${SYSTEM} != "MacOSX"] &&
   [${SYSTEM} != "Linux64"] &&
   [${SYSTEM} != "Win64"]; then
     echo "Error: Unavailable Platform $MACHINE"
   exit
fi

echo "Preparing release $RELEASE..."
echo "--------------------------"


case $SYSTEM in

    MacOSX)
        mkdir -p ${RELEASE}
        mkdir -p ${RELEASE}/bin
        mkdir -p ${RELEASE}/lib
        mkdir -p ${RELEASE}/include
        mkdir -p ${RELEASE}/material
        cp ${PREFIX}/lib/libofeli.a ${RELEASE}/lib/.
        cp ${PREFIX}/bin/cmesh ${RELEASE}/bin/.
        cp ${PREFIX}/bin/cfield ${RELEASE}/bin/.
        cp ${PREFIX}/bin/g2m ${RELEASE}/bin/.
        cp ${PREFIX}/bin/vmesh ${RELEASE}/bin/.
        cp ${PREFIX}/share/ofeli/material/* ${RELEASE}/material/.
        cp -rf ${PREFIX}/include/ofeli/* ${RELEASE}/include/.
        cp ofeli/install.sh ${RELEASE}/install.sh
        cp ofeli/README-sdk.md ${RELEASE}/README.md
        find ${RELEASE}/ -name "CMakeLists.txt" -exec rm {} \;
        tar czf ${RELEASE}-${SYSTEM}.tar.gz ${RELEASE}/
        rm -rf ${RELEASE}
        ;;

    Linux64)
        mkdir -p ${RELEASE}
        mkdir -p ${RELEASE}/bin
        mkdir -p ${RELEASE}/lib
        mkdir -p ${RELEASE}/include
        mkdir -p ${RELEASE}/material
        cp ${PREFIX}/lib/libofeli.a ${RELEASE}/lib/.
        cp ${PREFIX}/bin/cmesh ${RELEASE}/bin/.
        cp ${PREFIX}/bin/cfield ${RELEASE}/bin/.
        cp ${PREFIX}/bin/g2m ${RELEASE}/bin/.
        cp ${PREFIX}/bin/vmesh ${RELEASE}/bin/.
        cp ${PREFIX}/share/ofeli/material/* ${RELEASE}/material/.
        cp -rf ${PREFIX}/include/ofeli/* ${RELEASE}/include/.
        cp ofeli/install.sh ${RELEASE}/install.sh
        cp ofeli/README-sdk.md ${RELEASE}/README.md
        find ${RELEASE}/ -name "CMakeLists.txt" -exec rm {} \;
        tar czf ${RELEASE}-${SYSTEM}.tar.gz ${RELEASE}/
        rm -rf ${RELEASE}
        ;;

    Win64)
        ;;

esac

echo "Packing complete"
