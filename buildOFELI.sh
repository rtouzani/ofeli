#! /bin/sh
# buildOFELI.sh 
# A script build OFELI distribution

if [ $# -ne 2 ]; then
   echo 1>&2 "Usage: buildOFELI.sh <release> <MacOSX|Linux64|Linux32|Win32|Win64>"
   exit 2
fi

RELEASE=$1
REL=$1
PREFIX="/usr/local"
SYSTEM=$2

if [${SYSTEM} != "MacOSX"] &&
   [${SYSTEM} != "Linux64"] &&
   [${SYSTEM} != "Linux32"] &&
   [${SYSTEM} != "Win32"] &&
   [${SYSTEM} != "Win64"]; then
     echo "Error: Unavailable Platform $MACHINE"
   exit
fi

echo "Preparing release $RELEASE..."
echo "--------------------------"

RELEASE=ofeli-${RELEASE}


case $SYSTEM in

    MacOSX)
        echo "Creating Source Package ..."
        tar --exclude='./build'  -czf ofeli.tar.gz ofeli
        mv ofeli ${RELEASE}
        tar --exclude='./build' --exclude='.git' --exclude='.vscode' --exclude='./buildOFELI.sh' --exclude='install.sh' -czf ${RELEASE}-src.tar.gz ${RELEASE}
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

        echo "Creating Software Development Kit Package ..."
        mkdir -p ${RELEASE}-${SYSTEM}
        mkdir -p ${RELEASE}-${SYSTEM}/bin
        mkdir -p ${RELEASE}-${SYSTEM}/lib
        mkdir -p ${RELEASE}-${SYSTEM}/include
        mkdir -p ${RELEASE}-${SYSTEM}/material
        cp ${PREFIX}/lib/libofeli.a ${RELEASE}-${SYSTEM}/lib/.
        cp ${PREFIX}/bin/cmesh ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/bin/cfield ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/bin/g2m ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/bin/vmesh ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/share/ofeli/material/* ${RELEASE}-${SYSTEM}/material/.
        cp -rf ${PREFIX}/include/ofeli/* ${RELEASE}-${SYSTEM}/include/.
        cp ofeli/install.sh ${RELEASE}-${SYSTEM}/install.sh
        cp ofeli/README-sdk.md ${RELEASE}-${SYSTEM}/README.md
        find ${RELEASE}-${SYSTEM}/ -name "CMakeLists.txt" -exec rm {} \;
        tar czf ${RELEASE}-${SYSTEM}.tar.gz ${RELEASE}-${SYSTEM}/
        rm -rf ${RELEASE}-${SYSTEM}
        ;;

    Linux64)
        echo "Creating Software Development Kit Package ..."
        mkdir -p ${RELEASE}-${SYSTEM}
        mkdir -p ${RELEASE}-${SYSTEM}/bin
        mkdir -p ${RELEASE}-${SYSTEM}/lib
        mkdir -p ${RELEASE}-${SYSTEM}/include
        mkdir -p ${RELEASE}-${SYSTEM}/material
        cp ${PREFIX}/lib/libofeli.a ${RELEASE}-${SYSTEM}/lib/.
        cp ${PREFIX}/bin/cmesh ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/bin/cfield ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/bin/g2m ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/bin/vmesh ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/share/ofeli/material/* ${RELEASE}-${SYSTEM}/material/.
        cp -rf ${PREFIX}/include/ofeli/* ${RELEASE}-${SYSTEM}/include/.
        cp ofeli/install.sh ${RELEASE}-${SYSTEM}/install.sh
        cp ofeli/README-sdk.md ${RELEASE}-${SYSTEM}/README.md
        find ${RELEASE}-${SYSTEM}/ -name "CMakeLists.txt" -exec rm {} \;
        tar czf ${RELEASE}-${SYSTEM}.tar.gz ${RELEASE}-${SYSTEM}/
        rm -rf ${RELEASE}-${SYSTEM}
        ;;

    Linux32)
        echo "Creating Software Development Kit Package ..."
        mkdir -p ${RELEASE}-${SYSTEM}
        mkdir -p ${RELEASE}-${SYSTEM}/bin
        mkdir -p ${RELEASE}-${SYSTEM}/lib
        mkdir -p ${RELEASE}-${SYSTEM}/include
        mkdir -p ${RELEASE}-${SYSTEM}/material
        cp ${PREFIX}/lib/libofeli.a ${RELEASE}-${SYSTEM}/lib/.
        cp ${PREFIX}/bin/cmesh ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/bin/cfield ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/bin/g2m ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/bin/vmesh ${RELEASE}-${SYSTEM}/bin/.
        cp ${PREFIX}/share/ofeli/material/* ${RELEASE}-${SYSTEM}/material/.
        cp -rf ${PREFIX}/include/ofeli/* ${RELEASE}-${SYSTEM}/include/.
        cp ofeli/install.sh ${RELEASE}-${SYSTEM}/install.sh
        cp ofeli/README-sdk.md ${RELEASE}-${SYSTEM}/README.md
        find ${RELEASE}-${SYSTEM}/ -name "CMakeLists.txt" -exec rm {} \;
        tar czf ${RELEASE}-${SYSTEM}.tar.gz ${RELEASE}-${SYSTEM}/
        rm -rf ${RELEASE}-${SYSTEM}
        ;;

    Win64)
        ;;

    Win32)
        ;;

esac

echo "Packing complete"
