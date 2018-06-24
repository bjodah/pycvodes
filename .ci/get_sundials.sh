#!/bin/bash -eu
VERSION=$1
PREFIX=$2
if [ -d "$PREFIX" ]; then 2>&1 echo "Directory already exists: $PREFIX"; exit 1; fi
if [[ $VERSION == "3.1.1" ]]; then
    SUNDIALS_FNAME="sundials-3.1.1.tar.gz"
    SUNDIALS_MD5="e63f4de0be5be97f750b30b0fa11ef34"
elif [[ $VERSION == "2.7.0" ]]; then
    SUNDIALS_FNAME="sundials-2.7.0.tar.gz"
    SUNDIALS_MD5="c304631b9bc82877d7b0e9f4d4fd94d3"    
else
    2>&1 echo "Unknown sundials version $VERSION"
fi

SUNDIALS_URLS=(\
    "http://hera.physchem.kth.se/~repo/${SUNDIALS_MD5}/${SUNDIALS_FNAME}" \
    "http://computation.llnl.gov/projects/sundials/download/${SUNDIALS_FNAME}" \
)
TIMEOUT=60  # 60 seconds

for URL in "${SUNDIALS_URLS[@]}"; do
    if echo $SUNDIALS_MD5 $SUNDIALS_FNAME | md5sum -c --; then
        echo "Found ${SUNDIALS_FNAME} with matching checksum, using this file."
    else
        echo "Downloading ${URL}..."
        timeout $TIMEOUT wget --quiet --tries=2 --timeout=$TIMEOUT $URL -O $SUNDIALS_FNAME || continue
    fi
    if echo $SUNDIALS_MD5 $SUNDIALS_FNAME | md5sum -c --; then
        tar xzf $SUNDIALS_FNAME
        mkdir sundials_build
        cd sundials_build
        cmake -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX \
              -DCMAKE_BUILD_TYPE:STRING="Release" \
              -DBUILD_SHARED_LIBS:BOOL=ON \
              -DBUILD_STATIC_LIBS:BOOL=OFF \
              -DEXAMPLES_ENABLE_C:BOOL=OFF \
              -DEXAMPLES_INSTALL:BOOL=OFF \
              -DOPENMP_ENABLE:BOOL=OFF \
              "${@:3}" ../sundials-${VERSION}/
        make -j 2 >/dev/null 2>&1
        make install
        if [ $? -ne 0 ]; then
            2>&1 echo "Build/install of sundials-${VERSION} failed."
            exit 1
        fi
        cd ..
        rm -r sundials*
        exit 0
    fi
done
exit 1
