#!/bin/bash
TIMEOUT=60  # 60 seconds
SUNDIALS_FNAME="sundials-3.1.0.tar.gz"
SUNDIALS_MD5="1a84ca41c7f71067e03d519ddbcd9dae"
SUNDIALS_URLS=(\
"http://hera.physchem.kth.se/~repo/${SUNDIALS_MD5}/${SUNDIALS_FNAME}" \
"http://computation.llnl.gov/projects/sundials/download/${SUNDIALS_FNAME}" \
)
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
        cmake -DCMAKE_BUILD_TYPE=Release \
              -DBUILD_SHARED_LIBS=ON \
              -DBUILD_STATIC_LIBS=OFF \
              -DEXAMPLES_ENABLE_C=OFF \
              -DEXAMPLES_INSTALL=OFF \
              -DOPENMP_ENABLE=OFF \
              -DLAPACK_ENABLE=ON \
              -DSUNDIALS_INDEX_TYPE=int32_t \
              ../sundials-*/
        make install
        cd ..
        rm -r sundials*
        exit 0
    fi
done
exit 1
