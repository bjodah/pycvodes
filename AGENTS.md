# pycvodes

This is the source code for "pycvodes"—the python wrapper of SUNDIAL's time stepping ODE-integrator "cvodes". The current task at hand (branch "sundials-7") is to add support for sundials-7.4.0.

```console
(cpython-v3.13-apt-deb) root@7f565d9c483b:/work# which python
/opt-3/cpython-v3.13-apt-deb/bin/python
(cpython-v3.13-apt-deb) root@7f565d9c483b:/work# ls /opt-3/ | grep sundials | grep release
sundials-6.7.0-release
sundials-7.4.0-release
```

## Strategy
pycvodes currently supports multiple versions of sundials by strategic use of conditional preprocessor macros which inspects the sundials (major) version.

Here's an example of how the macros are used:
```console
(cpython-v3.13-apt-deb) 17:36 root@7f565d9c483b:/work# git grep SUNDIALS_VERSION_MAJOR | wc -l
161
(cpython-v3.13-apt-deb) 17:36 root@7f565d9c483b:/work# git grep -A4 -B4 SUNDIALS_VERSION_MAJOR | head -20
pycvodes/include/cvodes_anyode.hpp-}
pycvodes/include/cvodes_anyode.hpp-
pycvodes/include/cvodes_anyode.hpp-template <class OdeSys>
pycvodes/include/cvodes_anyode.hpp-int jac_dense_cb(
pycvodes/include/cvodes_anyode.hpp:#if SUNDIALS_VERSION_MAJOR < 3
pycvodes/include/cvodes_anyode.hpp-    long int N,
pycvodes/include/cvodes_anyode.hpp-#endif
pycvodes/include/cvodes_anyode.hpp-    realtype t,
pycvodes/include/cvodes_anyode.hpp-    N_Vector y, N_Vector fy,
pycvodes/include/cvodes_anyode.hpp:#if SUNDIALS_VERSION_MAJOR < 3
pycvodes/include/cvodes_anyode.hpp-    DlsMat Jac,
pycvodes/include/cvodes_anyode.hpp-#else
pycvodes/include/cvodes_anyode.hpp-    SUNMatrix Jac,
pycvodes/include/cvodes_anyode.hpp-#endif
pycvodes/include/cvodes_anyode.hpp-    void *user_data,
pycvodes/include/cvodes_anyode.hpp-    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3
pycvodes/include/cvodes_anyode.hpp-    ){
pycvodes/include/cvodes_anyode.hpp-    // callback of req. signature wrapping OdeSys method.
pycvodes/include/cvodes_anyode.hpp:#if SUNDIALS_VERSION_MAJOR < 3
pycvodes/include/cvodes_anyode.hpp-    AnyODE::ignore(N);
```


The CI-config (.woodpecker.yaml) now contains explicit testing of sundials-7.4.0 (in addition to sundials-6.7.0), these tests need to pass before we can merge this branch to master.
```console
# bash -l .ci/test-case.sh --tmp --python /opt-3/cpython-v3.11.*-release/bin/python3 /opt-3/sundials-6.7.0-release
...
==954== For lists of detected and suppressed errors, rerun with: -s
==954== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
+ cd -
/tmp/tmp.VzXyW2X2ak
+ rm -rf -- /tmp/tmp.VzXyW2X2ak
(cpython-v3.13-apt-deb) 17:22 root@7f565d9c483b:/work# echo $?
0
(cpython-v3.13-apt-deb) 17:40 root@7f565d9c483b:/work# bash -l .ci/test-case.sh --tmp --python /opt-3/cpython-v3.11.*-release/bin/python3 /opt-3/sundials-7.4.0-release
...
running build_ext
building 'pycvodes._cvodes' extension
creating build
creating build/temp.linux-x86_64-cpython-311
creating build/temp.linux-x86_64-cpython-311/pycvodes
g++ -Wsign-compare -DNDEBUG -g -O3 -Wall -Os -march=nehalem -mtune=skylake -Os -march=nehalem -mtune=skylake -isystem /opt-3/sundials-7.4.0-release/include -D_GLIBCXX_DEBUG -D_GLIBCXX_PEDANTIC -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION -fPIC -DPYCVODES_NO_KLU=0 -DPYCVODES_NO_LAPACK=0 -DANYODE_NO_LAPACK=0 -I/opt-3/cpython-v3.11.13-release/lib/python3.11/site-packages/numpy/core/include -Ipycvodes/include -Iexternal/anyode/include -I/opt-3/cpython-v3.11.13-release/include/python3.11 -c pycvodes/_cvodes.cpp -o build/temp.linux-x86_64-cpython-311/pycvodes/_cvodes.o -DVERSION_INFO=\"0.14.7.post1+gb0bdf35.dirty\" -std=c++17
In file included from pycvodes/_cvodes.cpp:1179:
pycvodes/include/cvodes_cxx.hpp:38:10: fatal error: cvodes/cvodes_spils.h: No such file or directory
   38 | #include <cvodes/cvodes_spils.h>
      |          ^~~~~~~~~~~~~~~~~~~~~~~
compilation terminated.
error: command '/usr/bin/g++' failed with exit code 1
+ rm -rf -- /tmp/tmp.MZ6JBiEko8
```


## Background
The missing old functions have previously been deprecated by SUNDIALS and is described in `SUNDIALS-CHANGELOG-60--74.md`. 
Some snippets of relevant parts of sundials docs for 7.4.0 is found in `sundials-docs-excerpts-740/`.


## Task

Use the two example commands above for how to run the tests to guide the work in upgrading support for `sundials-7.4.0`. Ensure that old tests with sundials-6.7.0 keeps passing as we add support for version 7+.
