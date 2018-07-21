CXX ?= g++
LIBS ?=$(shell python3 -c "from pycvodes._config import env; print(' '.join('-l%s' % l for l in filter( lambda x: x != '', env['SUNDIALS_LIBS'].split(',') + env['LAPACK'].split(','))))")
PYCVODES_NO_LAPACK=$(shell python3 -c "from pycvodes._config import env; print('1' if env['LAPACK'] == '' else '0')")
CXXFLAGS ?= -std=c++14 -Wall -Wextra -Werror -pedantic -g -ggdb -O0 -DPYCVODES_NO_LAPACK=$(PYCVODES_NO_LAPACK)
EXTRA_FLAGS ?=-D_GLIBCXX_DEBUG
CXXFLAGS += $(EXTRA_FLAGS)
INCLUDE ?= -I../pycvodes/include -I../external/anyode/include
DEFINES ?=
OPENMP_FLAG ?= -fopenmp
OPENMP_LIB ?= -lgomp
.PHONY: test clean

TARGETS=test_cvodes_anyode test_cvodes_anyode_autorestart test_cvodes_anyode_quad test_cvodes_cxx test_sundials_cxx

ifeq ($(PYCVODES_NO_LAPACK),0)
TARGETS += test_cvodes_anyode_svd test_tough_case test_cvodes_anyode_jtimes
endif

test: $(TARGETS) test_cvodes_anyode_parallel
	./test_cvodes_anyode --abortx 1
	./test_cvodes_anyode_autorestart --abortx 1
	./test_cvodes_anyode_quad --abortx 1
ifeq ($(PYCVODES_NO_LAPACK),0)
	./test_cvodes_anyode_svd --abortx 1
	./test_tough_case --abortx 1
	./test_cvodes_anyode_jtimes --abortx 1
endif
	./test_cvodes_cxx --abortx 1
	./test_sundials_cxx --abortx 1
ifeq ($(CXX),clang++-6.0)
	@echo "skipping openmp when using clang"
else
	./test_cvodes_anyode_parallel --abortx 1
endif
	python3 ./_test_cvodes_anyode.py

clean:
	rm catch.o
	rm $(TARGETS)
	if [ -f test_cvodes_anyode_parallel ] ; then rm test_cvodes_anyode_parallel ; fi;

test_%: test_%.cpp ../pycvodes/include/cvodes_cxx.hpp catch.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEFINES) -o $@ $< catch.o $(LIBS)

test_cvodes_anyode_parallel: test_cvodes_anyode_parallel.cpp ../pycvodes/include/cvodes_*.hpp catch.o
ifeq ($(CXX),clang++-6.0)
	@echo "skipping openmp when using clang"
else
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAG) $(INCLUDE) $(DEFINES) -o $@ $< catch.o $(LIBS) $(OPENMP_LIB)
endif

catch.hpp: catch.hpp.bz2
	bzcat $< >$@

catch.o: catch.cpp catch.hpp

