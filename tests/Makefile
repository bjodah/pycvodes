CXX ?= g++
LAPACK ?= -llapack -lblas
LIBS ?=-lsundials_cvodes -lsundials_nvecserial $(LAPACK)
LIBS += $(EXTRA_LIBS)
CXXFLAGS ?= -std=c++14 -Wall -Wextra -Werror -pedantic -g -ggdb -O0
# -O3 -ffast-math -funroll-loops -march=native
EXTRA_FLAGS ?=-D_GLIBCXX_DEBUG
CXXFLAGS += $(EXTRA_FLAGS)
INCLUDE ?= -I../pycvodes/include -I../external/anyode/include
DEFINES ?=
OPENMP_FLAG ?= -fopenmp
OPENMP_LIB ?= -lgomp

.PHONY: test clean

TARGETS=test_cvodes_anyode test_cvodes_anyode_autorestart test_cvodes_anyode_svd test_cvodes_cxx test_sundials_cxx test_tough_case

test: test_cvodes_anyode test_cvodes_anyode_autorestart test_cvodes_anyode_svd test_cvodes_cxx test_sundials_cxx test_tough_case test_cvodes_anyode_parallel
	./test_cvodes_anyode --abortx 1
	./test_cvodes_anyode_autorestart --abortx 1
	./test_cvodes_anyode_svd --abortx 1
	./test_cvodes_cxx --abortx 1
	./test_sundials_cxx --abortx 1
	./test_tough_case --abortx 1
ifeq ($(CXX),clang++-5.0)
	@echo "skipping openmp when using clang"
else
	./test_cvodes_anyode_parallel --abortx 1
endif
	python ./_test_cvodes_anyode.py

clean:
	rm catch.hpp
	rm $(TARGETS)
	if [ -f test_cvodes_anyode_parallel ] ; then rm test_cvodes_anyode_parallel ; fi;

test_%: test_%.cpp ../pycvodes/include/cvodes_cxx.hpp catch.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEFINES) -o $@ $< $(LIBS)

test_cvodes_anyode_parallel: test_cvodes_anyode_parallel.cpp ../pycvodes/include/cvodes_*.hpp catch.hpp
ifeq ($(CXX),clang++-5.0)
	@echo "skipping openmp when using clang"
else
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAG) $(INCLUDE) $(DEFINES) -o $@ $< $(LIBS) $(OPENMP_LIB)
endif

catch.hpp: catch.hpp.bz2
	bunzip2 -k -f $<
