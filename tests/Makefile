CXX ?= g++
LIBS ?=-lsundials_cvodes -lsundials_nvecserial -llapack -lblas
CXXFLAGS ?= -std=c++11 -Wall -Wextra -Werror -pedantic -g -ggdb -O0 -D_GLIBCXX_DEBUG
# -O3 -ffast-math -funroll-loops -march=native
CXXFLAGS += $(EXTRA_FLAGS)
INCLUDE ?= -I../include -I../external/anyode/include
DEFINES ?=
OPENMP_FLAG ?= -fopenmp
OPENMP_LIB ?= -lgomp

.PHONY: test clean

test: test_cvodes_anyode test_cvodes_anyode_autorestart test_cvodes_anyode_svd test_cvodes_cxx test_sundials_cxx test_cvodes_anyode_parallel test_tough_case
	./test_tough_case --abortx 1
	./test_cvodes_anyode --abortx 1
	./test_cvodes_anyode_autorestart --abortx 1
	./test_cvodes_anyode_svd --abortx 1
	./test_cvodes_cxx --abortx 1
	./test_sundials_cxx --abortx 1
	./test_cvodes_anyode_parallel --abortx 1
	python ./_test_cvodes_anyode.py

clean:
	rm catch.hpp
	rm test_cvodes_cxx

test_%: test_%.cpp ../include/cvodes_cxx.hpp catch.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEFINES) -o $@ $< $(LIBS)

test_cvodes_anyode_parallel: test_cvodes_anyode_parallel.cpp ../include/cvodes_*.hpp
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAG) $(INCLUDE) $(DEFINES) -o $@ $< $(LIBS) $(OPENMP_LIB)

catch.hpp: catch.hpp.bz2
	bunzip2 -k -f $<
