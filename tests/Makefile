CXX ?= g++
PYTHON ?= python3
PYCVODES_NO_LAPACK ?= 0
PYCVODES_NO_KLU ?= 0
CXXFLAGS ?= -std=c++14 -Wall -Wextra -Werror -pedantic -g -ggdb -O0 -DANYODE_NO_LAPACK=$(PYCVODES_NO_LAPACK) -DPYCVODES_NO_LAPACK=$(PYCVODES_NO_LAPACK) -DPYCVODES_NO_KLU=$(PYCVODES_NO_KLU)
EXTRA_FLAGS ?=-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
CXXFLAGS += $(EXTRA_FLAGS) $(CFLAGS)
INCLUDE ?= -I../pycvodes/include -I../external/anyode/include $(shell $(PYTHON)-config --includes)
DEFINES ?=
OPENMP_FLAG ?= -fopenmp
OPENMP_LIB ?= -lgomp
.PHONY: test clean

TARGETS=test_cvodes_anyode test_cvodes_anyode_autorestart test_cvodes_anyode_quad test_cvodes_cxx test_sundials_cxx

ifneq ($(PYCVODES_NO_LAPACK),1)
TARGETS += test_cvodes_anyode_svd test_tough_case test_cvodes_anyode_jtimes
endif
ifneq ($(PYCVODES_NO_KLU), 1)
TARGETS += test_cvodes_anyode_sparse
endif

test: $(TARGETS) test_cvodes_anyode_parallel
	./test_cvodes_anyode --abortx 1
	./test_cvodes_anyode_autorestart --abortx 1
	./test_cvodes_anyode_quad --abortx 1
ifneq ($(PYCVODES_NO_LAPACK),1)
	./test_cvodes_anyode_svd --abortx 1
	./test_tough_case --abortx 1
	./test_cvodes_anyode_jtimes --abortx 1
endif
ifneq ($(PYCVODES_NO_KLU), 1)
	./test_cvodes_anyode_sparse --abortx 1
endif
	./test_cvodes_cxx --abortx 1
	./test_sundials_cxx --abortx 1
ifeq ($(CXX),clang++-11)
	@echo "skipping openmp when using clang"
else
	./test_cvodes_anyode_parallel --abortx 1
endif
	$(PYTHON) ./_test_cvodes_anyode.py

clean:
	rm catch.o
	rm $(TARGETS)
	if [ -f test_cvodes_anyode_parallel ] ; then rm test_cvodes_anyode_parallel ; fi;

test_%: test_%.cpp ../pycvodes/include/cvodes_cxx.hpp catch.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEFINES) -o $@ $< catch.o $(LDFLAGS)

test_cvodes_anyode_parallel: test_cvodes_anyode_parallel.cpp ../pycvodes/include/cvodes_*.hpp catch.o
ifeq ($(CXX),clang++-6.0)
	@echo "skipping openmp when using clang"
else
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAG) $(INCLUDE) $(DEFINES) -o $@ $< catch.o $(LDFLAGS) $(OPENMP_LIB)
endif

catch.hpp: catch.hpp.bz2
	bzcat $< >$@

catch.o: catch.cpp catch.hpp

