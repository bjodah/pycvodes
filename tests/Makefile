CXX ?= g++
PYTHON ?= python3
PYCVODES_NO_LAPACK ?= 0
PYCVODES_NO_KLU ?= 0
CXXFLAGS := -std=c++20 -Wall -Wextra -Werror -pedantic -g -O0 -DANYODE_NO_LAPACK=$(PYCVODES_NO_LAPACK) -DPYCVODES_NO_LAPACK=$(PYCVODES_NO_LAPACK) -DPYCVODES_NO_KLU=$(PYCVODES_NO_KLU) $(CXXFLAGS)
ifneq (,$(findstring clang,$(CXX)))
	# Found "clang" in CXX
	CXXFLAGS += -glldb
	EXTRA_FLAGS ?=
else
	CXXFLAGS += -ggdb3
	EXTRA_FLAGS ?=-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
endif

CXXFLAGS += $(EXTRA_FLAGS) $(CFLAGS)
#INCLUDE ?= -I../pycvodes/include -I../external/anyode/include $(shell $(PYTHON)-config --includes)
INCLUDE ?= -I../pycvodes/include -I../external/anyode/include $(shell $(PYTHON) -c "import sysconfig; print('-I' + sysconfig.get_path('include'))")
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
	env CFLAGS="$(CXXFLAGS)" $(PYTHON) ./_test_cvodes_anyode.py
	$(CONTEXT) ./test_cvodes_anyode --abortx 1
	$(CONTEXT) ./test_cvodes_anyode_autorestart --abortx 1
	$(CONTEXT) ./test_cvodes_anyode_quad --abortx 1
ifneq ($(PYCVODES_NO_LAPACK),1)
	$(CONTEXT) ./test_cvodes_anyode_svd --abortx 1
	$(CONTEXT) ./test_tough_case --abortx 1
	$(CONTEXT) ./test_cvodes_anyode_jtimes --abortx 1
endif
ifneq ($(PYCVODES_NO_KLU), 1)
	$(CONTEXT) ./test_cvodes_anyode_sparse --abortx 1
endif
	$(CONTEXT) ./test_cvodes_cxx --abortx 1
	$(CONTEXT) ./test_sundials_cxx --abortx 1
ifeq ($(CXX),clang++-11)
	@echo "skipping openmp when using clang"
else
	$(CONTEXT) ./test_cvodes_anyode_parallel --abortx 1
endif

clean:
	$(RM) doctest.o
	$(RM) $(TARGETS)
	if [ -f test_cvodes_anyode_parallel ] ; then rm test_cvodes_anyode_parallel ; fi;

test_%: test_%.cpp ../pycvodes/include/cvodes_cxx.hpp doctest.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEFINES) -o $@ $< doctest.o $(LDFLAGS)

test_cvodes_anyode_parallel: test_cvodes_anyode_parallel.cpp ../pycvodes/include/cvodes_*.hpp doctest.o
ifeq ($(CXX),clang++-6.0)
	@echo "skipping openmp when using clang"
else
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAG) $(INCLUDE) $(DEFINES) -o $@ $< doctest.o $(LDFLAGS) $(OPENMP_LIB)
endif

doctest.h: doctest.h.bz2
	bzcat $< >$@

doctest.o: doctest.cpp doctest.h

