CXX ?= g++
LIBS ?=-lsundials_cvodes -llapack -lsundials_nvecserial
CXXFLAGS ?= -std=c++11 -Wall -Wextra -Werror -pedantic -g -ggdb -O0 -D_GLIBCXX_DEBUG
CXXFLAGS += $(EXTRA_FLAGS)
INCLUDE ?= -I../include
DEFINES ?=


.PHONY: test clean

test: test_cvodes_anyode test_cvodes_cxx test_sundials_cxx
	./test_cvodes_anyode --abortx 1
	./test_cvodes_cxx --abortx 1
	./test_sundials_cxx --abortx 1
	python ./_test_cvodes_anyode.py

clean:
	rm catch.hpp
	rm test_cvodes_cxx

test_%: test_%.cpp ../include/cvodes_cxx.hpp catch.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEFINES) -o $@ $< $(LIBS)

catch.hpp: catch.hpp.bz2
	bunzip2 -k -f $<
