CXX       ?= clang++
CXXFLAGS  ?= -std=c++20 -g -O3
CPPFLAGS  ?= -I. 
LDFLAGS   ?= -L${mkOpenblasLib}
LIBS      ?= -lopenblas
OPTS      ?=

all : main

main : main.o 
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

main.o : main.cpp Matrix.hpp Matrix_impl.hpp
	$(CXX) $(OPTS) $(CXXFLAGS) $(CPPFLAGS) -c $< 

clean :
	$(RM) *.o

distclean : clean
	$(RM) main


