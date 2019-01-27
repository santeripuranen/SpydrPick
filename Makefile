CXX = g++

BUILDDIR = build
SRCDIR = src

CXXFLAGS = -std=c++11 -pthread -O3 -Wall -pedantic

all: spydrpick 

spydrpick: 
		mkdir -p $(BUILDDIR)
			$(CXX) $(CXXFLAGS) $(SRCDIR)/ARACNE.cpp -o $(BUILDDIR)/aracne
			$(CXX) $(CXXFLAGS) $(SRCDIR)/mi_gwes.cpp -o $(BUILDDIR)/mi_gwes

.PHONY: clean

clean:
	\rm $(BUILDDIR)/aracne $(BUILDDIR)/mi_gwes
