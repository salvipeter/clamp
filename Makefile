all: clamp.o clamp-test

LIBGEOM=../libgeom

# CXXFLAGS=-std=c++20 -Wall -pedantic -I$(LIBGEOM) -O3
CXXFLAGS=-std=c++20 -Wall -pedantic -I$(LIBGEOM) -O0 -g -fsanitize=address
LDFLAGS=-L$(LIBGEOM)/release -lgeom

clamp.o: clamp.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clamp-test: clamp-test.cc clamp.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)
