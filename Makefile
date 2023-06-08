CPP=g++
CXXFLAGS=-std=c++17 -O3 -Wall -Werror -pedantic -funroll-loops -Iinclude -DNDEBUG 
DEBUGFLAGS=-std=c++17 -Og -ggdb -g -Wall -Werror -pedantic -funroll-loops -Iinclude 


SRC=$(wildcard src/*.cc)
OBJ=$(patsubst src/%.cc,bin/%.o,$(SRC))

bin/%.o: src/%.cc
	$(CPP) $(CXXFLAGS) $^ -c -o $@

sa: $(OBJ)
	$(CPP) $(CXXFLAGS) $^ -o $@

debug:CXXFLAGS=$(DEBUGFLAGS)
debug: $(OBJ)
	$(CPP) $(CXXFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -rf debug sa bin/*.o