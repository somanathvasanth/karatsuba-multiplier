SHELL := /bin/bash

CXX = g++
CXXFLAGS =-std=c++17
EXEC = helper
SRC = 23b0924.cpp


.PHONY: all test

all: $(EXEC) test

$(EXEC): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(SRC)

test: $(EXEC)
	@for i in $(shell seq 1 10); do \
	    echo "Running test $$i..."; \
	    time ./$(EXEC) < input$$i > temp_output$$i; \
	    if diff temp_output$$i output$$i; then \
	        echo -e "\e[32mTest $$i: PASSED\e[0m"; \
	    else \
	        echo -e "\e[31mTest $$i: FAILED\e[0m"; \
	    fi; \
	    echo ""; \
	done

clean:
	rm -f $(EXEC) temp_output*
