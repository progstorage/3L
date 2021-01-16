ifeq ($(shell uname),Linux)
	CC=g++
	CFLAGS= -fopenmp -std=c++17 -lstdc++fs
else ifeq ($(shell uname),Darwin)
	CC=clang++
	CFLAGS= -Xpreprocessor -fopenmp -lomp -std=c++17
endif
SOURCES=main.cpp
all:
	$(CC) $(CFLAGS) $(SOURCES)
