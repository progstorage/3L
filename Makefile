if(APPLE)
CC=clang++
CFLAGS= -Xpreprocessor -fopenmp -lomp -std=c++17
else()
CC=g++
CFLAGS= -fopenmp -std=c++17 -lstdc++fs
endif()
SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
all:
$(CC) $(CFLAGS) $(SOURCES) $(OBJECTS)
