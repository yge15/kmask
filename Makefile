# Makefile for kmask
CXX := g++
CXXFLAGS := -std=c++17 -pthread -O3 -march=core-avx2
LDFLAGS := -lstdc++fs

TARGET := kmask
SRCS := main.cpp kmask.cpp progress.cpp

all:
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

clean:
	rm -f $(TARGET)
