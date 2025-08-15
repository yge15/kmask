# Makefile for kmask
CXX := g++
CXXFLAGS := -std=c++17 -pthread
LDFLAGS := -lstdc++fs

TARGET := kmask
SRCS := main.cpp kmask.cpp

all:
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

clean:
	rm -f $(TARGET)