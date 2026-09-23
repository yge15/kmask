CXX ?= g++
CXXFLAGS ?= -O3
CXXFLAGS += -std=c++17 -pthread
LDFLAGS += -pthread
PREFIX ?= /usr/local
TARGET := kmask
SRCS := main.cpp kmask.cpp progress.cpp

# GCC 8 needs -lstdc++fs for std::filesystem; newer GCC and clang do not
GCC_MAJOR := $(shell $(CXX) -dumpversion 2>/dev/null | cut -d. -f1)
IS_CLANG  := $(findstring clang,$(shell $(CXX) --version 2>/dev/null))
ifeq ($(IS_CLANG),)
ifeq ($(GCC_MAJOR),8)
LDFLAGS += -lstdc++fs
endif
endif

all:
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

install: all
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $(TARGET) $(DESTDIR)$(PREFIX)/bin/

clean:
	rm -f $(TARGET)