CXX ?= g++
PKG_CONFIG ?= pkg-config
RNALIB_PKG ?= RNAlib2

CXXFLAGS ?= -O3 -std=c++0x
CPPFLAGS ?= -Iinclude
LDFLAGS ?=
LDLIBS ?= -lpthread
OPENMP ?= -fopenmp

# If pkg-config cannot find RNAlib2, either set PKG_CONFIG_PATH or override:
#   make BPmap_test RNALIB_CFLAGS="-I/path/to/ViennaRNA/include" RNALIB_LIBS="-L/path/to/lib -lRNA"
RNALIB_CFLAGS ?= $(shell $(PKG_CONFIG) --cflags $(RNALIB_PKG) 2>/dev/null)
RNALIB_LIBS ?= $(shell $(PKG_CONFIG) --libs $(RNALIB_PKG) 2>/dev/null)

COMMON_SOURCES = src/load_data.cpp src/buildmap.cpp src/msa.cpp
CORE_SOURCES = $(COMMON_SOURCES) src/align.cpp

.PHONY: all BPmap BPmap_test clean check-rnalib

all: BPmap BPmap_test

check-rnalib:
	@if [ -z "$(RNALIB_CFLAGS)$(RNALIB_LIBS)" ]; then \
		echo "ERROR: cannot find $(RNALIB_PKG) via pkg-config."; \
		echo "Try one of:"; \
		echo "  export PKG_CONFIG_PATH=/path/to/ViennaRNA/lib/pkgconfig:\$$PKG_CONFIG_PATH"; \
		echo "  make BPmap_test RNALIB_CFLAGS='-I/path/to/include' RNALIB_LIBS='-L/path/to/lib -lRNA'"; \
		exit 1; \
	fi

bin:
	mkdir -p bin

artifacts:
	mkdir -p artifacts

BPmap: check-rnalib bin
	$(CXX) $(OPENMP) $(CXXFLAGS) $(CPPFLAGS) $(RNALIB_CFLAGS) \
		$(CORE_SOURCES) apps/BPmap.cpp \
		$(LDFLAGS) $(RNALIB_LIBS) $(LDLIBS) \
		-o bin/BPmap

BPmap_test: check-rnalib artifacts
	$(CXX) $(OPENMP) $(CXXFLAGS) $(CPPFLAGS) $(RNALIB_CFLAGS) \
		$(CORE_SOURCES) apps/BPmap.cpp \
		$(LDFLAGS) $(RNALIB_LIBS) $(LDLIBS) \
		-o artifacts/BPmap_test

clean:
	rm -f bin/BPmap artifacts/BPmap_test
