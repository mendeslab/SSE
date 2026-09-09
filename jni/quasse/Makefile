# Independent of the retained MoSSE native build. All outputs stay under SSE/build/.
# Produces libsse_quasse.so for SSE.fft.FftwFFT; Ant's optional native target invokes this file.
CXX = g++-16
JAVA_HOME ?= $(shell dirname $(shell dirname $(shell readlink -f $(shell command -v javac))))
BUILD_DIR ?= ../../build/gcc-16-debug-O
CXXFLAGS ?= -g -Og
CPPFLAGS += -I$(BUILD_DIR)/include -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux $(shell pkg-config --cflags fftw3)
LDLIBS += $(shell pkg-config --libs fftw3)

.PHONY: all
all: $(BUILD_DIR)/libsse_quasse.so

# Generate signatures from Java. These classes only support header generation;
# Ant/IDE compilation supplies application classes, so this directory is not on its classpath.
$(BUILD_DIR)/include/SSE_fft_FftwFFT.h: ../../src/SSE/fft/ComplexFFT.java ../../src/SSE/fft/FftwFFT.java
	mkdir -p $(BUILD_DIR)/include $(BUILD_DIR)/classes
	$(JAVA_HOME)/bin/javac --release 17 -h $(BUILD_DIR)/include -d $(BUILD_DIR)/classes $^

# Link only QuaSSE's C++ code and FFTW; no MoSSE or GSL dependency.
$(BUILD_DIR)/libsse_quasse.so: fft.cpp $(BUILD_DIR)/include/SSE_fft_FftwFFT.h Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++20 -Wall -Wextra -Werror -fPIC -pthread -shared \
		$(LDFLAGS) -o $@ $< $(LDLIBS)
