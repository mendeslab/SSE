# Independent of the retained MoSSE native build. All outputs stay under SSE/build/.
CXX = g++-16
JAVA_HOME ?= $(shell dirname $(shell dirname $(shell readlink -f $(shell command -v javac))))
BUILD_DIR ?= ../../build/gcc-16-debug-O
CXXFLAGS ?= -g -Og
CPPFLAGS += -I$(BUILD_DIR)/include -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux $(shell pkg-config --cflags fftw3)
LDLIBS += $(shell pkg-config --libs fftw3)

.PHONY: all
all: $(BUILD_DIR)/libsse_quasse.so

$(BUILD_DIR)/include/SSE_FftwFFT.h: ../../src/SSE/ComplexFFT.java ../../src/SSE/FftwFFT.java
	mkdir -p $(BUILD_DIR)/include $(BUILD_DIR)/classes
	$(JAVA_HOME)/bin/javac --release 17 -h $(BUILD_DIR)/include -d $(BUILD_DIR)/classes $^

$(BUILD_DIR)/libsse_quasse.so: fft.cpp $(BUILD_DIR)/include/SSE_FftwFFT.h Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++20 -Wall -Wextra -Werror -fPIC -pthread -shared \
		$(LDFLAGS) -o $@ $< $(LDLIBS)
