# Independent of the retained MoSSE native build. All outputs stay under SSE/build/.
# Produces libsse_quasse.so for FFTs and integration; Ant's optional native target invokes it.
MODE ?= release
ifeq ($(MODE),release)
BUILD_DIR ?= ../../build/native
CXXFLAGS ?= -O3
else ifeq ($(MODE),debug-O)
BUILD_DIR ?= ../../build/native-debug-O
CXXFLAGS ?= -O -g
else ifeq ($(MODE),debug)
BUILD_DIR ?= ../../build/native-debug
CXXFLAGS ?= -g
else
$(error Unknown MODE '$(MODE)'; choose release, debug-O, or debug)
endif

# Use Make's default C++ compiler, or CXX supplied by the caller. Architecture flags are opt-in.
JAVA_HOME ?= $(shell dirname $(shell dirname $(shell readlink -f $(shell command -v javac))))

.PHONY: all clean
all: $(BUILD_DIR)/libsse_quasse.so

# Delete only this build's known outputs, preserving other modes and unrelated files.
# The obsolete configuration records can be dropped once all existing builds have been cleaned.
clean:
	rm -f -- "$(BUILD_DIR)/libsse_quasse.so" \
	    "$(BUILD_DIR)/include/SSE_fft_FftwFFT.h" "$(BUILD_DIR)/include/SSE_QuaSSENativeIntegrator.h" \
	    "$(BUILD_DIR)/classes/SSE/fft/ComplexFFT.class" "$(BUILD_DIR)/classes/SSE/fft/FftwFFT.class" \
	    '$(BUILD_DIR)/classes/SSE/fft/FftwFFT$$NativeState.class' \
	    "$(BUILD_DIR)/classes/SSE/QuaSSENativeIntegrator.class" \
	    '$(BUILD_DIR)/classes/SSE/QuaSSENativeIntegrator$$NativeState.class' \
	    "$(BUILD_DIR)/build-config" "$(BUILD_DIR)/build-config.tmp"

# Generate signatures from Java. These classes only support header generation;
# Ant/IDE compilation supplies application classes, so this directory is not on its classpath.
$(BUILD_DIR)/include/SSE_fft_FftwFFT.h: ../../src/SSE/fft/ComplexFFT.java ../../src/SSE/fft/FftwFFT.java
	mkdir -p $(BUILD_DIR)/include $(BUILD_DIR)/classes
	"$(JAVA_HOME)/bin/javac" --release 17 -h $(BUILD_DIR)/include -d $(BUILD_DIR)/classes $^

# Link only QuaSSE's C++ code and FFTW; no MoSSE or GSL dependency.
$(BUILD_DIR)/include/SSE_QuaSSENativeIntegrator.h: ../../src/SSE/QuaSSENativeIntegrator.java
	mkdir -p $(BUILD_DIR)/include $(BUILD_DIR)/classes
	"$(JAVA_HOME)/bin/javac" --release 17 -h $(BUILD_DIR)/include -d $(BUILD_DIR)/classes $<

$(BUILD_DIR)/libsse_quasse.so: fft.cpp integration.cpp native_support.h \
        $(BUILD_DIR)/include/SSE_fft_FftwFFT.h $(BUILD_DIR)/include/SSE_QuaSSENativeIntegrator.h \
        Makefile
	$(CXX) $(CPPFLAGS) -I$(BUILD_DIR)/include \
	    -I"$(JAVA_HOME)/include" -I"$(JAVA_HOME)/include/linux" $(shell pkg-config --cflags fftw3) \
	    $(CXXFLAGS) -std=c++20 -Wall -Wextra -Werror -fPIC -pthread -shared $(LDFLAGS) \
	    -o $@ fft.cpp integration.cpp $(LDLIBS) $(shell pkg-config --libs fftw3)
