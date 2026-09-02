CC = gcc
JAVA_HOME ?= $(shell dirname $(shell dirname $(shell readlink -f $(shell command -v javac))))
JNI_PLATFORM ?= linux

CPPFLAGS += -I. -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/$(JNI_PLATFORM)
CFLAGS += -g -O2 -fPIC -Werror
LDFLAGS += -shared
LDLIBS += -lfftw3 -lm

DEPS = mosse-eqs-fftC.h quasse-eqs-fftC.h rfftw.h quasse_wrappers.h
OBJ = main.o quasse_wrappers.o

%.o: %.c $(DEPS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

libtest.so: $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)


.PHONY: clean

clean:
	rm -f $(OBJ) libtest.so
