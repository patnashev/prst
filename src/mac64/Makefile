# Makefile for Linux

EXE       = prst
LIB_GWNUM = ../../framework/gwnum/mac64/gwnum.a

COMPOBJS_COMMON = md5.o arithmetic.o group.o giant.o lucas.o config.o inputnum.o integer.o logging.o file.o task.o exp.o fermat.o pocklington.o morrison.o proof.o testing.o support.o
COMPOBJS        = $(COMPOBJS_COMMON) prst.o

# Source directories

VPATH     = ..:../../framework:../../framework/arithmetic

# See README.md for how to get macOS SDK
LFLAGS    = # -isysroot MacOSX10.7.sdk/ -mmacosx-version-min=10.7

DUMMY2    =
LIBS2     = -lm -lpthread -lstdc++ /usr/local/lib/libgmp.a  # link against .a file directly to avoid using .dylib

#CC        = /usr/local/bin/gcc-12 -mmacosx-version-min=10.9
#CC        = /usr/local/bin/gcc-12
#CXX       = /usr/local/bin/g++-12
CC = gcc
CXX = g++
CFLAGS    = -I/usr/local/include/ -I.. -I../../framework/gwnum                                                -DGMP -DX86_64 -O2 -std=c99     -Wall -Wextra
CXXFLAGS  = -I/usr/local/include/ -I.. -I../../framework -I../../framework/arithmetic -I../../framework/gwnum -DGMP -DX86_64 -O2 -std=gnu++17 -Wall -Wextra

# the code is too noisy for all warnings
CXXFLAGS += -Wno-sign-compare

#########################################################################

all: $(EXE)

$(EXE): $(COMPOBJS)
	$(CC) $(LFLAGS) -o $@ $(COMPOBJS) $(LIB_GWNUM) $(LIBS2)

clean:
	rm -f $(EXE) $(COMPOBJS)
