# Makefile for Linux

EXE       = prst
LIB_GWNUM = /Users/daohe/CLionProjects/prst/framework/gwnum/mac64/__newbuild/gwnum/gwnum.a

COMPOBJS_COMMON = md5.o arithmetic.o group.o giant.o lucas.o config.o inputnum.o integer.o logging.o file.o task.o exp.o fermat.o pocklington.o morrison.o proof.o testing.o support.o
COMPOBJS        = $(COMPOBJS_COMMON) prst.o

# Source directories

VPATH    = ..:../../framework:../../framework/arithmetic

LIBS2    = -lm -lpthread -lstdc++
#CC       = /usr/local/bin/gcc-12 -isysroot MacOSX10.9.sdk -mmacosx-version-min=10.9 -Wl, -rpath, /usr/local/opt/gcc@10/lib/gcc/10/
CC       = /usr/local/bin/gcc-12 -mmacosx-version-min=10.9
#CXX      = /usr/local/bin/g++-12 -isysroot MacOSX10.9.sdk -mmacosx-version-min=10.9 -Wl, -rpath, /usr/local/opt/gcc@10/lib/gcc/10/
CXX      = /usr/local/bin/g++-12 -mmacosx-version-min=10.9
LDFLAGS  = -isysroot MacOSX10.9.sdk/
CFLAGS   = -I/usr/local/include/ -I.. -I../../framework/gwnum                                                -DX86_64 -O0 -g -std=c99     -Wall -Wextra
CXXFLAGS = -I/usr/local/include/ -I.. -I../../framework -I../../framework/arithmetic -I../../framework/gwnum -DX86_64 -O0 -g -std=gnu++17 -Wall -Wextra
CXXFLAGS += -Wno-sign-compare

#########################################################################

all: $(EXE)

$(EXE): $(COMPOBJS)
	$(CC) $(LDFLAGS) -o $@ $(COMPOBJS) $(LIB_GWNUM) $(LIBS2)

clean:
	rm -f $(EXE) $(COMPOBJS)
