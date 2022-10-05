# Makefile for Linux

EXE       = prst
LIB_GWNUM = ../../framework/gwnum/linux64/gwnum.a

COMPOBJS_COMMON = md5.o arithmetic.o group.o giant.o lucas.o inputnum.o integer.o logging.o file.o task.o exp.o fermat.o pocklington.o proof.o testing.o support.o
COMPOBJS   = $(COMPOBJS_COMMON) prst.o

# Source directories

VPATH     = ..:../../framework:../../framework/arithmetic

# insert proper dummy file to get _GWDATA data segment aligned
# on a 32 byte boundary.  Use proper dummyNN.o to do the alignment

DUMMY2 = ../../framework/gwnum/linux64/gwnum.ld
LIBS2  = -lm -lpthread -lstdc++ -static

CC  = gcc
CXX = g++
CFLAGS   = -I.. -I../../framework/gwnum                                                -DX86_64 -O2 -std=c99     -Wall -Wextra
CXXFLAGS = -I.. -I../../framework -I../../framework/arithmetic -I../../framework/gwnum -DX86_64 -O2 -std=gnu++17 -Wall -Wextra

#########################################################################

all: $(EXE)

$(EXE): $(COMPOBJS)
	$(CC) -o $@ $(COMPOBJS) $(LIB_GWNUM) $(DUMMY2) $(LIBS2)

clean:
	rm -f $(EXE) $(COMPOBJS)
