# Makefile for macOS

EXE       = prst
LIB_GWNUM = ../../framework/gwnum/mac64/gwnum.a

COMPOBJS_COMMON = md5.o arithmetic.o group.o giant.o lucas.o config.o inputnum.o integer.o logging.o file.o container.o task.o exp.o fermat.o order.o pocklington.o lucasmul.o morrison.o proof.o testing.o support.o batch.o
COMPOBJS        = $(COMPOBJS_COMMON) prst.o

# Source directories

VPATH     = ..:../../framework:../../framework/arithmetic

LFLAGS    = -ld_classic

DUMMY2    =
LIBS2     = -lm -lpthread -lstdc++ /usr/local/lib/libgmp.a  # link against .a file directly to avoid using .dylib

CC = clang
CXX = clang
CFLAGS    = -I/usr/local/include/ -I.. -I../../framework/gwnum                                                -DGMP -DX86_64 -O2 -fno-limit-debug-info -g -std=c99     -Wall -Wextra -target x86_64-apple-macos10.15
CXXFLAGS  = -I/usr/local/include/ -I.. -I../../framework -I../../framework/arithmetic -I../../framework/gwnum -DGMP -DX86_64 -O2 -fno-limit-debug-info -g -std=gnu++17 -Wall -Wextra -target x86_64-apple-macos10.15

# the code is too noisy for all warnings
CXXFLAGS += -Wno-sign-compare

#########################################################################

all: $(EXE)

$(EXE): $(COMPOBJS)
	$(CC) $(LFLAGS) -o $@ $(COMPOBJS) $(LIB_GWNUM) $(LIBS2)

clean::
	rm -f $(EXE) $(COMPOBJS)
