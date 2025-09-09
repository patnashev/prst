# Makefile for macOS

EXE       = prst
LIB_GWNUM = ../../framework/gwnum/mac64/gwnum.a

COMPOBJS_COMMON = md5.o arithmetic.o group.o giant.o lucas.o config.o inputnum.o integer.o logging.o file.o container.o task.o exp.o fermat.o order.o pocklington.o lucasmul.o morrison.o proof.o testing.o support.o batch.o
COMPOBJS        = $(COMPOBJS_COMMON) prst.o


# libc++ headers to use. Download and extract from https://github.com/llvm/llvm-project/releases/download/llvmorg-11.0.1/libcxx-11.0.1.src.tar.xz
LIBCXX_INCLUDE ?= /opt/libcxx-11/include/c++/v1
SDKROOT = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk
MINVER  = 10.14

# some tools use this env var as well, so set it as well
export MACOSX_DEPLOYMENT_TARGET := $(MINVER)

CC  = clang
CXX = clang++

VPATH     = ..:../../framework:../../framework/arithmetic

LFLAGS    =

# uncomment if running into linker bug
LFLAGS   += -ld_classic

DUMMY2    =
LIBS2     = -lm -lpthread /usr/local/lib/libgmp.a  # static link to gmp

COMMON_CFLAGS   = -target x86_64-apple-macos$(MINVER) \
                  -mmacosx-version-min=$(MINVER) \
                  -isysroot $(SDKROOT) \
                  -I/usr/local/include/ -I.. -I../../framework -I../../framework/arithmetic -I../../framework/gwnum \
                  -DGMP -DX86_64 -O2 -g -fno-limit-debug-info -Wall -Wextra

# Force libc++ headers from LIBCXX_INCLUDE up Xcode's new ones
CFLAGS   += -std=c99
CXXFLAGS += -std=gnu++17 -stdlib=libc++ -nostdinc++ -isystem $(LIBCXX_INCLUDE)

# asan flags
#LFLAGS   += -fsanitize=address -fno-omit-frame-pointer
#CFLAGS   += -fsanitize=address -fno-omit-frame-pointer
#CXXFLAGS += -fsanitize=address -fno-omit-frame-pointer

# the code is too noisy for all warnings
CXXFLAGS += -Wno-sign-compare -Wno-unused-parameter -Wno-unused-private-field

# Linker needs the same sysroot/min-version and to prefer the SDK's libs
COMMON_LDFLAGS = -stdlib=libc++ \
                 -mmacosx-version-min=$(MINVER) \
                 -isysroot $(SDKROOT) \
                 -Wl,-syslibroot,$(SDKROOT)

all: $(EXE)

# Compile rules (inherit the common flags)
%.o: %.c
	$(CC)  $(CPPFLAGS) $(COMMON_CFLAGS) $(CFLAGS)   -c $< -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(COMMON_CFLAGS) $(CXXFLAGS) -c $< -o $@

# IMPORTANT: link with CXX (C++ linker), include COMMON_LDFLAGS, and avoid -lstdc++
$(EXE): $(COMPOBJS)
	$(CXX) $(LFLAGS) $(COMMON_LDFLAGS) -o $@ $(COMPOBJS) $(LIB_GWNUM) $(LIBS2)

clean::
	rm -f $(EXE) $(COMPOBJS)
