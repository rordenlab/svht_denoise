# svht_denoise -- build.
#
# Targets
#   make             optimised build (default)
#   make debug       -O0 -g
#   make ubsan       UndefinedBehaviorSanitizer
#   make verbose     -Weverything, for exploration; NOT clean, and not a gate
#   make clean
#   make distclean   also sweep .dSYM bundles
#
# The only dependency beyond libm and pthreads is zlib, for .nii.gz via nifti_io.
#
# There is deliberately NO -ffast-math anywhere.  This code validates finite
# values and promises byte-identical output across thread counts; reassociation
# and non-finite shortcuts are incompatible with both.

CC      ?= cc
OPT     ?= -O3
STD      = -std=c99
WARN     = -Wall -Wextra -Wshadow -Wstrict-prototypes -Wmissing-prototypes -Wpointer-arith
# `override` so these survive a command-line assignment.  Without it,
# `make DEFS=-DNOTHING` silently drops the feature macros and HAVE_ZLIB, which
# changes behaviour rather than performance -- and the previous comment claimed
# that could not happen, which was simply untrue.
override DEFS = -D_POSIX_C_SOURCE=200809L
# Darwin hides the BSD types that <sys/sysctl.h> needs unless _DARWIN_C_SOURCE is
# also set; without it a strict _POSIX_C_SOURCE build cannot ask for the core count.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  override DEFS += -D_DARWIN_C_SOURCE
endif
# CFLAGS is yours to override (optimisation, extra warnings).  Everything the
# code REQUIRES in order to be correct lives in DEFS and REQ_CFLAGS, both marked
# `override`, so no command-line assignment can drop the language standard, the
# feature macros or -DHAVE_ZLIB -- each of which changes behaviour rather than
# just performance.
# User CPPFLAGS is appended after DEFS, so it can add but not remove.
#
# This is a plain `=`, NOT `?=`, and that is load-bearing.  make imports
# environment variables as already-defined, so with `?=` a developer or CI with
# CFLAGS in their shell profile silently won the assignment -- and because the
# sanitizer targets work by setting OPT in a sub-make, OPT never reached the
# compile line.  `CFLAGS=-O2 make ubsan` then built with no -fsanitize and no
# warning set, ran the tests, and reported green: a sanitizer gate that checked
# nothing.  A makefile `=` beats the environment but still loses to a
# command-line `make CFLAGS=...`, which is the override that was actually
# intended.
CFLAGS    = $(OPT) $(WARN)
override REQ_CFLAGS = $(STD)
ALL_CFLAGS = $(REQ_CFLAGS) $(CFLAGS) $(DEFS) $(CPPFLAGS)
LDLIBS   = -lm -lz -lpthread

# .nii.gz support is compiled into nifti_io only when HAVE_ZLIB is defined; without
# it every compressed read fails at run time with a "recompile for compressed data"
# message, which is a confusing way to discover a build flag is missing.
override DEFS += -DHAVE_ZLIB

# nifti_io.c is vendored verbatim and is not ours to tidy; build it without the
# strict warning set so a clean build stays meaningful for the code we own.
NIFTI_CFLAGS = $(REQ_CFLAGS) $(OPT) -w $(DEFS) $(CPPFLAGS)

# `force` below is the first real target in the file, which would otherwise
# become the default goal; say what the default is instead of relying on order.
.DEFAULT_GOAL := all

OBJS  = dn.o dn_coef.o dn_lambda_star.o dn_eig.o dn_patch.o dn_nii.o dn_run.o
MAIN  = svht_denoise.o
NIFTI = nifti_io.o

# Flag stamp.  Object files carry no record of the flags that built them, so a
# plain `make` after `make ubsan` happily relinks sanitiser-instrumented objects
# -- which fails at link time if you are lucky, and silently ships instrumented
# code if you are not.  Rebuild everything whenever the flags change.
# (niimath solves the same problem the same way; see its .al_obj_flags.)
FLAGSTAMP := .build_flags
CURFLAGS  := $(CC) $(ALL_CFLAGS) $(LDLIBS)

# The stamp is a real file that every object depends on, and it is REWRITTEN
# (so its mtime moves) only when the flags actually change.  Deleting objects
# from a recipe does not work: make plans the whole graph from the file state it
# saw at startup, so a target it already judged up to date stays that way.
.PHONY: force
force: ;
$(FLAGSTAMP): force
	@if [ "$$(cat $@ 2>/dev/null)" != "$(CURFLAGS)" ]; then \
		printf '%s' "$(CURFLAGS)" > $@; \
	fi

BIN     = svht_denoise

.PHONY: all debug ubsan verbose clean

all: $(BIN)

$(BIN): $(OBJS) $(MAIN) $(NIFTI)
	$(CC) $(ALL_CFLAGS) -o $@ $^ $(LDLIBS)

$(NIFTI): nifti_io.c nifti_io.h $(FLAGSTAMP)
	$(CC) $(NIFTI_CFLAGS) -c -o $@ nifti_io.c

# -MMD -MP makes the compiler emit each object's header dependencies, so editing
# dn.h or dn_eig.h rebuilds what includes it.  Without this the pattern rule
# depended only on the .c and the flag stamp, and a header change could link
# objects built against different declarations or constants.
%.o: %.c $(FLAGSTAMP)
	$(CC) $(ALL_CFLAGS) -MMD -MP -c -o $@ $<

-include $(OBJS:.o=.d) $(MAIN:.o=.d)

debug:
	$(MAKE) OPT="-O0 -g" all

ubsan:
	$(MAKE) OPT="-O1 -g -fsanitize=undefined -fno-sanitize-recover=all" all

verbose:
	$(MAKE) WARN="-Weverything -Wno-padded -Wno-declaration-after-statement -Wno-unsafe-buffer-usage -Wno-covered-switch-default -Wno-cast-align -Wno-double-promotion -Wno-float-equal -Wno-switch-default -Wno-pre-c11-compat" all

clean:
	rm -f *.o *.d $(BIN) $(FLAGSTAMP)

# Debug builds leave a .dSYM bundle beside the binary; sweep it too.
distclean: clean
	rm -rf *.dSYM

.PHONY: distclean
