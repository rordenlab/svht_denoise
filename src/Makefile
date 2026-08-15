# svht_denoise -- build.
#
# Targets
#   make             optimised build (default)
#   make debug       -O0 -g
#   make test        regression tests (uninstrumented)
#   make ubsan       UndefinedBehaviorSanitizer   (also: asan, tsan)
#   make ubsan-test  build AND run the tests instrumented (also: asan-test, tsan-test)
#                    Use these, NOT `make ubsan && make test` -- see the note below.
#   make verbose     -Weverything, for exploration; NOT clean, and not a gate
#   make clean
#   make distclean   also sweep .dSYM bundles
#   make macos-release   signed, notarized, stapled .pkg (see the section below)
#
# The only dependency beyond libm and pthreads is zlib, for .nii.gz via nifti_io;
# on macOS that is the bundled static zlib-ng by default -- see the zlib-ng
# section below for how to turn it off or point it elsewhere.
#
# There is deliberately NO -ffast-math anywhere.  This code validates finite
# values and promises byte-identical output across thread counts; reassociation
# and non-finite shortcuts are incompatible with both.

CC      ?= cc
OPT     ?= -O3
override STD = -std=c99
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

# ------------------------------------------------------------------ zlib-ng
#
# Writing .nii.gz is deflate-bound, and deflate is the largest serial block in a
# gz->gz run: ~0.9 s of a 2.4 s run on 104x104x72x21, single-threaded, while
# every other core idles.  zlib-ng compresses the same data ~2.9x faster than
# the system zlib 1.2.12 (0.23 s against 0.66 s for 77 MB) and produces a
# slightly SMALLER file, so there is no size/speed trade to weigh.  READING is
# unaffected -- inflate already runs at over 1 GB/s and was never the
# bottleneck.  Measured end to end: 2.13 s -> 1.68 s.
#
# On macOS it is therefore the default, linked STATICALLY, which also leaves the
# shipped binary with no libz dependency at all.  Sourced from the
# third_party/zlib-ng submodule pinned to an exact commit, so a given revision
# of this repo always produces the same binary.
#
#   make ZLIBNG=0              force the system -lz
#   make ZLIBNG=1              force zlib-ng (also enables it off macOS)
#   make ZLIBNG_ROOT=<dir>     use a prebuilt ZLIB_COMPAT tree instead
#
# A NON-recursive clone leaves the submodule directory empty.  That falls back
# to the system -lz rather than failing, and `make` says which one it used --
# the difference is speed, not behaviour, and refusing to build a working tool
# over it would be worse.  Populate it with, from this directory:
#   git submodule update --init ../third_party/zlib-ng
#
# ZNG_SRC and ZNG_BUILD are RELATIVE, and that is load-bearing rather than tidy.
# Both appear as make targets and prerequisites, and make separates those on
# whitespace, so an absolute path is unusable the moment the checkout lives at
# "/Users/me/My Project" -- `rm -rf $(ZNG_BUILD)` in `clean` then expanded to
# two arguments and DELETED "/Users/me/My".  Demonstrated, with a canary file.
# Quoting the shell uses is not a fix, because a target containing a space
# cannot be written in a makefile at all.  Relative paths sidestep it entirely:
# make always runs in src/, and neither ".zlibng" nor "../third_party/zlib-ng"
# can contain a space no matter where the repository sits.
ZNG_SRC   := ../third_party/zlib-ng
ZNG_BUILD := .zlibng
ZLIBNG    ?= $(if $(filter Darwin,$(UNAME_S)),1,0)
# Reject any other value rather than treating it as "off".  The gate used to be
# `ZLIBNG` equal to the literal `auto`, so `ZLIBNG=1` -- the obvious way to ask
# for it -- silently selected the SYSTEM zlib, and silently is the whole problem:
# the result still builds, still passes, and is only ~3x slower at writing.
ifeq ($(filter 0 1,$(ZLIBNG)),)
  $(error ZLIBNG must be 0 or 1 (got '$(ZLIBNG)'))
endif

ifdef ZLIBNG_ROOT
  # Prebuilt tree: must be a ZLIB_COMPAT build (drop-in libz.a + zlib.h).  The
  # packaged zlib-ng builds (Homebrew, distros) are NOT compat -- they install a
  # zng_-prefixed zlib-ng.h -- so they cannot satisfy nifti_io's <zlib.h>.
  ZNG_LIB := $(firstword $(wildcard $(ZLIBNG_ROOT)/libz.a $(ZLIBNG_ROOT)/lib/libz.a))
  ZNG_INC := $(dir $(firstword $(wildcard $(ZLIBNG_ROOT)/zlib.h $(ZLIBNG_ROOT)/include/zlib.h)))
  ifeq ($(ZNG_LIB),)
    $(error ZLIBNG_ROOT is set but no libz.a under $(ZLIBNG_ROOT) or $(ZLIBNG_ROOT)/lib -- run ./configure --zlib-compat && make there)
  endif
  ifeq ($(ZNG_INC),)
    $(error ZLIBNG_ROOT is set but no zlib.h under $(ZLIBNG_ROOT) or $(ZLIBNG_ROOT)/include -- use the same ZLIB_COMPAT build for headers and libz.a)
  endif
else ifeq ($(ZLIBNG),1)
  # Our own paths are relative and so immune to the checkout location, but
  # zlib-ng's build is NOT: its configure resolves an absolute SRCDIR and pastes
  # it unquoted into its generated recipes.  A space merely breaks it ("No rule
  # to make target arch/generic/adler32_c.o"), but a shell metacharacter is
  # executed -- a directory named with backticks around a command ran that
  # command, demonstrated with a canary.  A CI checkout path derived from a
  # branch or tag name is the realistic way that becomes someone else's input.
  # So the test is an ALLOW-LIST of characters known to survive being pasted
  # into a recipe, not a hunt for the ones that do not, and anything else is
  # treated exactly like a missing submodule: fall back, say why, keep building.
  # Checking the shell's $$PWD rather than $(CURDIR) so the value is never make
  # syntax either.  Upstream's quoting is not ours to fix.
  ZNG_PATH_OK := $(shell case "$$PWD" in (*[!A-Za-z0-9._/-]*) echo no;; (*) echo yes;; esac)
  ifneq ($(ZNG_PATH_OK),yes)
    ZNG_SKIP := unsafe character in the build path
  else ifeq ($(wildcard $(ZNG_SRC)/configure),)
    ZNG_SKIP := submodule not checked out
  else
    ZNG_LIB   := $(ZNG_BUILD)/libz.a
    ZNG_INC   := $(ZNG_BUILD)
  endif
endif

ifeq ($(ZNG_LIB),)
  ZLIB_LIB    = -lz
  ZLIB_WHICH  = system zlib$(if $(ZNG_SKIP), -- zlib-ng skipped: $(ZNG_SKIP))
else
  # zlib-ng's configure GENERATES zlib.h and zconf.h into the build directory,
  # so the include path is the build tree, not the source tree.
  ZLIB_CFLAGS = -I$(ZNG_INC)
  ZLIB_LIB    = $(ZNG_LIB)
  ZLIB_WHICH  = zlib-ng (static)
endif

# ------------------------------------------------------------ Apple Accelerate
#
# The split path's Householder reduction and back-transform come from LAPACK's
# dsytrd and dormtr, out of Accelerate.  DEFAULT ON for macOS, `ACCELERATE=0`
# to opt out; off everywhere else, since no other platform has the framework.
#
# Measured on 100x100x58x138: -14.4% wall and -16.2% CPU (58.6 -> 50.2 s,
# 741 -> 621 s).  The full solver path (n < DN_SPLIT_MIN) is untouched, so the
# small-N benchmark series is unaffected -- this is a large-N win only.
#
# It is a framework, not a library to install or bundle: present on every macOS
# that this builds for, resolved at load like libSystem, and it adds nothing to
# the executable.  That is why it is admitted where a BLAS is not.
#
# It changes the numbers, and the divergence needs a big sample to see: measured
# at ~3e-10 relative L2 with the rank map identical on a 100x100x58 x 138 uint16
# series (580k voxels), while the repository's own fixtures and a synthetic
# 16^3 x 138 all come out byte-identical between the two builds.  So "not
# bit-identical" is true but not reproducible on anything shipped here -- do not
# conclude from a green fixture that the two agree.  The Gram differs at every n
# rather than only above LAPACK's blocking crossover, though on integer-valued
# input (any uint16/int16 scan) the products and their sums are exact in double
# and dsyrk agrees bit for bit anyway.  Byte-identity ACROSS THREAD COUNTS holds
# within either build, which is the promise actually made.
ACCELERATE ?= $(if $(filter Darwin,$(UNAME_S)),1,0)
# Same hard error ZLIBNG gets, for the same reason: `ACCELERATE=yes` would
# otherwise select the portable reduction silently and ship the slower kernel.
ifeq ($(filter 0 1,$(ACCELERATE)),)
  $(error ACCELERATE must be 0 or 1 (got '$(ACCELERATE)'))
endif
ifeq ($(ACCELERATE),1)
  override REQ_CFLAGS += -DDN_USE_ACCELERATE -DACCELERATE_NEW_LAPACK
  ACCEL_LIB = -framework Accelerate
  ACCEL_WHICH = Accelerate (dsytrd/dormtr)
else
  ACCEL_WHICH = the built-in reduction
endif

ALL_CFLAGS = $(REQ_CFLAGS) $(CFLAGS) $(DEFS) $(ZLIB_CFLAGS) $(CPPFLAGS)
LDLIBS   = -lm $(ZLIB_LIB) -lpthread $(ACCEL_LIB)

# .nii.gz support is compiled into nifti_io only when HAVE_ZLIB is defined; without
# it every compressed read fails at run time with a "recompile for compressed data"
# message, which is a confusing way to discover a build flag is missing.
override DEFS += -DHAVE_ZLIB

# nifti_io.c is vendored verbatim and is not ours to tidy; build it without the
# strict warning set so a clean build stays meaningful for the code we own.
# ZLIB_CFLAGS must be here too: nifti_io.c is the translation unit that actually
# includes <zlib.h>, so leaving it out builds the compressed path against the
# SYSTEM header while linking zlib-ng's library.
NIFTI_CFLAGS = $(REQ_CFLAGS) $(OPT) -w $(DEFS) $(ZLIB_CFLAGS) $(CPPFLAGS)

# `force` below is the first real target in the file, which would otherwise
# become the default goal; say what the default is instead of relying on order.
.DEFAULT_GOAL := all

OBJS  = dn.o dn_coef.o dn_lambda_star.o dn_eig.o dn_patch.o dn_nii.o dn_phase.o dn_run.o
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
# (so its mtime moves) only when the flags actually change.
#
# The stale objects are removed at PARSE TIME, immediately below, and that is
# load-bearing rather than tidy.  A recipe cannot do it: make plans the whole
# graph from the file state it saw at startup, so a target it has already judged
# up to date stays that way even after the recipe deletes it.  Depending on the
# stamp's MTIME instead is defeated by `-j`, because GNU Make 3.81 -- which is
# what /usr/bin/make is on macOS -- plans each object against the stamp as it
# stood BEFORE the stamp's own recipe rewrote it.  Both failures are silent and
# both were measured:
#
#   make -j8 && make -j8 ACCELERATE=0     link never runs; the binary still has
#                                         Accelerate in it while the build says
#                                         "the built-in reduction"
#   make -j8 && make -j8 ubsan-test       78 passed, with the eigensolver and
#                                         patch kernel UNINSTRUMENTED
#
# The second is the sanitizer gate this file already believes it closed, still
# open under -j.  A MISSING object cannot be mis-dated, which is why deleting
# beats any amount of dependency ordering here.
ifneq ($(strip $(CURFLAGS)),$(strip $(shell cat $(FLAGSTAMP) 2>/dev/null)))
  $(shell rm -f *.o svht_denoise dn_testgen)
endif
.PHONY: force
force: ;
$(FLAGSTAMP): force
	@if [ "$$(cat $@ 2>/dev/null)" != "$(CURFLAGS)" ]; then \
		printf '%s' "$(CURFLAGS)" > $@; \
	fi

BIN     = svht_denoise

.PHONY: all debug ubsan asan tsan ubsan-test asan-test tsan-test verbose clean test \
        macos-pkg-adhoc macos-pkg macos-notarize macos-notary-profile \
        macos-release check-notary-profile distclean

# The echo lives HERE, not in the $(BIN) recipe, so it prints on an
# already-built tree as well as a fresh one.  That is the state you are in right
# before a benchmark, which is precisely when you need to know which library
# went in -- and in the link recipe it printed nothing at all.
all: $(BIN)
	@echo "  $(BIN): compression via $(ZLIB_WHICH), reduction via $(ACCEL_WHICH)"

# Regression tests for the cases that were once silent-failure bugs.  The fixture
# generator is deliberately standalone -- it writes NIfTI by hand rather than
# through nifti_io, so a fault there cannot hide by appearing on both sides of a
# comparison.  To run it instrumented use `make ubsan-test` and friends, never
# `make ubsan && make test` -- the reason is on those targets below.
TESTGEN = dn_testgen

$(TESTGEN): dn_testgen.c $(FLAGSTAMP)
	$(CC) $(ALL_CFLAGS) -o $@ dn_testgen.c -lm

test: $(BIN) $(TESTGEN)
	@BIN=./$(BIN) GEN=./$(TESTGEN) sh test.sh

# Say which zlib went in.  The two differ in .nii.gz WRITE speed by about 3x, so
# silently falling back to the system one -- which is what a non-recursive clone
# does -- would otherwise look like the tool had simply got slower.
# `$(filter %.o,$^)`, not `$^`: the zlib-ng archive is a prerequisite so that
# rebuilding it relinks, but it reaches the link line through $(LDLIBS).  Passing
# $^ wholesale listed it twice and drew "ld: warning: ignoring duplicate
# libraries", and a warning nobody can fix is how a build stops being read.
$(BIN): $(OBJS) $(MAIN) $(NIFTI) $(ZNG_LIB)
	$(CC) $(ALL_CFLAGS) -o $@ $(filter %.o,$^) $(LDLIBS)

# $(ZNG_LIB) is a prerequisite of BOTH rules above, not just a link input.  For
# this one because zlib-ng's configure is what GENERATES zlib.h, so nifti_io.c
# -- the only translation unit that includes it -- cannot be compiled before
# configure has run, and must be recompiled when the header changes underneath
# it.  Empty when the system zlib is in use, so it costs nothing there.
$(NIFTI): nifti_io.c nifti_io.h $(FLAGSTAMP) $(ZNG_LIB)
	$(CC) $(NIFTI_CFLAGS) -c -o $@ nifti_io.c

# Ties the archive to the submodule's CHECKED-OUT REVISION.  Without it the
# archive has no prerequisites at all, so `git submodule update` to a different
# commit silently relinks the PREVIOUS library: the pinned-revision
# reproducibility this whole arrangement exists for would hold for releases,
# which start with `make clean`, and quietly not hold for ordinary builds and
# benchmarks.  Same rewrite-only-when-it-changes trick as $(FLAGSTAMP), for the
# same reason.  `no-git` when git is unavailable, as in an exported tarball.
# Kept OUTSIDE $(ZNG_BUILD) so the recipe below can delete that directory
# wholesale without destroying its own prerequisite.
ZNG_REV := .zlibng-rev
$(ZNG_REV): force
	@rev=`git -C $(ZNG_SRC) rev-parse HEAD 2>/dev/null || echo no-git`; \
	if [ "`cat $@ 2>/dev/null`" != "$$rev" ]; then printf '%s' "$$rev" > $@; fi

# Wipes the build directory first, so a revision change gives a FROM-SCRATCH
# build.  Re-running configure over the old tree would leave objects belonging
# to the previous revision and rebuild only those whose sources moved, which is
# not the same thing as "this commit produces this binary" -- and that claim is
# the entire reason the submodule is pinned.  This is why $(ZNG_REV) lives
# outside the directory: the recipe would otherwise delete its own prerequisite.
# `cd`s into the build directory, so the source path needs one more `../`.
$(ZNG_BUILD)/libz.a: $(ZNG_REV)
	@echo "  building zlib-ng (submodule, static)"
	@rm -rf $(ZNG_BUILD) && mkdir -p $(ZNG_BUILD)
	@cd $(ZNG_BUILD) && ../$(ZNG_SRC)/configure --zlib-compat --static >config.log 2>&1 || { \
		echo "zlib-ng configure failed; last 20 lines of $(ZNG_BUILD)/config.log:" >&2; \
		tail -20 config.log >&2; exit 1; }
	@$(MAKE) --no-print-directory -C $(ZNG_BUILD) libz.a >$(ZNG_BUILD)/build.log 2>&1 || { \
		echo "zlib-ng build failed; last 20 lines of $(ZNG_BUILD)/build.log:" >&2; \
		tail -20 $(ZNG_BUILD)/build.log >&2; exit 1; }
	@touch $@

# -MMD -MP makes the compiler emit each object's header dependencies, so editing
# dn.h or dn_eig.h rebuilds what includes it.  Without this the pattern rule
# depended only on the .c and the flag stamp, and a header change could link
# objects built against different declarations or constants.
%.o: %.c $(FLAGSTAMP)
	$(CC) $(ALL_CFLAGS) -MMD -MP -c -o $@ $<

-include $(OBJS:.o=.d) $(MAIN:.o=.d)

debug:
	$(MAKE) OPT="-O0 -g" all

# Sanitizer flags live in one variable each, so the build target and the test
# target below cannot drift apart.
UBSAN_OPT = -O1 -g -fsanitize=undefined -fno-sanitize-recover=all
ASAN_OPT  = -O1 -g -fsanitize=address
TSAN_OPT  = -O1 -g -fsanitize=thread

ubsan:
	$(MAKE) OPT="$(UBSAN_OPT)" all
asan:
	$(MAKE) OPT="$(ASAN_OPT)" all
tsan:
	$(MAKE) OPT="$(TSAN_OPT)" all

# Run the suite AGAINST INSTRUMENTED CODE.  These exist because
# `make ubsan && make test` does not do that: the second command runs at the
# default OPT, so the flag stamp no longer matches, everything is rebuilt
# uninstrumented, and the suite reports green having exercised nothing.  It is
# the same trap the FLAGSTAMP comment above describes, one level up -- and worse,
# it did not reproduce reliably, so it could be missed by a spot check.  Building
# and testing inside ONE recursive make keeps a single OPT for both halves.
ubsan-test:
	$(MAKE) OPT="$(UBSAN_OPT)" test
asan-test:
	$(MAKE) OPT="$(ASAN_OPT)" test
tsan-test:
	$(MAKE) OPT="$(TSAN_OPT)" test

verbose:
	$(MAKE) WARN="-Weverything -Wno-padded -Wno-declaration-after-statement -Wno-unsafe-buffer-usage -Wno-covered-switch-default -Wno-cast-align -Wno-double-promotion -Wno-float-equal -Wno-switch-default -Wno-pre-c11-compat" all

# ---------------------------------------------------------------- macOS release
#
#   make macos-notary-profile APPLE_ID='you@example.com' TEAM_ID='ABCDE12345'
#   make macos-release VERSION=<version> \
#        MACOS_SIGN_IDENTITY='Developer ID Application: Name (TEAMID)' \
#        MACOS_INSTALLER_IDENTITY='Developer ID Installer: Name (TEAMID)'
#
# The profile is stored once per machine; the release step builds, tests, signs,
# packages, notarises, staples and verifies. `make macos-pkg-adhoc` is the
# unsigned local-testing path.
# Written out in full, and DELIBERATELY not derived from svht_denoise.c's
# DN_VERSION even though the two must agree.  Deriving it was tried: the
# `override ... := $(value ...)` freeze below exists to guarantee this variable
# holds literal characters and never make syntax, so a default of
# `$(DN_VERSION_H)` is frozen as those nine literal characters and the release
# ships an empty version.  Keeping the number here is the safe half of that
# trade; bump both when it changes.
VERSION ?= 0.1.20260808
MACOS_SIGN_IDENTITY ?=
MACOS_INSTALLER_IDENTITY ?=
NOTARY_PROFILE ?= svht-denoise-notary

# EXPORTED, and referenced below as shell variables ($$VERSION), never as make
# text ($(VERSION)).  That is a security property, not a style choice: make
# pastes a variable's value into the recipe as literal shell source, so a value
# containing a quote closes it and the rest runs as a command.  Demonstrated:
#   make macos-notarize VERSION='x";touch /tmp/PWNED;"'
# created the file.  package_macos.sh validates VERSION, but far too late -- the
# shell has already run the injected command before the script starts.  It
# matters because the documented flow is `make macos-release VERSION=<tag>`, and
# a CI job feeding a tag or branch name into that would be arbitrary code
# execution on the machine holding the Developer ID signing key.  Passed through
# the environment, a value can never become syntax.
# `export` alone is NOT enough.  A command-line assignment is still expanded as
# MAKE text when the child environment is built, so `$(shell ...)` inside the
# value runs -- on every target, including a plain `make all`.  Demonstrated:
#   make all VERSION='$$(shell touch /tmp/PWNED)'
# `$(value X)` takes the raw unexpanded text and `:=` freezes it, so what ships
# is the literal characters, which package_macos.sh's own validation rejects.
# Values arriving through the environment are not re-expanded; only
# command-line assignments are.
override VERSION                 := $(value VERSION)
override MACOS_SIGN_IDENTITY     := $(value MACOS_SIGN_IDENTITY)
override MACOS_INSTALLER_IDENTITY := $(value MACOS_INSTALLER_IDENTITY)
override NOTARY_PROFILE          := $(value NOTARY_PROFILE)
override APPLE_ID                := $(value APPLE_ID)
override TEAM_ID                 := $(value TEAM_ID)
export VERSION MACOS_SIGN_IDENTITY MACOS_INSTALLER_IDENTITY NOTARY_PROFILE APPLE_ID TEAM_ID

# Unsigned, for local testing. Gatekeeper will refuse it on another machine.
macos-pkg-adhoc:
	env MACOS_SIGN_IDENTITY=- MACOS_INSTALLER_IDENTITY= ../scripts/package_macos.sh "$$VERSION"

# Two DIFFERENT certificates, both from the same Apple Developer account:
# "Developer ID Application" signs the executable, "Developer ID Installer"
# signs the .pkg. Signing a package with the Application certificate produces a
# package the installer rejects, so both are required here rather than one.
macos-pkg:
	@test -n "$$MACOS_SIGN_IDENTITY" || { \
		echo "Set MACOS_SIGN_IDENTITY to a Developer ID Application identity or SHA-1." >&2; \
		echo "List them with: security find-identity -v" >&2; \
		exit 2; \
	}
	@test -n "$$MACOS_INSTALLER_IDENTITY" || { \
		echo "Set MACOS_INSTALLER_IDENTITY to a Developer ID Installer identity or SHA-1." >&2; \
		echo "List them with: security find-identity -v" >&2; \
		exit 2; \
	}
	../scripts/package_macos.sh "$$VERSION"

macos-notarize:
	../scripts/notarize_macos.sh "../dist/svht_denoise-$$VERSION-macos-arm64.pkg" "$$NOTARY_PROFILE"

macos-notary-profile:
	@case "$$NOTARY_PROFILE" in -*|"") \
		echo "NOTARY_PROFILE must not be empty or start with '-': notarytool would read it as an option." >&2; \
		exit 2;; esac
	@test -n "$$APPLE_ID" || { echo "Set APPLE_ID to your Developer Apple ID." >&2; exit 2; }
	@test -n "$$TEAM_ID" || { echo "Set TEAM_ID to your 10-character Developer Team ID." >&2; exit 2; }
	../scripts/store_notary_credentials.sh "$$NOTARY_PROFILE" "$$APPLE_ID" "$$TEAM_ID"

# Sub-makes rather than prerequisites, so the order holds under `make -j`:
# prerequisites may run concurrently, which would defeat the whole point of
# checking the credentials before spending minutes on a build and a signature.
macos-release:
	$(MAKE) check-notary-profile
	$(MAKE) macos-pkg
	../scripts/notarize_macos.sh "../dist/svht_denoise-$$VERSION-macos-arm64.pkg" "$$NOTARY_PROFILE"

# Ask notarytool, NOT `security find-generic-password`. `notarytool
# store-credentials` writes to the data-protection keychain, while `security`
# only queries the legacy file-based one, so that lookup fails for every
# profile including working ones. There is no offline query for that keychain.
# Asking notarytool costs a network round trip, which notarization needs
# anyway, and is strictly stronger: it proves the credentials AUTHENTICATE
# rather than merely that an item exists.
check-notary-profile:
	@xcrun notarytool history --keychain-profile "$$NOTARY_PROFILE" \
		>/dev/null 2>&1 || { \
		echo "No WORKING notarization credentials for profile '$$NOTARY_PROFILE'." >&2; \
		echo "(Stored but rejected credentials fail here too -- so do network outages.)" >&2; \
		echo "Store them once with:" >&2; \
		echo "  make macos-notary-profile APPLE_ID='you@example.com' TEAM_ID='ABCDE12345'" >&2; \
		echo "That prompts for an APP-SPECIFIC password (appleid.apple.com ->" >&2; \
		echo "Sign-In and Security -> App-Specific Passwords), not your Apple ID password." >&2; \
		exit 2; \
	}
	@echo "notarization credentials found for profile '$$NOTARY_PROFILE'"

# The zlib-ng build goes too, and that is deliberate rather than thorough.  It
# carries a deployment target inherited from the environment at the time it was
# configured, and package_macos.sh's first step is `make clean`: keeping a dev
# build here would let a library configured for the host's macOS be linked into
# a release that advertises 11.0.  The symptom is a wall of link warnings and an
# executable whose own minos still reads 11.0, so the package check passes and
# nothing fails.  Rebuilding it costs a couple of seconds.
clean:
	rm -f *.o *.d $(BIN) $(TESTGEN) $(FLAGSTAMP)
	rm -rf $(ZNG_BUILD) $(ZNG_REV)

# Debug builds leave a .dSYM bundle beside the binary; sweep it too.
distclean: clean
	rm -rf *.dSYM
