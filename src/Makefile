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

.PHONY: all debug ubsan asan tsan ubsan-test asan-test tsan-test verbose clean test \
        macos-pkg-adhoc macos-pkg macos-notarize macos-notary-profile \
        macos-release check-notary-profile distclean

all: $(BIN)

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
#   make macos-release VERSION=1.1.0 \
#        MACOS_SIGN_IDENTITY='Developer ID Application: Name (TEAMID)' \
#        MACOS_INSTALLER_IDENTITY='Developer ID Installer: Name (TEAMID)'
#
# The profile is stored once per machine; the release step builds, tests, signs,
# packages, notarises, staples and verifies. `make macos-pkg-adhoc` is the
# unsigned local-testing path.
VERSION ?= 1.1.0
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
	../scripts/notarize_macos.sh "../dist/svht_denoise-$$VERSION-macos-universal.pkg" "$$NOTARY_PROFILE"

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
	../scripts/notarize_macos.sh "../dist/svht_denoise-$$VERSION-macos-universal.pkg" "$$NOTARY_PROFILE"

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

clean:
	rm -f *.o *.d $(BIN) $(TESTGEN) $(FLAGSTAMP)

# Debug builds leave a .dSYM bundle beside the binary; sweep it too.
distclean: clean
	rm -rf *.dSYM
