#!/bin/sh
# Regression tests.  Run with `make test`.
#
# These are the cases that were once BUGS, not a general test suite: each one
# either silently destroyed a file, silently corrupted the data, or silently did
# nothing while reporting success.  Guard-level checks that already fail loudly
# are deliberately not duplicated here.
#
# POSIX sh, no framework, no network, no Python.  Fixtures come from dn_testgen,
# which writes NIfTI by hand so that a fault in nifti_io cannot hide by being
# present on both sides of a comparison.

set -u

BIN=${BIN:-./svht_denoise}
GEN=${GEN:-./dn_testgen}
WORK=$(mktemp -d "${TMPDIR:-/tmp}/svht_test.XXXXXX") || exit 1
trap 'rm -rf "$WORK"' EXIT INT TERM

pass=0; fail=0
ok()   { pass=$((pass+1)); printf '  ok   %s\n' "$1"; }
bad()  { fail=$((fail+1)); printf '  FAIL %s\n' "$1"; }
# check <name> <expected-exit> <cmd...>
check() {
	name=$1; want=$2; shift 2
	"$@" >"$WORK/stdout" 2>"$WORK/stderr"; got=$?
	if [ "$got" -eq "$want" ]; then ok "$name"; else
		bad "$name (exit $got, wanted $want)"; sed 's/^/       /' "$WORK/stderr"; fi
}

# check_collision <name> <cmd...>: exit 1 AND the collision refusal.
# An exit-code-only assertion is not enough for these: a dangling or cyclic
# output ALSO fails at the write with exit 1, so the test would pass with the
# collision check entirely disabled -- measured, it did.
check_collision() {
	name=$1; shift
	"$@" >"$WORK/stdout" 2>"$WORK/stderr"; got=$?
	if [ "$got" -eq 1 ] && grep -q "both resolve to" "$WORK/stderr"
	then ok "$name"
	else bad "$name (exit $got, and no collision refusal)"; sed 's/^/       /' "$WORK/stderr"; fi
}

for k in mag mask phase-rad phase-int phase-turns phase-deg phase-affine phase-const phase-siemens; do
	"$GEN" mk "$k" "$WORK/$k.nii" || { echo "fixture $k failed"; exit 1; }
done
M=$WORK/mag.nii

echo "path safety"
# An output must never land on an input or on another output.  Each of these
# once ran to completion and destroyed a file.
check_collision "output over input" "$BIN" "$M" "$WORK/o.nii" -quiet -noise "$M"
check_collision "two outputs, same name" "$BIN" "$M" "$WORK/o.nii" -quiet -noise "$WORK/n.nii" -rank "$WORK/n.nii"
check_collision "input reused as -phase" "$BIN" "$M" "$WORK/o.nii" -quiet -phase "$M"
# Case-only collision: real on macOS, not on a case-sensitive filesystem, so the
# expectation follows the filesystem rather than being hardcoded.
: > "$WORK/CaseProbe"
if [ -e "$WORK/caseprobe" ]; then want=1; else want=0; fi
rm -f "$WORK/CaseProbe"
if [ "$want" -eq 1 ]
then check_collision "case-only output clash" "$BIN" "$M" "$WORK/o.nii" -quiet -noise "$WORK/N.nii" -rank "$WORK/n.nii"
else check "case-only names are distinct here" 0 "$BIN" "$M" "$WORK/o.nii" -quiet -noise "$WORK/N.nii" -rank "$WORK/n.nii"; fi
# Unicode-equivalent names. Apple filesystems fold U+212A KELVIN onto "k" and
# treat NFC and NFD as one name -- the latter even when formatted case-SENSITIVE.
# Probed rather than assumed, exactly like the case test: elsewhere these really
# are different files and must be allowed.
sv_nfc=$(printf '\303\251')          # e-acute, composed
sv_nfd=$(printf 'e\314\201')         # e + combining acute
: > "$WORK/$sv_nfc.probe"
if [ -e "$WORK/$sv_nfd.probe" ]; then sv_uni=1; else sv_uni=0; fi
rm -f "$WORK/$sv_nfc.probe"
if [ "$sv_uni" -eq 1 ]; then
	check_collision "NFC/NFD outputs collide" "$BIN" "$M" "$WORK/o.nii" -quiet \
		-noise "$WORK/$sv_nfc.nii" -rank "$WORK/$sv_nfd.nii"
	check_collision "KELVIN folds onto k" "$BIN" "$M" "$WORK/o.nii" -quiet \
		-noise "$WORK/k.nii" -rank "$WORK/$(printf '\342\204\252').nii"
else
	check "NFC/NFD are distinct here" 0 "$BIN" "$M" "$WORK/o.nii" -quiet \
		-noise "$WORK/$sv_nfc.nii" -rank "$WORK/$sv_nfd.nii"
fi
rm -f "$WORK/$sv_nfc.nii" "$WORK/$sv_nfd.nii" "$WORK/k.nii" "$WORK/o.nii"
# ...and a non-ASCII DIRECTORY must not make every pair of outputs collide.
mkdir -p "$WORK/$sv_nfc dir" && cp "$M" "$WORK/$sv_nfc dir/in.nii"
check "non-ASCII directory" 0 "$BIN" "$WORK/$sv_nfc dir/in.nii" "$WORK/$sv_nfc dir/out.nii" \
	-quiet -noise "$WORK/$sv_nfc dir/noise.nii"
if [ -s "$WORK/$sv_nfc dir/out.nii" ] && [ -s "$WORK/$sv_nfc dir/noise.nii" ]
then ok "non-ASCII directory wrote both outputs"; else bad "non-ASCII directory lost an output"; fi

# The .hdr/.img form, which is the ONLY shape that reaches the input-resolution
# bug: the reader expands the prefix "p" to the existing p.hdr/p.img, while the
# rule for OUTPUT names expands it to p.nii.  Comparing the two spellings by the
# output rule once let a noise map overwrite the phase input and exit 0.
# A single-.nii fixture cannot reach this, so the pair is written explicitly.
"$GEN" mk "phase-rad@" "$WORK/p" || { echo "pair fixture failed"; exit 1; }
check_collision "hdr/img input aliasing" "$BIN" "$M" "$WORK/o.nii" -quiet -phase "$WORK/p" -noise "$WORK/p.hdr"
# ...and the pair must still be usable as an ordinary input.
check "hdr/img pair as -phase"   0 "$BIN" "$M" "$WORK/o.nii" -quiet -phase "$WORK/p"
# A DANGLING symlink: stat() fails on both sides, so the fallback compared the
# link's own name with its target's and let one output overwrite the other.
rm -f "$WORK/alias.nii" "$WORK/tgt.nii"
ln -s tgt.nii "$WORK/alias.nii"
check_collision "dangling symlink output" "$BIN" "$M" "$WORK/alias.nii" -quiet -noise "$WORK/tgt.nii"
rm -f "$WORK/alias.nii" "$WORK/tgt.nii"
# A chain longer than one link, still dangling, still an alias for the other
# output.  40 because Linux follows that many, so a resolver stopping earlier
# would leave it unresolved here and let the kernel resolve it at write time.
sv_i=0
while [ "$sv_i" -lt 40 ]; do
	sv_next=$((sv_i + 1))
	[ "$sv_next" -eq 40 ] && ln -s tgt.nii "$WORK/ln$sv_i.nii" || ln -s "ln$sv_next.nii" "$WORK/ln$sv_i.nii"
	sv_i=$sv_next
done
check_collision "40-link symlink chain" "$BIN" "$M" "$WORK/ln0.nii" -quiet -noise "$WORK/tgt.nii"
rm -f "$WORK"/ln* "$WORK/tgt.nii"
# A cycle cannot be resolved at all.  The names in hand are then not the names
# the kernel will use, so the only safe answer is to refuse.
ln -s cyc_b.nii "$WORK/cyc_a.nii"; ln -s cyc_a.nii "$WORK/cyc_b.nii"
check_collision "symlink cycle as output" "$BIN" "$M" "$WORK/cyc_a.nii" -quiet -noise "$WORK/cyc_b.nii"
rm -f "$WORK/cyc_a.nii" "$WORK/cyc_b.nii"
# The other direction: a long chain that collides with NOTHING must still be
# accepted.  Without this, a hop limit set too low looks fine -- exhaustion
# refuses conservatively, so every refusal test above still passes.  20 links,
# not 40, because macOS SYMLOOP_MAX is 32 and the kernel itself cannot write
# through a 40-link chain; the 40-link ACCEPT case is Linux-only.
sv_i=0
while [ "$sv_i" -lt 20 ]; do
	sv_next=$((sv_i + 1))
	[ "$sv_next" -eq 20 ] && ln -s okchain.nii "$WORK/ok$sv_i.nii" || ln -s "ok$sv_next.nii" "$WORK/ok$sv_i.nii"
	sv_i=$sv_next
done
check "20-link chain accepted"   0 "$BIN" "$M" "$WORK/ok0.nii" -quiet -noise "$WORK/sep.nii"
[ -s "$WORK/okchain.nii" ] && ok "20-link chain wrote through to its target" \
	|| bad "20-link chain wrote nothing"
rm -f "$WORK"/ok*.nii "$WORK/okchain.nii" "$WORK/sep.nii"

# -mask is an INPUT, so an output named for it must be refused too. This row of
# the collision table had no coverage, though the fixture was already generated.
# -real is the only output that writes a full 4D series, so it is the one that can
# overwrite a magnitude input outright -- and its row of the collision table had
# no coverage: deleting the row left the whole suite green while the run
# destroyed its own input at exit 0.
check_collision "input reused as -real" "$BIN" "$M" "$WORK/o.nii" -quiet -phase "$WORK/phase-rad.nii" -real "$M"
check_collision "output over -mask" "$BIN" "$M" "$WORK/o.nii" -quiet -mask "$WORK/mask.nii" -noise "$WORK/mask.nii"
check "-mask accepted normally"  0 "$BIN" "$M" "$WORK/masked.nii" -quiet -mask "$WORK/mask.nii"
# ...and it must actually mask: an ignored -mask would still exit 0.
check "unmasked run"             0 "$BIN" "$M" "$WORK/unmasked.nii" -quiet
"$GEN" cmp "$WORK/masked.nii" "$WORK/unmasked.nii" 0 >/dev/null 2>&1; sv_rc=$?
if [ "$sv_rc" -eq 1 ]; then ok "-mask actually restricts the output"
elif [ "$sv_rc" -eq 0 ]; then bad "-mask changed nothing; it may be ignored"
else bad "-mask comparison failed to run (cmp exit $sv_rc)"; fi

# same_file() CREATES not-yet-existing names to compare them by inode, and must
# remove only what it created.  Keying that off "the open succeeded" instead
# deletes an existing INPUT -- pinned directly here rather than relying on later
# tests happening to notice the file is gone.
cp "$M" "$WORK/keep.nii"
check "ordinary run" 0 "$BIN" "$M" "$WORK/keep_out.nii" -quiet -noise "$WORK/keep_n.nii"
if cmp -s "$M" "$WORK/keep.nii"
then ok "the input survives the collision check"; else bad "the input was modified or removed"; fi
rm -f "$WORK/keep.nii" "$WORK/keep_out.nii" "$WORK/keep_n.nii"

echo "argument validation"
check "\"-\" rejected as output"  1 "$BIN" "$M" - -quiet
check "\"-\" rejected as -noise"  1 "$BIN" "$M" "$WORK/o.nii" -quiet -noise -
check "-real needs -phase"       1 "$BIN" "$M" "$WORK/o.nii" -quiet -real "$WORK/r.nii"
check "-phaseunits needs -phase" 1 "$BIN" "$M" "$WORK/o.nii" -quiet -phaseunits auto
check "bad -phaseunits value"    1 "$BIN" "$M" "$WORK/o.nii" -quiet -phase "$WORK/phase-rad.nii" -phaseunits furlongs
check "missing -phase file"      1 "$BIN" "$M" "$WORK/o.nii" -quiet -phase "$WORK/absent.nii"

echo "output name resolution"
# nifti_set_filenames() rewrites nim->nifti_type from the names it builds, and all
# outputs share one nim.  A .hdr output therefore once flipped every LATER
# extensionless output from .nii to .hdr/.img -- names the check never modelled --
# which landed the denoised series on the -phase input at exit 0.
"$GEN" mk "phase-rad@" "$WORK/q" || exit 1
cp "$WORK/q.img" "$WORK/q.img.orig"
rm -f "$WORK/q.nii" "$WORK/rr.hdr"
# Exit code AND the expected outputs, not just "the input survived": without
# these the test also passes when the run fails before writing anything, which
# it did under an unrelated injected regression.
check "type-flip run succeeds"   0 "$BIN" "$M" "$WORK/q" -quiet -phase "$WORK/q.hdr" -real "$WORK/rr.hdr"
if [ -f "$WORK/q.nii" ] && [ -f "$WORK/rr.hdr" ]
then ok "type-flip run wrote both outputs"; else bad "type-flip run did not write its outputs"; fi
if cmp -s "$WORK/q.img" "$WORK/q.img.orig"
then ok "a .hdr output does not redirect a later one onto an input"
else bad "the -phase input was overwritten by a type-flipped output name"; fi
# An extensionless output must still resolve to .nii when nothing flipped the type.
rm -f "$WORK/e.nii" "$WORK/e.hdr"
check "extensionless output"     0 "$BIN" "$M" "$WORK/e" -quiet
[ -f "$WORK/e.nii" ] && ok "extensionless output became .nii" || bad "extensionless output did not become .nii"
# stdin: nifti accepts "-", and there is no file for an output to collide with.
if "$BIN" - "$WORK/si.nii" -quiet < "$M" >/dev/null 2>&1 && [ -s "$WORK/si.nii" ]
then ok "input on stdin"; else bad "input on stdin"; fi

echo "geometry"
# Equal dimensions are not the same picture: this pair differs only in voxel size.
check "phase off the grid"       1 "$BIN" "$M" "$WORK/o.nii" -quiet -phase "$WORK/phase-affine.nii"
# phase-affine is the same 9x9x3 shape, so it trips the VOLUME-count guard;
# mask-wrongdim is single-volume on a different grid, which is the only way to
# reach the dimension guard.  Both are real rejections and both are checked.
check "-mask wrong volumes"      1 "$BIN" "$M" "$WORK/o.nii" -quiet -mask "$WORK/phase-affine.nii"
"$GEN" mk mask-wrongdim "$WORK/mask-bad.nii" || { echo "fixture mask-wrongdim failed"; exit 1; }
check "-mask wrong dimensions"   1 "$BIN" "$M" "$WORK/o.nii" -quiet -mask "$WORK/mask-bad.nii"
# A NaN in the sform must FAIL CLOSED.  It did not: NaN comparisons are false, so
# the running maximum never updated and the offset came back as a clean 0 mm.
"$GEN" mk phase-nan "$WORK/phase-nan.nii" || exit 1
check "phase with a NaN sform"   1 "$BIN" "$M" "$WORK/o.nii" -quiet -phase "$WORK/phase-nan.nii"

echo "phase units"
# Every encoding is the SAME phase field, so all must give the same rotation.
for u in rad int turns deg; do
	check "rotate phase-$u"      0 "$BIN" "$M" "$WORK/d_$u.nii" -quiet \
		-phase "$WORK/phase-$u.nii" -real "$WORK/rot_$u.nii"
done
for u in int turns deg; do
	if "$GEN" cmp "$WORK/rot_rad.nii" "$WORK/rot_$u.nii" 0.01 >/dev/null 2>&1
	then ok "phase-$u agrees with phase-rad"
	else bad "phase-$u disagrees with phase-rad (max diff $("$GEN" cmp "$WORK/rot_rad.nii" "$WORK/rot_$u.nii" 0 2>/dev/null))"; fi
done
# Explicit units must reach the same answer as auto on full-turn data.
check "-phaseunits turns"        0 "$BIN" "$M" "$WORK/d_x.nii" -quiet \
	-phase "$WORK/phase-turns.nii" -phaseunits turns -real "$WORK/rot_x.nii"
if "$GEN" cmp "$WORK/rot_rad.nii" "$WORK/rot_x.nii" 0.01 >/dev/null 2>&1
then ok "explicit turns agrees with auto"; else bad "explicit turns disagrees with auto"; fi
# The rotation must actually do something: signed output, unlike a magnitude.
# `cmp` exits 2 for an unreadable or mis-sized file, so test for exactly 1
# ("readable, and differs") rather than "anything but 0".
"$GEN" cmp "$WORK/rot_rad.nii" "$M" 1.0 >/dev/null 2>&1; sv_rc=$?
if [ "$sv_rc" -eq 1 ]; then ok "rotation changed the data"
elif [ "$sv_rc" -eq 0 ]; then bad "rotation was a no-op (output still equals the magnitude)"
else bad "rotation comparison failed to run (cmp exit $sv_rc)"; fi

echo "numerical anchor"
# Everything above compares the tool against ITSELF under a different encoding, so
# a change to the arithmetic -- the binomial kernel, the reference-volume rule --
# appears on both sides and cancels.  Mutation testing confirmed that: altering a
# kernel tap broke nothing.  These two pin the actual numbers.  RMS rather than a
# byte comparison, at 1e-5 relative, so last-bit platform differences do not fire.
if "$GEN" rms "$WORK/rot_rad.nii" 507.123958 1e-5 >/dev/null
then ok "rotated series matches its recorded RMS"
else bad "rotated RMS moved: $("$GEN" rms "$WORK/rot_rad.nii" 0 0 2>/dev/null) vs 507.123958"; fi
if "$GEN" rms "$WORK/d_rad.nii" 402.549863 1e-5 >/dev/null
then ok "denoised series matches its recorded RMS"
else bad "denoised RMS moved: $("$GEN" rms "$WORK/d_rad.nii" 0 0 2>/dev/null) vs 402.549863"; fi

echo "no-op rotation"
# An explicit unit can still make every phase an exact multiple of a turn, which
# rotates by 1 everywhere and hands back the magnitude while reporting success.
# Detected from the phase itself (max |sin| == 0), not from the sign of the
# output: a high-SNR series can rotate correctly and stay entirely positive.
# Status and output captured separately from the message: the warning is emitted
# before the denoise, so piping the run into grep would let a later failure pass.
rm -f "$WORK/noop.nii"
"$BIN" "$M" "$WORK/noop.nii" -quiet -phase "$WORK/phase-siemens.nii" -phaseunits turns \
	>"$WORK/stdout" 2>"$WORK/stderr"; sv_rc=$?
if [ "$sv_rc" -ne 0 ]; then bad "no-op run failed (exit $sv_rc)"
elif [ ! -s "$WORK/noop.nii" ]; then bad "no-op run wrote no output"
elif grep -q "whole number of turns" "$WORK/stderr"; then ok "a no-op rotation is reported"
else bad "a no-op rotation passed silently"; fi
# ...but half-integer turns are exp(i*k*pi) = +/-1, a real sign modulation. A
# sin-only test passes that and would claim "the magnitude unchanged", which is
# false. The warning must stay silent here, and the rotation must really differ.
"$GEN" mk phase-halfturns "$WORK/phase-half.nii" || exit 1
# The exit status and the output are checked SEPARATELY from the message: piping
# into grep discards the denoiser's status, and treating "any nonzero comparator
# exit" as success would let exit 2 (unreadable/missing/mis-sized) print ok.
rm -f "$WORK/rot_half.nii"
"$BIN" "$M" "$WORK/o.nii" -quiet -phase "$WORK/phase-half.nii" -phaseunits turns \
	-real "$WORK/rot_half.nii" >"$WORK/stdout" 2>"$WORK/stderr"; sv_rc=$?
if [ "$sv_rc" -ne 0 ]; then bad "half-turn run failed (exit $sv_rc)"
elif [ ! -s "$WORK/rot_half.nii" ]; then bad "half-turn run wrote no rotated series"
elif grep -q "whole number of turns" "$WORK/stderr"; then
	bad "half-turn phase wrongly reported as a no-op"
else
	"$GEN" cmp "$WORK/rot_half.nii" "$M" 1.0 >/dev/null 2>&1; sv_rc=$?
	if [ "$sv_rc" -eq 1 ]; then ok "half-turn phase is a real rotation, not reported as a no-op"
	elif [ "$sv_rc" -eq 0 ]; then bad "half-turn phase did not actually change the data"
	else bad "half-turn comparison failed to run (cmp exit $sv_rc)"; fi
fi

# -noise and -rank go through the same nim-mutating writer as the main output and
# nothing else inspects them: either write could be deleted and the suite stayed
# green.  This run also pins dn_write_*'s "safe to call in any order" claim,
# which no other test reaches -- nothing else asks for more than two outputs.
check "all four outputs in one run" 0 "$BIN" "$M" "$WORK/all.nii" -quiet \
	-phase "$WORK/phase-rad.nii" -noise "$WORK/all_n.nii" \
	-rank "$WORK/all_r.nii" -real "$WORK/all_re.nii"
if [ -s "$WORK/all_n.nii" ] && [ -s "$WORK/all_r.nii" ] && [ -s "$WORK/all_re.nii" ]
then ok "-noise, -rank and -real were all written"
else bad "an auxiliary output was not written"; fi
"$GEN" cmp "$WORK/all.nii" "$WORK/d_rad.nii" 0 >/dev/null 2>&1; sv_rc=$?
if [ "$sv_rc" -eq 0 ]; then ok "the auxiliary outputs do not change the denoised series"
elif [ "$sv_rc" -eq 1 ]; then bad "asking for -noise/-rank/-real changed the denoised output"
else bad "auxiliary-output comparison failed to run (cmp exit $sv_rc)"; fi

echo "constant phase"
# A constant phase is a constant rotation, removed whatever its scale, so the
# rotated series must come back as exactly the magnitude.  This was once refused.
check "constant phase accepted"  0 "$BIN" "$M" "$WORK/d_c.nii" -quiet \
	-phase "$WORK/phase-const.nii" -real "$WORK/rot_c.nii"
if "$GEN" cmp "$WORK/rot_c.nii" "$M" 0 >/dev/null 2>&1
then ok "constant phase reproduces the magnitude exactly"
else bad "constant phase altered the magnitude (max diff $("$GEN" cmp "$WORK/rot_c.nii" "$M" 0 2>/dev/null))"; fi

echo "determinism"
# Byte-identical output across thread counts is a promise the project makes.
# The worker count is READ BACK, not assumed.  dn_effective_threads caps by the
# core count and by DN_CHUNK-sized chunks, so on this small fixture "-nthreads 8"
# already means 4 -- and on a one-core host it means 1, at which point the
# comparison below is a serial run against itself.  Measured: forcing
# dn_default_threads() to 1 left the whole suite green.
check "run, 1 thread"            0 "$BIN" "$M" "$WORK/t1.nii" -nthreads 1 -phase "$WORK/phase-rad.nii"
check "run, many threads"        0 "$BIN" "$M" "$WORK/t8.nii" -nthreads 8 -phase "$WORK/phase-rad.nii"
sv_n=$(sed -n 's/.*threads *: *//p' "$WORK/stderr")
if [ "${sv_n:-1}" -gt 1 ]
then ok "the multi-thread run really used $sv_n workers"
else bad "the multi-thread run used ${sv_n:-?} worker; the comparison below proves nothing"; fi
if cmp -s "$WORK/t1.nii" "$WORK/t8.nii"
then ok "1 thread == $sv_n threads, byte for byte"; else bad "thread count changed the output"; fi

echo
printf '%d passed, %d failed\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
