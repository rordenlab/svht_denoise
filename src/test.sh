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

# check_msg <name> <pattern> <cmd...>: exit 0 AND <pattern> on stderr.  Needed
# because stderr text is otherwise asserted by nothing: the -pF factor in the run
# summary, and the -phase warning, could both be deleted with this suite green.
check_msg() {
	name=$1; pat=$2; shift 2
	"$@" >"$WORK/stdout" 2>"$WORK/stderr"; got=$?
	if [ "$got" -eq 0 ] && grep -q "$pat" "$WORK/stderr"
	then ok "$name"
	else bad "$name (exit $got, or no '$pat' on stderr)"; sed 's/^/       /' "$WORK/stderr"; fi
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

# Relative L2 error of an output against its noise-free truth.
err() { "$GEN" rl2 "$1" "$2" 1 | awk '{print $2}'; }

for k in mag mask phase-rad phase-int phase-turns phase-deg phase-affine phase-const phase-siemens; do
	"$GEN" mk "$k" "$WORK/$k.nii" || { echo "fixture $k failed"; exit 1; }
done
M=$WORK/mag.nii
# The original pipeline.  Sections pinned before the defaults changed run with
# it, then override the one setting they test (the last value of an option wins).
LEGACY="-filter optthresh -aggregator exclusive -shape cube -vst n"

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
	check "rotate phase-$u"      0 "$BIN" "$M" "$WORK/d_$u.nii" $LEGACY -quiet \
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
check "all four outputs in one run" 0 "$BIN" "$M" "$WORK/all.nii" $LEGACY -quiet \
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

echo "compressed I/O"
# The .nii.gz path had NO coverage at all, which stopped being acceptable when
# the default macOS build switched from the system zlib to a statically linked
# zlib-ng: the whole suite passed either way, so nothing here would have noticed
# the compression library changing underneath it -- or failing to.
# The reference is "gz1plain", NOT "gz1.nii" beside "gz1.nii.gz", and the two
# must never share a basename.  nifti_findimgname() resolves the data file by
# stripping the extension and trying base+".nii" BEFORE base+".nii.gz", so with
# matching names the read below took its header from the .gz and its voxels from
# the plain sibling -- the comparison then ran two passes over identical bytes.
# Measured: with the names matched, corrupting a byte in gzread OR in gzwrite
# left the whole suite green; with them distinct, both are caught.
check "write .nii.gz"           0 "$BIN" "$M" "$WORK/gz1.nii.gz" -quiet
check "write .nii for reference" 0 "$BIN" "$M" "$WORK/gz1plain.nii"   -quiet
# Exit 0 is NOT enough, and this is not hypothetical: built without HAVE_ZLIB
# the tool returns success having created no .gz at all -- nifti_io rewrites the
# name and writes the uncompressed file beside it.  Assert the requested path
# exists before asserting anything about its contents.
if [ -f "$WORK/gz1.nii.gz" ]
then ok "the .gz output file exists"
else bad "exit 0 but no .gz file was created"; fi
# Actually a gzip stream, not an uncompressed file wearing a .gz name.
if [ "$(od -An -tx1 -N2 < "$WORK/gz1.nii.gz" | tr -d ' \n')" = "1f8b" ]
then ok "the .gz output is a real gzip stream"
else bad "the .gz output has no gzip magic"; fi
if [ "$(wc -c < "$WORK/gz1.nii.gz")" -lt "$(wc -c < "$WORK/gz1plain.nii")" ]
then ok "the .gz output is smaller than the uncompressed one"
else bad "the .gz output did not actually compress"; fi
# Read both back THROUGH THE TOOL, so this needs no gzip(1) and stays honest on
# any host: if the write or the read dropped or altered a byte, the two
# second-pass results diverge.  Comparing the .gz against the .nii directly
# would not work -- they are different file formats holding the same data.
check "read .nii.gz back in"    0 "$BIN" "$WORK/gz1.nii.gz" "$WORK/gz2_c.nii" -quiet
check "read .nii back in"       0 "$BIN" "$WORK/gz1plain.nii"    "$WORK/gz2_u.nii" -quiet
if "$GEN" cmp "$WORK/gz2_c.nii" "$WORK/gz2_u.nii" 0 >/dev/null 2>&1
then ok "compressed and uncompressed round trips agree exactly"
else bad "the .gz round trip changed the data (max diff $("$GEN" cmp "$WORK/gz2_c.nii" "$WORK/gz2_u.nii" 0 2>/dev/null))"; fi

echo "split eigensolver path"
# Every fixture above is 8 volumes, which is below DN_SPLIT_MIN (48), so all of
# them take the FULL solver.  This one is 48 -- the smallest count that takes the
# split path: values-only tql2, inverse iteration on the tridiagonal, and the
# back-transform through the stored reflectors.  Measured with llvm-cov: on the
# 8-volume fixtures tred_reduce, tridiag_solve, tridiag_eigvec and tred_apply_q
# have ZERO coverage, and on this one they have 81-92% while tred2 has none.
# The two fixtures are complementary, which is why both stay.
"$GEN" mk big48 "$WORK/big48.nii" || { echo "fixture big48 failed"; exit 1; }
check "48-volume run" 0 "$BIN" "$WORK/big48.nii" "$WORK/b48.nii" $LEGACY -quiet -rank "$WORK/b48_r.nii"
# Ranks are integers, so a NON-INTEGER rank-map RMS is by itself proof that the
# map varies.  That matters: a run that retained nothing anywhere would exit 0
# with dn_eig_vectors never entered, and one constant rank would satisfy any
# anchor that still produced a constant.  Here the rank runs 3 to 7.
if "$GEN" rms "$WORK/b48_r.nii" 5.74348799 1e-5 >/dev/null
then ok "split-path rank map matches its recorded RMS"
else bad "split-path rank map moved: $("$GEN" rms "$WORK/b48_r.nii" 0 0 2>/dev/null) vs 5.74348799"; fi
if "$GEN" rms "$WORK/b48.nii" 1061.54322 1e-5 >/dev/null
then ok "split-path denoised series matches its recorded RMS"
else bad "split-path denoised RMS moved: $("$GEN" rms "$WORK/b48.nii" 0 0 2>/dev/null) vs 1061.54322"; fi
# The byte-identity promise, on the path where the eigenvectors are built one at
# a time and re-orthogonalised against their predecessors.  The 8-volume case
# below never reaches that code.
check "48-volume, 1 thread"      0 "$BIN" "$WORK/big48.nii" "$WORK/b48_t1.nii" $LEGACY -nthreads 1 -quiet
# NOT -quiet, and the worker count is read back, for the same reason the
# determinism section below does it: dn_effective_threads caps by DN_CHUNK-sized
# chunks and by the core count, so on a one-core host "-nthreads 8" means 1 and
# the comparison below would be a serial run against itself, reporting success.
check "48-volume, many threads"  0 "$BIN" "$WORK/big48.nii" "$WORK/b48_t8.nii" $LEGACY -nthreads 8
sv_n=$(sed -n 's/.*threads *: *//p' "$WORK/stderr")
if [ "${sv_n:-1}" -gt 1 ]
then ok "the multi-thread split run really used $sv_n workers"
else bad "the split run used ${sv_n:-?} worker; the comparison below proves nothing"; fi
if cmp -s "$WORK/b48_t1.nii" "$WORK/b48_t8.nii"
then ok "split path is byte-identical across thread counts"; else bad "split path output depends on the thread count"; fi

# dn_testgen's rl2 is the metric an optimisation is judged on, so its own
# arithmetic needs pinning: it is twenty lines that could silently start
# returning zero, and every use of it is a comparison it would then pass.
if "$GEN" rl2 "$WORK/b48_t1.nii" "$WORK/b48_t8.nii" 0 >/dev/null 2>&1
then ok "rl2 reports zero for identical images"
else bad "rl2 reported a difference between identical images"; fi
"$GEN" rl2 "$WORK/b48.nii" "$WORK/big48.nii" 1e-9 >/dev/null 2>&1; sv_rc=$?
if [ "$sv_rc" -eq 1 ]; then ok "rl2 detects a real difference"
elif [ "$sv_rc" -eq 0 ]; then bad "rl2 called a denoised series identical to its input"
else bad "rl2 failed to run (exit $sv_rc)"; fi

echo "determinism"
# Byte-identical output across thread counts is a promise the project makes.
# The worker count is READ BACK, not assumed.  dn_effective_threads caps by the
# core count and by DN_CHUNK-sized chunks, so on this small fixture "-nthreads 8"
# already means 4 -- and on a one-core host it means 1, at which point the
# comparison below is a serial run against itself.  Measured: forcing
# dn_default_threads() to 1 left the whole suite green.
check "run, 1 thread"            0 "$BIN" "$M" "$WORK/t1.nii" $LEGACY -nthreads 1 -phase "$WORK/phase-rad.nii"
check "run, many threads"        0 "$BIN" "$M" "$WORK/t8.nii" $LEGACY -nthreads 8 -phase "$WORK/phase-rad.nii"
sv_n=$(sed -n 's/.*threads *: *//p' "$WORK/stderr")
if [ "${sv_n:-1}" -gt 1 ]
then ok "the multi-thread run really used $sv_n workers"
else bad "the multi-thread run used ${sv_n:-?} worker; the comparison below proves nothing"; fi
if cmp -s "$WORK/t1.nii" "$WORK/t8.nii"
then ok "1 thread == $sv_n threads, byte for byte"; else bad "thread count changed the output"; fi

# -degibbs is a build option, so ask the binary rather than assuming.  A
# DEGIBBS=0 build must run the rest of this suite green and skip only this.
#
# Ask the BINARY WHAT IT DOES, not what its help text looks like.  Grepping -help
# for "-degibbs" once matched prose in the -phase entry on a DEGIBBS=0 build, ran
# this whole block against a binary that refuses every option in it, and reported
# 27 failures; grepping for the option entry instead would still break on a
# rewrap.  A rejected value is answered either "must be y, n or o" or "not
# compiled into this build", and only the second means the feature is absent.
if ! "$BIN" -degibbs z 2>&1 | grep -q 'not compiled into this build'; then
echo "degibbs"
# "only" does no denoising, so the denoiser's own options have nothing to act on.
# Ignoring them silently would write fewer files than were asked for, at exit 0.
check "-degibbs o rejects -noise" 1 "$BIN" "$M" "$WORK/o.nii" -quiet -degibbs o -noise "$WORK/n.nii"
check "-degibbs o rejects -mask"  1 "$BIN" "$M" "$WORK/o.nii" -quiet -degibbs o -mask "$WORK/mask.nii"
check "bad -degibbs value"        1 "$BIN" "$M" "$WORK/o.nii" -quiet -degibbs z
# -mask with -degibbs y once wrote NON-ZERO data outside the mask -- 647 of 648
# voxels, peaking at a quarter of the data range -- because ringing removal ran
# over the hard zero edge the mask leaves, which also rings INWARDS into voxels
# that are inside it.  Refused rather than re-zeroed: re-zeroing would restore
# the invariant this suite checks while leaving the in-mask corruption.
check "-mask refused with -degibbs y" 1 "$BIN" "$M" "$WORK/o.nii" -quiet -mask "$WORK/mask.nii" -degibbs y
check "-degibbs y runs"           0 "$BIN" "$M" "$WORK/dgy.nii" -quiet -degibbs y
# ...and it must HONOUR -pF.  Every other -pF test uses "only", so dropping the
# factor on the denoise-then-degibbs path -- the one a real user runs -- left the
# whole suite green.
check "-degibbs y honours -pF"    0 "$BIN" "$M" "$WORK/dgy78.nii" -quiet -degibbs y -pF 0.875
"$GEN" cmp "$WORK/dgy78.nii" "$WORK/dgy.nii" 0 >/dev/null 2>&1; sv_rc=$?
if [ "$sv_rc" -eq 1 ]; then ok "-degibbs y + -pF differs from -degibbs y"
elif [ "$sv_rc" -eq 0 ]; then bad "-pF was ignored on the -degibbs y path"
else bad "comparison failed to run (cmp exit $sv_rc)"; fi
# The run summary names the factor, and the -phase warning survives -quiet --
# both are stderr-only, so nothing else here would notice them being deleted.
check_msg "the summary names the factor" "partial Fourier 7/8" \
	"$BIN" "$M" "$WORK/rep.nii" -degibbs o -pF 0.875
check_msg "-phase + -degibbs y warns, even under -quiet" "truncates it to zero" \
	"$BIN" "$M" "$WORK/warn.nii" -quiet -phase "$WORK/phase-rad.nii" -degibbs y
# -pF changes how ringing is removed, so it is meaningless without a stage that
# removes any; and the pipeline differs PER factor rather than varying
# continuously, so an unimplemented factor must be refused, never rounded to the
# nearest one we do have.
check "-pF needs -degibbs"        1 "$BIN" "$M" "$WORK/o.nii" -quiet -pF 0.75
# 1.0 is a valid factor and still means nothing without -degibbs.  It was the one
# value accepted there, because the guard tested the factor instead of whether
# the option was given.
check "-pF 1.0 needs -degibbs"    1 "$BIN" "$M" "$WORK/o.nii" -quiet -pF 1.0
check "-pF 7/8 runs"              0 "$BIN" "$M" "$WORK/pf875.nii" -quiet -degibbs o -pF 0.875
# PINNED, for the same reason the full-Fourier path is: exit 0 says nothing about
# the arithmetic, and a 7/8 run that silently returned the conventional answer
# passed every check here before this line existed.
check "-pF 7/8 RMS is unchanged"  0 "$GEN" rms "$WORK/pf875.nii" 614.525105 1e-6
# The full-Fourier answer for the same fixture, made HERE rather than reused from
# further down: it is both the "-pF absent changes nothing" anchor and the thing
# 7/8 has to differ from, and a file produced later would order this wrongly.
check "no -pF still runs"         0 "$BIN" "$M" "$WORK/dgnopf.nii" -quiet -degibbs o
"$GEN" cmp "$WORK/pf875.nii" "$WORK/dgnopf.nii" 0 >/dev/null 2>&1; sv_rc=$?
if [ "$sv_rc" -eq 1 ]; then ok "-pF 7/8 differs from full k-space"
elif [ "$sv_rc" -eq 0 ]; then bad "-pF 7/8 returned the full-Fourier answer"
else bad "7/8 comparison failed to run (cmp exit $sv_rc)"; fi
check "-pF rejects out of range"  1 "$BIN" "$M" "$WORK/o.nii" -quiet -degibbs o -pF 1.5
check "-pF rejects junk"          1 "$BIN" "$M" "$WORK/o.nii" -quiet -degibbs o -pF abc
# In range, plausible, and NOT implemented.  Same branch and same message as 1.5
# -- there is no separate range check, deliberately, since dn_degibbs_check owns
# which factors exist.  Both are kept because they are different value CLASSES a
# user actually types, not because the refusals differ.
check "-pF rejects unimplemented" 1 "$BIN" "$M" "$WORK/o.nii" -quiet -degibbs o -pF 0.9
# The 9x9x3 fixture has an ODD y dimension, which the 6/8 pipeline cannot split
# into odd and even columns of equal length.
check "-pF refuses odd y"         1 "$BIN" "$M" "$WORK/o.nii" -quiet -degibbs o -pF 0.75
# The OTHER geometry refusal: y is big enough for the plain method but its -pF
# interleave is not.  A distinct branch from the odd-y and minimum-dimension ones
# above, reachable with real data at either factor, and reached by no other
# fixture -- 9x6x3 clears DG_MIN_DIM and then splits to 3 (6/8) and 4 (7/8).
"$GEN" mk mag-shorty "$WORK/mag-shorty.nii" || exit 1
check "-pF 6/8 refuses a short interleave" 1 "$BIN" "$WORK/mag-shorty.nii" "$WORK/o.nii" -quiet -degibbs o -pF 0.75
check "-pF 7/8 refuses a short interleave" 1 "$BIN" "$WORK/mag-shorty.nii" "$WORK/o.nii" -quiet -degibbs o -pF 0.875
# ...and the same fixture must still degibbs fine at full k-space, so the
# refusals above are about the interleave and not about the fixture.
check "short-y fixture runs full k-space"  0 "$BIN" "$WORK/mag-shorty.nii" "$WORK/shorty.nii" -quiet -degibbs o
# ...so 6/8 needs an even-y fixture, and without one the whole 6/8 SUCCESS
# path -- half-length axis, strided sub-image -- was executed by no test at all.
"$GEN" mk mag-even "$WORK/mag-even.nii" || exit 1
check "-pF 6/8 runs"              0 "$BIN" "$WORK/mag-even.nii" "$WORK/pf75.nii" -quiet -degibbs o -pF 0.75
check "-pF 6/8 RMS is unchanged"  0 "$GEN" rms "$WORK/pf75.nii" 421.277497 1e-6
check "even-y full k-space runs"  0 "$BIN" "$WORK/mag-even.nii" "$WORK/e_full.nii" -quiet -degibbs o
"$GEN" cmp "$WORK/pf75.nii" "$WORK/e_full.nii" 0 >/dev/null 2>&1; sv_rc=$?
if [ "$sv_rc" -eq 1 ]; then ok "-pF 6/8 differs from full k-space"
elif [ "$sv_rc" -eq 0 ]; then bad "-pF 6/8 returned the full-Fourier answer"
else bad "6/8 comparison failed to run (cmp exit $sv_rc)"; fi
# The method assumes MAGNITUDE data, so negative input is TRUNCATED on the way
# in -- truncated, not folded to +|v|, which would turn zero-mean noise into a
# Rician-like floor.  The witness is byte-identity: degibbs of a signed field
# must equal degibbs of the same field with the truncation already applied.  It
# is a no-op on magnitude input, which is why every pin above is untouched by it.
"$GEN" mk mag-signed  "$WORK/mag-signed.nii"  || exit 1
"$GEN" mk mag-clamped "$WORK/mag-clamped.nii" || exit 1
check "signed input runs"         0 "$BIN" "$WORK/mag-signed.nii"  "$WORK/sg_full.nii" -quiet -degibbs o
check "pre-clamped input runs"    0 "$BIN" "$WORK/mag-clamped.nii" "$WORK/cl_full.nii" -quiet -degibbs o
if cmp -s "$WORK/sg_full.nii" "$WORK/cl_full.nii"
then ok "negative input is truncated before unringing"; else bad "signed and pre-clamped input gave different output"; fi
check "signed input runs at 7/8"  0 "$BIN" "$WORK/mag-signed.nii"  "$WORK/sg78.nii" -quiet -degibbs o -pF 0.875
check "pre-clamped input at 7/8"  0 "$BIN" "$WORK/mag-clamped.nii" "$WORK/cl78.nii" -quiet -degibbs o -pF 0.875
if cmp -s "$WORK/sg78.nii" "$WORK/cl78.nii"
then ok "7/8 truncates negative input too"; else bad "7/8 signed and pre-clamped input differ"; fi
check "signed 7/8 RMS is unchanged" 0 "$GEN" rms "$WORK/sg78.nii" 413.200122 1e-6
# Clamping the INPUT does not make the OUTPUT non-negative, and must not be
# mistaken for doing so: ringing correction undershoots, and the conventional
# path keeps that undershoot exactly as mrdegibbs does.  7/8 is the one path that
# also clamps its output, matching its own reference.
check "conventional output may go negative" 0 "$GEN" neg "$WORK/sg_full.nii"
check "7/8 output has no negatives" 1 "$GEN" neg "$WORK/sg78.nii"
# "only" is the one mode that takes a 3D image: the >= 2 volume rule belongs to
# the denoiser, which is not running.  The mask fixture is the single-volume one.
check "-degibbs o accepts 3D"     0 "$BIN" "$WORK/mask.nii" "$WORK/dg3d.nii" -quiet -degibbs o
# "no" must be EXACTLY the build without the option; anything else means adding a
# second algorithm has disturbed the default path.
check "-degibbs n runs"           0 "$BIN" "$M" "$WORK/dgn.nii" -quiet -degibbs n
check "no -degibbs runs"          0 "$BIN" "$M" "$WORK/dgoff.nii" -quiet
if cmp -s "$WORK/dgn.nii" "$WORK/dgoff.nii"
then ok "-degibbs n == no -degibbs, byte for byte"; else bad "-degibbs n disturbed the default path"; fi
# Its own pool, so its own byte-identity check: the denoiser's above proves
# nothing about this one.
check "-degibbs o, 1 thread"      0 "$BIN" "$M" "$WORK/dg_t1.nii" -quiet -degibbs o -nthreads 1
check "-degibbs o, many threads"  0 "$BIN" "$M" "$WORK/dg_t8.nii" -quiet -degibbs o -nthreads 8
if cmp -s "$WORK/dg_t1.nii" "$WORK/dg_t8.nii"
then ok "degibbs: 1 thread == 8 threads, byte for byte"; else bad "degibbs thread count changed the output"; fi
# A PINNED number, because nothing else here would notice a changed twiddle: the
# kernel stays self-consistent under any coefficient, so only a recorded value
# fails when the arithmetic moves.  This one is not merely recorded from our own
# output -- mrdegibbs produces the same image BIT FOR BIT on this fixture, and
# the same RMS to the tolerance below.  The fixture is 9x9, deliberately ODD,
# which is the branch that skips the Nyquist zeroing.
check "degibbs RMS is unchanged"  0 "$GEN" rms "$WORK/dg_t1.nii" 614.68488 1e-6
# ...and it has to actually do something; an identity transform would pin fine.
"$GEN" cmp "$WORK/dg_t1.nii" "$M" 0 >/dev/null 2>&1; sv_rc=$?
if [ "$sv_rc" -eq 1 ]; then ok "-degibbs o changed the image"
elif [ "$sv_rc" -eq 0 ]; then bad "-degibbs o returned its input unchanged"
else bad "degibbs comparison failed to run (cmp exit $sv_rc)"; fi
fi

echo "per-shell demeaning"
# The sidecar is found by name.  Shells: b=0 x2, 1000 x3, 2000 x2; the lone 3000
# is a shell of one and must NOT be demeaned, so G = 3.
sv_d=$WORK/dm; mkdir -p "$sv_d"; cp "$M" "$sv_d/mag.nii"
echo "0 5 995 1000 1005 2000 2010 3000" > "$sv_d/mag.bval"
check_msg "-demean y finds the sidecar, clusters 3 shells" "b=1000 x3 b=2005 x2" \
	"$BIN" "$sv_d/mag.nii" "$WORK/dm.nii" $LEGACY -demean y -nthreads 1 -noise "$WORK/dm_s.nii"
check "-demean y, many threads"   0 "$BIN" "$sv_d/mag.nii" "$WORK/dm8.nii" $LEGACY -quiet -demean y -nthreads 8
if cmp -s "$WORK/dm.nii" "$WORK/dm8.nii"
then ok "demean: 1 thread == 8 threads, byte for byte"; else bad "demean: thread count changed the output"; fi
# Pinned: beta must count N-G columns, and the means must be restored.  Getting
# either wrong moves this, where a self-consistency check would not notice.
check "demeaned RMS is unchanged" 0 "$GEN" rms "$WORK/dm.nii" 614.498087 1e-6
# The output alone does not notice beta counting N instead of N-G on this
# fixture (mutation-tested); sigma depends on beta directly, so pin it too.
check "demeaned noise RMS is unchanged" 0 "$GEN" rms "$WORK/dm_s.nii" 30.3373923 1e-6
check "-demean n ignores a sidecar" 0 "$BIN" "$sv_d/mag.nii" "$WORK/dmn.nii" $LEGACY -quiet -demean n
"$BIN" "$M" "$WORK/plain.nii" $LEGACY -quiet
if cmp -s "$WORK/dmn.nii" "$WORK/plain.nii"
then ok "-demean n == no demeaning, byte for byte"; else bad "-demean n changed the output"; fi
echo "0 0 1000" > "$sv_d/short.bval"
# MRtrix3 shares a shell only BELOW 80 s/mm^2 apart, so 1000 and 1080 split.
echo "0 0 1000 1000 1080 1080 2000 2000" > "$sv_d/edge.bval"
check_msg "a gap of exactly 80 starts a shell" "4 shells" \
	"$BIN" "$sv_d/mag.nii" "$WORK/edge.nii" $LEGACY -demean y -bval "$sv_d/edge.bval"
check "-bval with the wrong count"  1 "$BIN" "$sv_d/mag.nii" "$WORK/x.nii" -quiet -demean y -bval "$sv_d/short.bval"
check "-demean y with no b-values"  1 "$BIN" "$M" "$WORK/x.nii" -quiet -demean y
check "-bval without -demean y"     1 "$BIN" "$sv_d/mag.nii" "$WORK/x.nii" -quiet -bval "$sv_d/mag.bval"

echo "sigma-based filters"
# big48 is low rank plus noise ramped across x, and big48-clean is its truth, so
# a filter's error can be measured rather than only pinned.
"$GEN" mk big48 "$WORK/f48.nii" && "$GEN" mk big48-clean "$WORK/f48c.nii" || { echo "fixture big48-clean failed"; exit 1; }
for sv_f in optthresh truncate optshrink; do
	check "-filter $sv_f runs" 0 "$BIN" "$WORK/f48.nii" "$WORK/f_$sv_f.nii" $LEGACY -quiet -filter "$sv_f" \
		-nthreads 1 -noise "$WORK/fs_$sv_f.nii"
done
# The point of shrinkage: closer to the truth than the hard threshold it replaces.
sv_e1=$(err "$WORK/f_optthresh.nii" "$WORK/f48c.nii")
sv_e2=$(err "$WORK/f_optshrink.nii" "$WORK/f48c.nii")
if awk "BEGIN{exit !($sv_e2 < $sv_e1)}"
then ok "optshrink beats optthresh against the truth ($sv_e2 < $sv_e1)"
else bad "optshrink error $sv_e2 is not below optthresh's $sv_e1"; fi
# Pinned, because the comparison above survives many wrong shrinkers.  The noise
# map is the pass-1 estimate, shared by both sigma-based filters.
check "optshrink RMS is unchanged"  0 "$GEN" rms "$WORK/f_optshrink.nii" 1050.76197 1e-6
check "truncate RMS is unchanged"   0 "$GEN" rms "$WORK/f_truncate.nii" 1076.92601 1e-6
check "sigma map RMS is unchanged"  0 "$GEN" rms "$WORK/fs_optshrink.nii" 326.854857 1e-6
check "-filter optshrink, many threads" 0 "$BIN" "$WORK/f48.nii" "$WORK/f_os8.nii" $LEGACY -quiet -filter optshrink -nthreads 8
if cmp -s "$WORK/f_optshrink.nii" "$WORK/f_os8.nii"
then ok "optshrink: 1 thread == 8 threads, byte for byte"; else bad "optshrink: thread count changed the output"; fi
check "-filter with a bad name"      1 "$BIN" "$M" "$WORK/x.nii" -quiet -filter median

echo "spherical patch"
# Measured in mm: on 1 x 1 x 2.5 mm voxels the 27-voxel ball reaches two voxels
# in-plane but one through-plane, so it closes at sqrt(2.5^2 + 1) with 31.
"$GEN" mk mag-aniso "$WORK/aniso.nii" || { echo "fixture mag-aniso failed"; exit 1; }
check_msg "sphere on isotropic voxels"   "radius 1.73 mm = 27 voxels" "$BIN" "$M" "$WORK/sp.nii" $LEGACY -shape sphere -nthreads 1
check_msg "sphere on anisotropic voxels" "radius 2.69 mm = 31 voxels" "$BIN" "$WORK/aniso.nii" "$WORK/spa.nii" $LEGACY -shape sphere
check "sphere RMS is unchanged"          0 "$GEN" rms "$WORK/sp.nii" 614.250794 1e-6
check "sphere on 8 threads"              0 "$BIN" "$M" "$WORK/sp8.nii" $LEGACY -quiet -shape sphere -nthreads 8
if cmp -s "$WORK/sp.nii" "$WORK/sp8.nii"
then ok "sphere: 1 thread == 8 threads, byte for byte"; else bad "sphere: thread count changed the output"; fi
check "sphere + optshrink on big48"      0 "$BIN" "$WORK/f48.nii" "$WORK/spf.nii" $LEGACY -quiet -shape sphere \
	-filter optshrink -nthreads 1 -noise "$WORK/spf_s.nii"
check "sphere optshrink RMS is unchanged" 0 "$GEN" rms "$WORK/spf.nii" 1050.45481 1e-6
check "sphere sigma map RMS is unchanged" 0 "$GEN" rms "$WORK/spf_s.nii" 330.26322 1e-6

echo "overlapping-patch averaging"
# Pinned configurations on the 9x9 fixture; threading is tested on mag-wide below.
# Neither proves the tile width safe -- see AGENTS.md on why tsan-test does not.
check "gaussian, sphere, stride 2" 0 "$BIN" "$M" "$WORK/ag1.nii" $LEGACY -quiet -aggregator gaussian -shape sphere -stride 2
check "uniform, stride 2"          0 "$BIN" "$M" "$WORK/agu.nii" $LEGACY -quiet -aggregator uniform -stride 2
# Every 9x9 fixture is one tile per colour, so runs on it are serial at any
# -nthreads.  mag-wide is wide enough for several; the worker count is read back.
"$GEN" mk mag-wide "$WORK/wide.nii" || { echo "fixture mag-wide failed"; exit 1; }
check "defaults on a wide image, 1 thread" 0 "$BIN" "$WORK/wide.nii" "$WORK/w1.nii" -quiet -nthreads 1
check "defaults on a wide image, 8 threads" 0 "$BIN" "$WORK/wide.nii" "$WORK/w8.nii" -nthreads 8
sv_n=$(sed -n 's/.*threads *: *//p' "$WORK/stderr")
if [ "${sv_n:-1}" -gt 1 ]
then ok "the averaging run really used $sv_n workers"
else bad "the averaging run used ${sv_n:-?} worker; the comparison below proves nothing"; fi
if cmp -s "$WORK/w1.nii" "$WORK/w8.nii"
then ok "averaging: 1 thread == $sv_n threads, byte for byte"; else bad "averaging: thread count changed the output"; fi
check "sphere stride-2 RMS is unchanged" 0 "$GEN" rms "$WORK/ag1.nii" 614.241667 1e-6
check "gaussian stride 2 on big48" 0 "$BIN" "$WORK/f48.nii" "$WORK/ag48.nii" $LEGACY -quiet -aggregator gaussian \
	-stride 2 -filter optshrink -nthreads 1 -rank "$WORK/ag48_r.nii" -noise "$WORK/ag48_s.nii"
sv_e3=$(err "$WORK/ag48.nii" "$WORK/f48c.nii")
if awk "BEGIN{exit !($sv_e3 < $sv_e2)}"
then ok "averaging beats one patch per voxel against the truth ($sv_e3 < $sv_e2)"
else bad "averaged error $sv_e3 is not below exclusive optshrink's $sv_e2"; fi
check "averaged RMS is unchanged"      0 "$GEN" rms "$WORK/ag48.nii" 1049.44369 1e-6
check "averaged rank RMS is unchanged" 0 "$GEN" rms "$WORK/ag48_r.nii" 8.37028843 1e-6
check "averaged noise RMS is unchanged" 0 "$GEN" rms "$WORK/ag48_s.nii" 326.568655 1e-6
check "-stride 2 refused under exclusive" 1 "$BIN" "$M" "$WORK/x.nii" -quiet -aggregator exclusive -stride 2
# Off-centre weights a float accumulator cannot carry once counted as coverage and
# zeroed voxels at exit 0: a tiny FWHM, or 20 mm slices at stride 2.  Both must now
# fall back to the voxel's own patch.
check "-aggregator_fwhm 0.07"          0 "$BIN" "$M" "$WORK/fw07.nii" -quiet -aggregator_fwhm 0.07
check "tiny FWHM zeroes no voxel"      0 "$GEN" nozero "$WORK/fw07.nii"
"$GEN" mk mag-slab "$WORK/slab.nii" || { echo "fixture mag-slab failed"; exit 1; }
check "20 mm slices, stride 2"         0 "$BIN" "$WORK/slab.nii" "$WORK/slab_o.nii" -quiet -aggregator_fwhm 0.5 -nthreads 1
check "20 mm slices zero no voxel"     0 "$GEN" nozero "$WORK/slab_o.nii"
# ...and weight TIMES value: data near 1e-16 on 11 mm slices underflowed float
# even where the weight alone did not, so no weight cutoff could catch it.
"$GEN" mk mag-tiny "$WORK/tiny.nii" || { echo "fixture mag-tiny failed"; exit 1; }
check "data near 1e-16, stride 2"      0 "$BIN" "$WORK/tiny.nii" "$WORK/tiny_o.nii" -quiet -shape cube -vst n -aggregator_fwhm 0.5
check "tiny data zeroes no voxel"      0 "$GEN" nozero "$WORK/tiny_o.nii"
# The -noise map is a running mean too; summed, its weight*sigma underflowed even
# in double at FWHM 0.156 (weights near 1e-304) and wrote 0 for covered voxels.
check "tiny data, -noise, FWHM 0.156"  0 "$BIN" "$WORK/tiny.nii" "$WORK/tiny_o2.nii" -quiet -shape cube -vst n \
	-aggregator_fwhm 0.156 -noise "$WORK/tiny_n.nii"
check "tiny data's noise map has no zero" 0 "$GEN" nozero "$WORK/tiny_n.nii"
check "20 mm slices RMS is unchanged"  0 "$GEN" rms "$WORK/slab_o.nii" 613.319384 1e-6
# Each option must do something: pinned, and different from what it replaces.
check "uniform RMS is unchanged"       0 "$GEN" rms "$WORK/agu.nii" 614.244489 1e-6
check "-aggregator_fwhm 1"             0 "$BIN" "$M" "$WORK/fw1.nii" -quiet -aggregator_fwhm 1
check "FWHM 1 RMS is unchanged"        0 "$GEN" rms "$WORK/fw1.nii" 613.309137 1e-6
check "-aggregator_fwhm needs gaussian"   1 "$BIN" "$M" "$WORK/x.nii" -quiet -aggregator uniform -aggregator_fwhm 3
check "-stride 3 is refused"              1 "$BIN" "$M" "$WORK/x.nii" -quiet -aggregator gaussian -stride 3

echo "variance stabilisation"
# Rician noise of sigma 100 on a signal at SNR 0.5-3, against its noise-free
# truth: the noise floor dominates, so a debiased output must score better.
"$GEN" mk rician "$WORK/ric.nii" && "$GEN" mk rician-clean "$WORK/ricc.nii" || { echo "fixture rician failed"; exit 1; }
check "-noise_dof 2 on Rician data"    0 "$BIN" "$WORK/ric.nii" "$WORK/dof2.nii" -quiet -noise_dof 2
check "-noise_dof 2 RMS is unchanged"  0 "$GEN" rms "$WORK/dof2.nii" 88.060271 1e-6
# A channel count far above the truth puts every sample under the modelled floor,
# and the unbiased inverse then returns 0 everywhere: that must not pass silently.
check_msg "a wrong -noise_dof is warned about" "of the output is zero" \
	"$BIN" "$WORK/ric.nii" "$WORK/dof4.nii" -noise_dof 4
check "magnitude, no -vst"   0 "$BIN" "$WORK/ric.nii" "$WORK/v0.nii" $LEGACY -quiet -filter optshrink
check "-vst y"               0 "$BIN" "$WORK/ric.nii" "$WORK/v1.nii" $LEGACY -quiet -filter optshrink -vst y \
	-nthreads 1 -noise "$WORK/v1_s.nii"
check "-vst y, keeping the bias" 0 "$BIN" "$WORK/ric.nii" "$WORK/v2.nii" $LEGACY -quiet -filter optshrink -vst y -preserve_noise_bias
sv_v0=$(err "$WORK/v0.nii" "$WORK/ricc.nii")
sv_v1=$(err "$WORK/v1.nii" "$WORK/ricc.nii")
sv_v2=$(err "$WORK/v2.nii" "$WORK/ricc.nii")
if awk "BEGIN{exit !($sv_v1 < $sv_v2 && $sv_v2 < $sv_v0)}"
then ok "debiased < biased VST < no VST against the truth ($sv_v1 < $sv_v2 < $sv_v0)"
else bad "VST ordering broken: debiased $sv_v1, biased $sv_v2, none $sv_v0"; fi
# The map must come out near the true 100: raw magnitudes alone read ~91.
check "VST output RMS is unchanged" 0 "$GEN" rms "$WORK/v1.nii" 182.797846 1e-6
check "VST sigma RMS is unchanged"  0 "$GEN" rms "$WORK/v1_s.nii" 101.459479 1e-6
check "-vst y, many threads" 0 "$BIN" "$WORK/ric.nii" "$WORK/v8.nii" $LEGACY -quiet -filter optshrink -vst y -nthreads 8
if cmp -s "$WORK/v1.nii" "$WORK/v8.nii"
then ok "vst: 1 thread == 8 threads, byte for byte"; else bad "vst: thread count changed the output"; fi
check "-vst y with -phase is refused"       1 "$BIN" "$M" "$WORK/x.nii" -quiet -vst y -phase "$WORK/phase-rad.nii"
check "-noise_dof with -vst n is refused"   1 "$BIN" "$M" "$WORK/x.nii" -quiet -vst n -noise_dof 2
check "-preserve_noise_bias with -vst n"    1 "$BIN" "$M" "$WORK/x.nii" -quiet -vst n -preserve_noise_bias

echo "defaults"
# The defaults are sphere + optshrink + gaussian stride 2 + VST.  Everything above
# pins each piece against the original pipeline; these pin the combination.
check "defaults on big48"        0 "$BIN" "$WORK/f48.nii" "$WORK/def.nii" -quiet -nthreads 1 \
	-noise "$WORK/def_s.nii" -rank "$WORK/def_r.nii"
# big48's noise is uniform and goes negative, so the VST stays off here.
check "defaults RMS is unchanged"       0 "$GEN" rms "$WORK/def.nii" 1049.27494 1e-6
check "defaults noise RMS is unchanged" 0 "$GEN" rms "$WORK/def_s.nii" 329.844936 1e-6
check "defaults rank RMS is unchanged"  0 "$GEN" rms "$WORK/def_r.nii" 8.90914359 1e-6
check "defaults, no side outputs" 0 "$BIN" "$WORK/f48.nii" "$WORK/def_plain.nii" -quiet -nthreads 8
if cmp -s "$WORK/def.nii" "$WORK/def_plain.nii"
then ok "defaults: -noise/-rank and the thread count leave the output alone"
else bad "defaults: -noise/-rank or the thread count changed the output"; fi
check "defaults on rician"       0 "$BIN" "$WORK/ric.nii" "$WORK/def_ric.nii" -quiet
sv_d1=$(err "$WORK/def_ric.nii" "$WORK/ricc.nii")
if awk "BEGIN{exit !($sv_d1 < $sv_v0)}"
then ok "defaults beat the old pipeline's optshrink on Rician data ($sv_d1 < $sv_v0)"
else bad "defaults error $sv_d1 is not below $sv_v0"; fi
# -help promises a masked run equals an unmasked one inside the mask.  Under
# averaging that holds only because a centre is kept whenever its patch touches
# the mask; keeping centres by any other rule breaks it at the mask edge.
check "masked == unmasked inside the mask" 0 "$GEN" cmpmask "$WORK/masked.nii" "$WORK/unmasked.nii" "$WORK/mask.nii"
# Signed input cannot be magnitude, so the AUTOMATIC VST stays off, exactly as
# an explicit -vst n would; that is what keeps a series rotated to real elsewhere
# (benchmark-large feeds one) from being stabilised as Rician.
"$GEN" mk mag-signed "$WORK/signed.nii" || { echo "fixture mag-signed failed"; exit 1; }
check_msg "signed input turns the VST off" "negative values so is not magnitude" \
	"$BIN" "$WORK/signed.nii" "$WORK/sg_auto.nii"
if ! "$BIN" -degibbs z 2>&1 | grep -q 'not compiled into this build'; then
	check_msg "signed input warns before -degibbs truncates it" "input has negative values, and -degibbs" \
		"$BIN" "$WORK/signed.nii" "$WORK/sg_dg.nii" -degibbs y
fi
check "signed input, -vst n" 0 "$BIN" "$WORK/signed.nii" "$WORK/sg_n.nii" -quiet -vst n
if cmp -s "$WORK/sg_auto.nii" "$WORK/sg_n.nii"
then ok "signed input: automatic VST == -vst n, byte for byte"; else bad "signed input was stabilised"; fi
"$BIN" "$M" "$WORK/x.nii" -phase "$WORK/phase-rad.nii" >"$WORK/stdout" 2>"$WORK/stderr"
if ! grep -q "vst  " "$WORK/stderr"
then ok "phase input turns the VST off"; else bad "phase input was stabilised as magnitude"; fi

echo "audit regressions"
# Each of these once failed: an output replaced the b-values it read; a noiseless
# patch (sigma 0) made optshrink write NaN; a one-slice image made the sphere
# report a false eigensolver failure.
echo "0 0 1000 1000 1000 2000 2000 3000" > "$WORK/values.nii"; cp "$WORK/values.nii" "$WORK/values.keep"
check_collision "output over -bval" "$BIN" "$M" "$WORK/values.nii" -quiet -demean y -bval "$WORK/values.nii"
if cmp -s "$WORK/values.nii" "$WORK/values.keep"
then ok "-bval survived the refusal"; else bad "-bval was overwritten"; fi
check "noiseless, optshrink, exclusive" 0 "$BIN" "$WORK/f48c.nii" "$WORK/nl1.nii" $LEGACY -quiet -filter optshrink
check "noiseless output is the input"   0 "$GEN" rl2 "$WORK/nl1.nii" "$WORK/f48c.nii" 1e-6
check "noiseless, defaults"             0 "$BIN" "$WORK/f48c.nii" "$WORK/nl2.nii" -quiet
check "noiseless default output is the input" 0 "$GEN" rl2 "$WORK/nl2.nii" "$WORK/f48c.nii" 1e-6
"$GEN" mk mag-thin "$WORK/thin.nii" || { echo "fixture mag-thin failed"; exit 1; }
check "one slice, sphere, averaging"    0 "$BIN" "$WORK/thin.nii" "$WORK/th1.nii" -quiet
"$GEN" mk rician-zero "$WORK/ricz.nii" || { echo "fixture rician-zero failed"; exit 1; }
check "isolated zero samples, defaults" 0 "$BIN" "$WORK/ricz.nii" "$WORK/rz.nii" -quiet -nthreads 1
# Skipping a lone zero as background (the old rule) moves this to ~172.9.
check "isolated zeros are stabilised with their voxel" 0 "$GEN" rms "$WORK/rz.nii" 158.896087 1e-6
# 100 volumes on 243 voxels: the main patch fits, and the larger sigma-estimation
# patch (k^3 >= 2N) once did not, failing a run the user's settings allowed.
"$GEN" mk mag-tall "$WORK/tall.nii" || { echo "fixture mag-tall failed"; exit 1; }
# Two volumes leave too few noise columns for a median sigma: it read ~6x too high
# and the VST zeroed most of the output.  The run must fall back and say so.
"$GEN" mk mag-pair "$WORK/pair.nii" || { echo "fixture mag-pair failed"; exit 1; }
check_msg "two volumes fall back to optthresh" "cannot give a noise level" "$BIN" "$WORK/pair.nii" "$WORK/pair_o.nii"
check "two volumes zero no voxel" 0 "$GEN" nozero "$WORK/pair_o.nii"
# ...and a noise map from them would be meaningless, under any filter.
check "two volumes refuse -noise" 1 "$BIN" "$WORK/pair.nii" "$WORK/x.nii" -quiet $LEGACY -noise "$WORK/pair_n.nii"
check "100 volumes on a 9x9x3 image" 0 "$BIN" "$WORK/tall.nii" "$WORK/tall_o.nii" -quiet
check "one slice, sphere, exclusive"    0 "$BIN" "$WORK/thin.nii" "$WORK/th2.nii" -quiet -aggregator exclusive -extent 5

echo "tensor mode"
# Recursive mode-wise denoising of space x time x echo tensors, against the
# matrix mode's own median-estimator noise map, with ME-PCA's window rules.
# Given a KNOWN variance the kernel and windows are checked against that
# package by tools/tensor_parity.py, which needs Python.  These pin what can
# be pinned without it: the noise map against the fixtures' true SD of 20, the
# window geometry, the refusals, determinism and numerical anchors.
for k in me-1 me-2 me-3 me-clean-1 me-clean-2 me-clean-3 me-shift-2 me-tr-2 me-short-2 me-neg-2 me-mask; do
	"$GEN" mk "$k" "$WORK/$k.nii" || { echo "fixture $k failed"; exit 1; }
done
E1=$WORK/me-1.nii; E2=$WORK/me-2.nii; E3=$WORK/me-3.nii
TEN="-tensor gd -tensor_window 3,3,3"
check "tensor: three echoes" 0 "$BIN" "$E1" "$WORK/t.nii" -quiet $TEN -echo "$E2" -echo "$E3" \
	-noise "$WORK/t_noise.nii" -rank "$WORK/t_rank.nii" -snr_gain "$WORK/t_gain.nii"
for k in 1 2 3; do
	before=$(err "$WORK/me-$k.nii" "$WORK/me-clean-$k.nii")
	after=$(err "$WORK/t_echo-$k.nii" "$WORK/me-clean-$k.nii")
	if awk -v a="$after" -v b="$before" 'BEGIN { exit !(a < 0.7 * b) }'
	then ok "tensor: echo $k closer to the truth ($before -> $after)"
	else bad "tensor: echo $k not denoised ($before -> $after)"; fi
done
# The estimated noise level against the oracle: the same run at the TRUE
# variance (SD 20) must not be much better than the map's.
check "tensor: at the true variance" 0 "$BIN" "$E1" "$WORK/tk.nii" -quiet $TEN -echo "$E2" -echo "$E3" \
	-tensor_sigma2 400
gd_err=$(err "$WORK/t_echo-1.nii" "$WORK/me-clean-1.nii")
or_err=$(err "$WORK/tk_echo-1.nii" "$WORK/me-clean-1.nii")
if awk -v g="$gd_err" -v o="$or_err" 'BEGIN { exit !(g < 1.05 * o) }'
then ok "tensor: estimated noise within 5% of the oracle's error ($gd_err against $or_err)"
else bad "tensor: estimated noise costs too much against the oracle ($gd_err against $or_err)"; fi
# Anchors: the arithmetic, the noise map and the ranks, so a change to any of
# them cannot cancel out.  The maps read 20.80 against the true 20: the median
# estimator runs ~4% high on 12-volume echoes, which costs nothing above.
check "tensor: echo 1 RMS anchor" 0 "$GEN" rms "$WORK/t_echo-1.nii" 667.614906 1e-6
check "tensor: noise map anchor" 0 "$GEN" rms "$WORK/t_noise.nii" 20.8027458 1e-6
check "tensor: rank map anchor"  0 "$GEN" rms "$WORK/t_rank.nii" 1.65457336 1e-6
check "tensor: every voxel covered" 0 "$GEN" nan "$WORK/t_noise.nii" 0
"$BIN" "$E1" "$WORK/t1.nii" -quiet $TEN -echo "$E2" -echo "$E3" -nthreads 1 -rank "$WORK/t1_rank.nii"
"$BIN" "$E1" "$WORK/t8.nii" $TEN -echo "$E2" -echo "$E3" -nthreads 8 -rank "$WORK/t8_rank.nii" 2>"$WORK/t8.err"
same=1
for f in echo-1 echo-2 echo-3 rank; do cmp -s "$WORK/t1_$f.nii" "$WORK/t8_$f.nii" || same=0; done
# Only meaningful if more than one worker really ran, so a host that has the
# cores must have used them; a serial-against-serial match proves nothing.
nth=$(sed -n 's/^.*threads *: *//p' "$WORK/t8.err")
ncpu=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)
if [ "$same" -ne 1 ]; then bad "tensor: output depends on the thread count"
elif [ "${nth:-0}" -ge 2 ]; then ok "tensor: 1 and $nth threads byte-identical"
elif [ "$ncpu" -ge 2 ]; then bad "tensor: only ${nth:-0} worker ran on $ncpu cores, so threads were not compared"
else ok "tensor: byte-identical (single-core host: no parallel run to compare)"; fi
check "tensor: sigma2 0 keeps all but the numerical-rank floor" 0 "$BIN" "$E1" "$WORK/t0.nii" -quiet $TEN -echo "$E2" -echo "$E3" -tensor_sigma2 0
check "tensor: sigma2 0 returns the input" 0 "$GEN" rl2 "$WORK/t0_echo-2.nii" "$E2" 1e-6
# Stride 4 with a 3-wide window: corners 0,4,8,12 in-plane cover 12 of 16, and
# only corner 0 fits in z, covering 3 of 5 slices.  1280 - 12*12*3 = 848 voxels
# are covered by nothing, and must read NaN in every map: once per echo in
# the noise map, once per mode in the rank map.  They keep their input value.
check "tensor: stride 4" 0 "$BIN" "$E1" "$WORK/ts.nii" -quiet $TEN -echo "$E2" -echo "$E3" -tensor_stride 4 \
	-noise "$WORK/ts_noise.nii" -rank "$WORK/ts_rank.nii"
check "tensor: stride 4 leaves 848 voxels uncovered" 0 "$GEN" nan "$WORK/ts_noise.nii" 2544
check "tensor: ... which keep their input value" 0 "$GEN" cmpnan "$WORK/ts_echo-1.nii" "$E1" "$WORK/ts_noise.nii" 1e-5
check "tensor: ... in every rank volume" 0 "$GEN" nan "$WORK/ts_rank.nii" 2544
# A window is denoised when it TOUCHES the mask: corners x 1-10, y 2-11, z 0-2
# reach it, covering 12x12x5 voxels and leaving 560, in each of 3 echoes.
check "tensor: mask" 0 "$BIN" "$E1" "$WORK/tm.nii" -quiet $TEN -echo "$E2" -echo "$E3" -mask "$WORK/me-mask.nii" \
	-noise "$WORK/tm_noise.nii"
check "tensor: windows touching the mask cover all but 560" 0 "$GEN" nan "$WORK/tm_noise.nii" 1680
check "tensor: in-plane window" 0 "$BIN" "$E1" "$WORK/tp.nii" -quiet -tensor gd -tensor_window 4,4 \
	-echo "$E2" -echo "$E3"
if cmp -s "$WORK/tp_echo-1.nii" "$WORK/t_echo-1.nii"
then bad "tensor: an in-plane window gave the local window's answer"; else ok "tensor: in-plane window differs from local"; fi
# Each option must reach the kernel: an output that matches the default means
# it was ignored, and an anchor pins what it does.  Demeaning moves the RMS by
# only 6e-6, so its anchor alone would barely separate it; the cmp does.
check "tensor: center assignment" 0 "$BIN" "$E1" "$WORK/tc.nii" -quiet $TEN -echo "$E2" -echo "$E3" -aggregator exclusive
check "tensor: demeaned" 0 "$BIN" "$E1" "$WORK/td.nii" -quiet $TEN -echo "$E2" -echo "$E3" -tensor_demean y
check "tensor: mode order" 0 "$BIN" "$E1" "$WORK/to.nii" -quiet $TEN -echo "$E2" -echo "$E3" -tensor_order time,space,echo
for opt in tc:-aggregator_exclusive td:-tensor_demean to:-tensor_order; do
	f=${opt%%:*}; name=${opt#*:}
	if cmp -s "$WORK/${f}_echo-1.nii" "$WORK/t_echo-1.nii"
	then bad "tensor: $name gave the default answer"; else ok "tensor: $name changes the output"; fi
done
check "tensor: center assignment RMS anchor" 0 "$GEN" rms "$WORK/tc_echo-1.nii" 668.046525 1e-6
check "tensor: demeaned RMS anchor" 0 "$GEN" rms "$WORK/td_echo-1.nii" 667.610777 1e-6
check "tensor: mode order RMS anchor" 0 "$GEN" rms "$WORK/to_echo-1.nii" 667.628282 1e-6
# -filter in tensor mode.  truncate is the default's cut without the shrinkage;
# optthresh cuts every mode HIGHER (lambda_* exceeds 1 + sqrt(beta)), so its
# space-mode rank can never exceed truncate's in any window, and must differ
# from both while still denoising.
check "tensor: -filter truncate" 0 "$BIN" "$E1" "$WORK/tf.nii" -quiet $TEN -echo "$E2" -echo "$E3" \
	-filter truncate -rank "$WORK/tf_rank.nii"
check "tensor: -filter optthresh" 0 "$BIN" "$E1" "$WORK/tt.nii" -quiet $TEN -echo "$E2" -echo "$E3" \
	-filter optthresh -rank "$WORK/tt_rank.nii"
for f in tf tt; do
	if cmp -s "$WORK/${f}_echo-1.nii" "$WORK/t_echo-1.nii"
	then bad "tensor: $f gave the default answer"; else ok "tensor: $f differs from optshrink"; fi
done
if cmp -s "$WORK/tt_echo-1.nii" "$WORK/tf_echo-1.nii"
then bad "tensor: optthresh gave truncate's answer"; else ok "tensor: optthresh differs from truncate"; fi
rf=$("$GEN" volrms "$WORK/tf_rank.nii" 0 1 0 || true); rt=$("$GEN" volrms "$WORK/tt_rank.nii" 0 1 0 || true)
if awk -v t="$rt" -v f="$rf" 'BEGIN { exit !(t > 0 && t < f) }'
then ok "tensor: optthresh keeps fewer space components ($rt against $rf)"
else bad "tensor: optthresh does not cut above the edge ($rt against $rf)"; fi
for k in 1 2 3; do
	before=$(err "$WORK/me-$k.nii" "$WORK/me-clean-$k.nii"); after=$(err "$WORK/tt_echo-$k.nii" "$WORK/me-clean-$k.nii")
	if awk -v a="$after" -v b="$before" 'BEGIN { exit !(a < 0.7 * b) }'
	then ok "tensor: optthresh echo $k closer to the truth ($before -> $after)"
	else bad "tensor: optthresh echo $k not denoised ($before -> $after)"; fi
done
check "tensor: truncate RMS anchor" 0 "$GEN" rms "$WORK/tf_echo-1.nii" 667.625707 1e-6
check "tensor: optthresh RMS anchor" 0 "$GEN" rms "$WORK/tt_echo-1.nii" 667.586047 1e-6
check "tensor: unknown -filter" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -filter optimal
check "-tensor_shrinker without -tensor" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet -tensor_shrinker gd
# -tensor_shrinker gd: Gavish & Donoho's observed dimensions.  It shrinks every
# kept value at least as hard as olesen's (which adds 4 r sigma^2), so it must
# differ from the default and still denoise.
check "tensor: -tensor_shrinker gd" 0 "$BIN" "$E1" "$WORK/tg.nii" -quiet $TEN -echo "$E2" -echo "$E3" \
	-tensor_shrinker gd
if cmp -s "$WORK/tg_echo-1.nii" "$WORK/t_echo-1.nii"
then bad "tensor: -tensor_shrinker gd gave olesen's answer"; else ok "tensor: -tensor_shrinker gd differs from olesen"; fi
for k in 1 2 3; do
	before=$(err "$WORK/me-$k.nii" "$WORK/me-clean-$k.nii"); after=$(err "$WORK/tg_echo-$k.nii" "$WORK/me-clean-$k.nii")
	if awk -v a="$after" -v b="$before" 'BEGIN { exit !(a < 0.7 * b) }'
	then ok "tensor: gd shrinker echo $k closer to the truth ($before -> $after)"
	else bad "tensor: gd shrinker echo $k not denoised ($before -> $after)"; fi
done
check "tensor: gd shrinker RMS anchor" 0 "$GEN" rms "$WORK/tg_echo-1.nii" 667.600445 1e-6
check "tensor: explicit olesen is the default" 0 "$BIN" "$E1" "$WORK/to2.nii" -quiet $TEN -echo "$E2" -echo "$E3" \
	-tensor_shrinker olesen
if cmp -s "$WORK/to2_echo-1.nii" "$WORK/t_echo-1.nii"
then ok "tensor: -tensor_shrinker olesen is byte-identical to the default"
else bad "tensor: -tensor_shrinker olesen differs from the default"; fi
check "tensor: a shrinker with nothing to shrink" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" \
	-filter truncate -tensor_shrinker gd
check "tensor: unknown shrinker" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -tensor_shrinker mp
# The flags tensor mode shares with matrix mode.  -aggregator gaussian is
# matrix mode's only; -preserve_noise_bias keeps the floor, so it must differ;
# -noise_dof describes the VST and is refused without it; and matrix mode
# takes one -phase, not one per echo.
check "tensor: -aggregator gaussian" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -aggregator gaussian
check "tensor: -aggregator uniform is the default" 0 "$BIN" "$E1" "$WORK/tu.nii" -quiet $TEN -echo "$E2" -echo "$E3" \
	-aggregator uniform
if cmp -s "$WORK/tu_echo-1.nii" "$WORK/t_echo-1.nii"
then ok "tensor: -aggregator uniform is byte-identical to the default"
else bad "tensor: -aggregator uniform differs from the default"; fi
check "tensor: -preserve_noise_bias" 0 "$BIN" "$E1" "$WORK/tb.nii" -quiet $TEN -echo "$E2" -echo "$E3" \
	-preserve_noise_bias
if cmp -s "$WORK/tb_echo-1.nii" "$WORK/t_echo-1.nii"
then bad "tensor: -preserve_noise_bias changed nothing"; else ok "tensor: -preserve_noise_bias keeps the floor"; fi
check "tensor: -noise_dof 4" 0 "$BIN" "$E1" "$WORK/tn.nii" -quiet $TEN -echo "$E2" -echo "$E3" -noise_dof 4
if cmp -s "$WORK/tn_echo-1.nii" "$WORK/t_echo-1.nii"
then bad "tensor: -noise_dof 4 changed nothing"; else ok "tensor: -noise_dof reaches the VST"; fi
check "tensor: -noise_dof without the VST" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -vst n -noise_dof 4
check "tensor: -aggregator_fwhm is matrix mode's" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" \
	-aggregator_fwhm 2
check "matrix: two -phase images" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet -phase "$E2" -phase "$E3"
check "tensor: echo on another grid"   1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$WORK/me-shift-2.nii"
check "tensor: echo with another TR"   1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$WORK/me-tr-2.nii"
check "tensor: echo with fewer volumes" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$WORK/me-short-2.nii"
check "tensor: negative values"        1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$WORK/me-neg-2.nii"
# Complex input: each echo rotated with its own phase.  The phase is smooth, so
# the rotation removes it and leaves the signal on the real axis under Gaussian
# noise, which must denoise towards the truth -- and differently from the same
# echoes taken as magnitude, or the rotation did nothing.
for k in 1 2 3; do
	for p in cmag cphase; do
		"$GEN" mk "me-$p-$k" "$WORK/me-$p-$k.nii" || { echo "fixture me-$p-$k failed"; exit 1; }
	done
done
C1=$WORK/me-cmag-1.nii; C2=$WORK/me-cmag-2.nii; C3=$WORK/me-cmag-3.nii
P1=$WORK/me-cphase-1.nii; P2=$WORK/me-cphase-2.nii; P3=$WORK/me-cphase-3.nii
check "tensor: complex" 0 "$BIN" "$C1" "$WORK/tz.nii" -quiet $TEN -echo "$C2" -echo "$C3" \
	-phase "$P1" -phase "$P2" -phase "$P3" -phaseunits radians
check "tensor: the same echoes as magnitude" 0 "$BIN" "$C1" "$WORK/tzm.nii" -quiet $TEN -echo "$C2" -echo "$C3"
if cmp -s "$WORK/tz_echo-3.nii" "$WORK/tzm_echo-3.nii"
then bad "tensor: -phase changed nothing"; else ok "tensor: -phase rotates the echoes"; fi
before=$(err "$C3" "$WORK/me-clean-3.nii"); after=$(err "$WORK/tz_echo-3.nii" "$WORK/me-clean-3.nii")
if awk -v a="$after" -v b="$before" 'BEGIN { exit !(a < 0.7 * b) }'
then ok "tensor: complex echo 3 closer to the truth ($before -> $after)"
else bad "tensor: complex echo 3 not denoised ($before -> $after)"; fi
# Echoes with their own noise levels (SD 40, 20, 10), as echo-specific
# reconstructions produce.  Each echo's noise volume must recover its own level,
# and taking each echo in its own units must beat one level pooled over all.
for k in 1 2 3; do "$GEN" mk "me-hetero-$k" "$WORK/me-hetero-$k.nii" || { echo "fixture me-hetero-$k failed"; exit 1; }; done
H1=$WORK/me-hetero-1.nii; H2=$WORK/me-hetero-2.nii; H3=$WORK/me-hetero-3.nii
check "tensor: echoes with their own noise" 0 "$BIN" "$H1" "$WORK/th.nii" -quiet $TEN -echo "$H2" -echo "$H3" \
	-noise "$WORK/th_noise.nii"
check "tensor: echo 1's noise volume reads its SD of 40" 0 "$GEN" volrms "$WORK/th_noise.nii" 0 40 0.1
check "tensor: echo 3's noise volume reads its SD of 10" 0 "$GEN" volrms "$WORK/th_noise.nii" 2 10 0.1
check "tensor: the same echoes, pooled" 0 "$BIN" "$H1" "$WORK/thp.nii" -quiet $TEN -echo "$H2" -echo "$H3" \
	-tensor_noise pooled
for k in 1 3; do
	own=$(err "$WORK/th_echo-$k.nii" "$WORK/me-clean-$k.nii"); pool=$(err "$WORK/thp_echo-$k.nii" "$WORK/me-clean-$k.nii")
	if awk -v a="$own" -v b="$pool" 'BEGIN { exit !(a < b) }'
	then ok "tensor: echo $k, own noise beats pooled ($own against $pool)"
	else bad "tensor: echo $k, pooling did as well as each echo's own noise ($pool against $own)"; fi
done
check "tensor: -vst n" 0 "$BIN" "$E1" "$WORK/tv.nii" -quiet $TEN -echo "$E2" -vst n
check "tensor: -vst y with -phase" 1 "$BIN" "$C1" "$WORK/x.nii" -quiet $TEN -echo "$C2" \
	-phase "$P1" -phase "$P2" -phaseunits radians -vst y
check "tensor: -vst y with -tensor_sigma2" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" \
	-tensor_sigma2 400 -vst y
check "tensor: -tensor_noise with -tensor_sigma2" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" \
	-tensor_sigma2 400 -tensor_noise pooled
check "tensor: unknown -tensor_noise" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -tensor_noise mp
# A noise-free echo gives no noise level.  Beside noisy echoes that is broken
# input, and once made the whole run "noiseless", returning every echo
# unchanged at exit 0; it must be refused.  When EVERY echo is noise-free there
# is nothing to remove, and the input comes back.
check "tensor: one noiseless echo among noisy ones" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN \
	-echo "$WORK/me-clean-2.nii" -echo "$E3"
if grep -q "echo 2 gives no noise level" "$WORK/stderr"
then ok "tensor: ... names the echo"; else bad "tensor: the noiseless echo was not named"; sed 's/^/       /' "$WORK/stderr"; fi
check "tensor: every echo noiseless" 0 "$BIN" "$WORK/me-clean-1.nii" "$WORK/tq.nii" -quiet $TEN \
	-echo "$WORK/me-clean-2.nii" -echo "$WORK/me-clean-3.nii"
check "tensor: ... returns the input" 0 "$GEN" rl2 "$WORK/tq_echo-2.nii" "$WORK/me-clean-2.nii" 1e-6
check "tensor: one phase for two echoes" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -phase "$P1"
check "tensor: -phaseunits needs -phase" 1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -phaseunits radians
check_collision "tensor: a phase that is an echo" \
	"$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -phase "$P1" -phase "$E2"
check "tensor: one echo"               1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN
check "-echo without -tensor"          1 "$BIN" "$E1" "$WORK/x.nii" -quiet -echo "$E2"
check "tensor: matrix option refused"  1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -shape cube
check "tensor: window too large"       1 "$BIN" "$E1" "$WORK/x.nii" -quiet -tensor gd -echo "$E2" -tensor_window 17,3,3
check "tensor: one-element window"     1 "$BIN" "$E1" "$WORK/x.nii" -quiet -tensor gd -echo "$E2" -tensor_window 3
check "tensor: bad mode order"         1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -tensor_order space,time,time
check "tensor: stride list length"     1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -tensor_stride 1,2
check "tensor: negative sigma2"        1 "$BIN" "$E1" "$WORK/x.nii" -quiet $TEN -echo "$E2" -tensor_sigma2 -1
# No Marchenko-Pastur fit is shipped, under any name.
check "tensor: mppca is not a method"  1 "$BIN" "$E1" "$WORK/x.nii" -quiet -tensor mppca -echo "$E2"
# The median estimator needs 3 noise columns; 2 volumes must be refused, not run.
"$GEN" mk mag-pair "$WORK/pair-e1.nii" && "$GEN" mk mag-pair "$WORK/pair-e2.nii"
check "tensor: 2 volumes cannot give a noise level" 1 "$BIN" "$WORK/pair-e1.nii" "$WORK/x.nii" -quiet \
	-tensor gd -echo "$WORK/pair-e2.nii" -tensor_window 3,3,3
check "tensor: ... but run with a known variance" 0 "$BIN" "$WORK/pair-e1.nii" "$WORK/tpair.nii" -quiet \
	-tensor gd -echo "$WORK/pair-e2.nii" -tensor_window 3,3,3 -tensor_sigma2 100
# A stride the parser accepts but no image can use once overflowed the tile
# arithmetic.  It must leave corner 0 alone: 1280 - 27 voxels uncovered, x 2 echoes.
check "tensor: stride INT32_MAX" 0 "$BIN" "$E1" "$WORK/tbig.nii" -quiet $TEN -echo "$E2" \
	-tensor_stride 2147483647 -noise "$WORK/tbig_noise.nii"
check "tensor: ... covers corner 0 only" 0 "$GEN" nan "$WORK/tbig_noise.nii" 2506
# Echo outputs are DERIVED names, so an input can sit exactly where one lands.
cp "$E2" "$WORK/c_echo-2.nii"; cp "$E2" "$WORK/c2.keep"
check_collision "tensor: echo output over an echo input" \
	"$BIN" "$E1" "$WORK/c.nii" -quiet $TEN -echo "$WORK/c_echo-2.nii"
if cmp -s "$WORK/c_echo-2.nii" "$WORK/c2.keep"
then ok "tensor: the echo input survived"; else bad "tensor: an echo input was overwritten"; fi
check_collision "tensor: -noise over an echo output" \
	"$BIN" "$E1" "$WORK/c.nii" -quiet $TEN -echo "$E2" -noise "$WORK/c_echo-1.nii"
# Echo 2 as a .hdr/.img pair, under the very name an extensionless output
# derives for it.  Outputs take echo 1's type, so echo 2's goes to m_echo-2.nii
# and the pair survives.  Each echo once kept its own type: the output was
# checked as m_echo-2.nii but written as the pair, over the input, at exit 0.
"$GEN" mk me-2@ "$WORK/m_echo-2" || { echo "fixture me-2@ failed"; exit 1; }
cp "$WORK/m_echo-2.img" "$WORK/m2img.keep"
check "tensor: echoes of mixed file types" 0 "$BIN" "$E1" "$WORK/m" -quiet $TEN -echo "$WORK/m_echo-2.hdr"
if cmp -s "$WORK/m_echo-2.img" "$WORK/m2img.keep" && [ -s "$WORK/m_echo-2.nii" ]
then ok "tensor: a pair-format echo survived, its output written as .nii"
else bad "tensor: a pair-format echo input was overwritten"; fi

echo
printf '%d passed, %d failed\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
