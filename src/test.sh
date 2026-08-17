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
check "48-volume run" 0 "$BIN" "$WORK/big48.nii" "$WORK/b48.nii" -quiet -rank "$WORK/b48_r.nii"
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
check "48-volume, 1 thread"      0 "$BIN" "$WORK/big48.nii" "$WORK/b48_t1.nii" -nthreads 1 -quiet
# NOT -quiet, and the worker count is read back, for the same reason the
# determinism section below does it: dn_effective_threads caps by DN_CHUNK-sized
# chunks and by the core count, so on a one-core host "-nthreads 8" means 1 and
# the comparison below would be a serial run against itself, reporting success.
check "48-volume, many threads"  0 "$BIN" "$WORK/big48.nii" "$WORK/b48_t8.nii" -nthreads 8
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
check "run, 1 thread"            0 "$BIN" "$M" "$WORK/t1.nii" -nthreads 1 -phase "$WORK/phase-rad.nii"
check "run, many threads"        0 "$BIN" "$M" "$WORK/t8.nii" -nthreads 8 -phase "$WORK/phase-rad.nii"
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

echo
printf '%d passed, %d failed\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
