#!/bin/sh
# Repeatable timing for the large-N tuning target.
#
#   scripts/bench.sh <input> <output> <binary> [binary...]
#
# Prints the MEDIAN of three runs per binary -- wall clock, CPU time and peak RSS
# -- and refuses to print anything on a busy machine.  That refusal is the point:
# several measurements during this work were spent on a loaded host and had to be
# thrown away, and a wrong number costs more than no number because it gets
# written down.
#
# With more than one binary the rounds ALTERNATE, and reverse on every other
# round, so the two never sit in the same thermal position twice.  Running one
# binary three times and then the other three times gives the second one a hotter
# machine, which is worth a few percent -- the same order as the changes being
# measured.
#
# Median rather than minimum: the minimum flatters whichever run happened to
# catch a cold machine, and three samples is not enough to trust a tail.
#
# Environment:
#   ARGS=...      extra arguments passed to every binary (default -quiet)
#   MAXLOAD=<n>   refuse above this 1-minute load average (default 2.0)
#   RUNS=<n>      rounds to run (default 3; odd, so the median is a real sample)

set -u

[ $# -ge 3 ] || { echo "usage: $0 <input> <output> <binary> [binary...]" >&2; exit 2; }
IN=$1; OUT=$2; shift 2

ARGS=${ARGS:--quiet}
MAXLOAD=${MAXLOAD:-2.0}
RUNS=${RUNS:-3}

[ -r "$IN" ] || { echo "cannot read: $IN" >&2; exit 2; }
for b in "$@"; do [ -x "$b" ] || { echo "not executable: $b" >&2; exit 2; }; done
[ -x /usr/bin/time ] || { echo "/usr/bin/time not found" >&2; exit 2; }

# The output name is positional and therefore easy to get wrong, and everything
# downstream WRITES to it.  Omit it -- `bench.sh in.nii ./old ./new` -- and OUT
# silently becomes the first binary, which the tool then overwrites with a NIfTI
# while the run reports a clean set of timings: the baseline you were measuring
# against is gone and nothing says so.  Reproduced, which is why this is here.
#
# Compared after resolving each path, not as the strings given: `./old` and `old`
# name one file and only one of the two spellings would ever match.  An earlier
# version tested `-x "$OUT"` instead, which is a poor proxy twice over -- it
# refuses any output on an exFAT or SMB mount, where every file is mode 0755, and
# it would not notice a non-executable file being clobbered.
resolve() { (cd "$(dirname -- "$1")" 2>/dev/null && printf '%s/%s\n' "$(pwd -P)" "$(basename -- "$1")"); }
sv_out=$(resolve "$OUT")
[ "$sv_out" = "$(resolve "$IN")" ] && { echo "output is the input: $OUT" >&2; exit 2; }
for b in "$@"; do
	[ "$sv_out" = "$(resolve "$b")" ] && { echo "output is one of the binaries: $OUT" >&2; exit 2; }
done
# A symlinked output is written THROUGH, so it destroys whatever it points at --
# and none of the checks above see that target.  This used to be covered by
# accident: the `rm -f "$OUT"` that once ran before each timed run unlinked the
# symlink and left the target alone.  Removing that `rm` fixed one hazard and
# opened this one, so it is named explicitly now rather than relied upon.
[ -L "$OUT" ] && { echo "output is a symlink, refusing: $OUT" >&2; exit 2; }

# Power state, before load, because it is the larger error and the quieter one.
# On Apple silicon, battery power and Low Power Mode each move work onto the
# efficiency cores: measured here, the identical pair of binaries took 95 s and
# 189 s of wall clock on mains and on battery.  It inflates CPU TIME as well as
# wall clock, so neither column gives it away, and it hits both candidates
# equally -- a RATIO survives it untouched (1.681x against 1.675x on the same
# pair) while every absolute number is out by two.  Nothing in the numbers says
# which of the two a timing was, which is exactly why this refuses rather than
# warns.
if [ "$(uname -s)" = Darwin ] && [ "${ALLOW_LOWPOWER:-0}" != 1 ]; then
	if pmset -g batt 2>/dev/null | head -1 | grep -q "Battery Power"; then
		echo "running on battery; refusing to time (set ALLOW_LOWPOWER=1 to override)" >&2
		exit 1
	fi
	# `powermode` on Apple silicon, `lowpowermode` on Intel.  Match either, and
	# only refuse on a value actually read: an unrecognised pmset would otherwise
	# yield the empty string, which is `!= 0`, and refuse for ever.
	lpm=$(pmset -g 2>/dev/null | awk '$1 == "powermode" || $1 == "lowpowermode" { print $2 }')
	# 0 normal, 1 Low Power, 2 High Power.  Refuse ONLY 1: `!= 0` also blocked
	# High Power Mode, which is the mode you would deliberately benchmark in.
	if [ "$lpm" = 1 ]; then
		echo "Low Power Mode is on; refusing to time (set ALLOW_LOWPOWER=1 to override)" >&2
		exit 1
	fi
fi

# 1-minute load average.  uptime's wording differs between platforms ("load
# averages:" against "load average:") but the three numbers are always the last
# three fields, so take them positionally.
#
# Refuse an unparseable one rather than continuing: awk compares an empty `l`
# against `m` as STRINGS, which is false, so a reshaped or missing uptime would
# sail through the gate below and print a banner with a hole in it.  A gate that
# fails open is worse than no gate, because the run still looks measured.
load=$(uptime | awk '{ n = split($0, f, " "); gsub(",", "", f[n-2]); print f[n-2] }')
# Both sides of the comparison, not just the one read from the system: awk
# compares two non-numeric strings LEXICALLY, so `MAXLOAD=two` disabled the gate
# entirely and still printed a banner.  Validating one and not the other is how
# the hole survived the round that added this check.
case $load in
''|*[!0-9.]*) echo "cannot read the load average from uptime; refusing to time" >&2; exit 1 ;;
esac
case $MAXLOAD in
''|*[!0-9.]*) echo "MAXLOAD must be numeric (got '$MAXLOAD')" >&2; exit 2 ;;
esac
if awk -v l="$load" -v m="$MAXLOAD" 'BEGIN { exit !(l > m) }'; then
	echo "machine is loaded (1-min load $load > $MAXLOAD); refusing to time" >&2
	echo "  set MAXLOAD to override, but then do not compare against a quiet run" >&2
	exit 1
fi

case $(uname -s) in
Darwin) TIMEFLAG=-l ;;   # BSD time: peak RSS in BYTES
*)      TIMEFLAG=-v ;;   # GNU time: peak RSS in KIBIBYTES
esac

tmp=$(mktemp -d "${TMPDIR:-/tmp}/svht_bench.XXXXXX") || exit 1
trap 'rm -rf "$tmp"' EXIT INT TERM

# Say which build each binary is, for the same reason `make` says which zlib and
# which reduction went in: two binaries can differ by a compile-time switch with
# nothing in the timings to show it, and this project has already lost one
# comparison that way.
echo "bench: $IN  [$ARGS]  load $load  ${RUNS} rounds"
for b in "$@"; do
	sv_kern=portable
	otool -L "$b" 2>/dev/null | grep -q 'Accelerate\.framework' && sv_kern=Accelerate
	otool -L "$b" 2>/dev/null | awk '/^\t/ { print $1 }' | grep -q 'libz\.' &&
		sv_zl="system zlib" || sv_zl="static zlib-ng"
	printf '  %-28s %s, %s\n' "$(basename "$b")" "$sv_kern" "$sv_zl"
done

# One timed run; appends "wall cpu rss" to that binary's sample file.
run_one() {
	# No `rm -f "$OUT"` here.  svht_denoise overwrites an existing output anyway,
	# so the removal bought nothing and cost a great deal when OUT was not what
	# the caller meant -- see the guards above.
	# shellcheck disable=SC2086 -- ARGS is deliberately word-split
	/usr/bin/time $TIMEFLAG "$1" "$IN" "$OUT" $ARGS >/dev/null 2>"$tmp/raw" || {
		echo "run failed: $1" >&2; sed 's/^/    /' "$tmp/raw" >&2; exit 1; }
	# Both time(1) formats in one pass.  BSD puts "real user sys" on a single line
	# and the peak RSS on its own; GNU labels each on its own line.
	#
	# An unset awk variable prints as 0, so a format this does not recognise --
	# another time(1), a locale, or the tool's own stderr landing in the same file
	# -- would otherwise contribute a perfectly plausible `0.00 0.00 0.0` to the
	# median.  The GNU branch in particular has never run on the machine these
	# numbers came from.  Exit non-zero instead and let the caller see it.
	awk -v flag="$TIMEFLAG" '
		flag == "-l" {
			if ($2 == "real") { wall = $1; cpu = $3 + $5 }
			if (/maximum resident set size/) rss = $1 / 1048576
		}
		flag == "-v" {
			if (/Elapsed .wall clock/) { split($NF, p, ":")
			                             wall = (length(p) == 3) ? p[1]*3600 + p[2]*60 + p[3] : p[1]*60 + p[2] }
			if (/User time/)                 user = $NF
			if (/System time/)               sys  = $NF
			if (/Maximum resident set size/)  rss = $NF / 1024
		}
		END { if (flag == "-v") { if (user == "" || sys == "") exit 1; cpu = user + sys }
		      if (wall == "" || cpu == "" || rss == "") exit 1
		      printf "%.2f %.2f %.1f\n", wall, cpu, rss }
	' "$tmp/raw" >>"$2" || {
		echo "could not parse $(basename /usr/bin/time) $TIMEFLAG output:" >&2
		sed 's/^/    /' "$tmp/raw" >&2; exit 1; }
	printf '    %-28s %s\n' "$(basename "$1")" "$(tail -1 "$2")"
}

r=0
while [ "$r" -lt "$RUNS" ]; do
	r=$((r + 1))
	echo "  round $r"
	# Reverse the order on even rounds so no binary is always last in the round.
	if [ $((r % 2)) -eq 1 ]; then
		i=0
		for b in "$@"; do i=$((i + 1)); run_one "$b" "$tmp/s$i"; done
	else
		# `\${$i}`, not `\$$i`: with ten or more binaries the latter expands as
		# $1 followed by a literal 0 under a POSIX sh.
		i=$#
		while [ "$i" -ge 1 ]; do
			eval "b=\${$i}"
			run_one "$b" "$tmp/s$i"
			i=$((i - 1))
		done
	fi
done

echo "  medians (wall s / cpu s / peak RSS MiB)"
med() { cut -d' ' -f"$2" "$1" | sort -n | awk 'NR == int((n+1)/2)' n="$RUNS"; }
i=0
for b in "$@"; do
	i=$((i + 1))
	printf '    %-28s wall %-8s cpu %-9s rss %s\n' \
		"$(basename "$b")" "$(med "$tmp/s$i" 1)" "$(med "$tmp/s$i" 2)" "$(med "$tmp/s$i" 3)"
done
