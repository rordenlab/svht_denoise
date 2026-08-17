#!/usr/bin/env python3
"""Build an animated grayscale GIF comparing denoising methods.

One column, three rows (top to bottom): input data, MRtrix3 dwidenoise,
svht_denoise. Each frame is one diffusion volume; frames are rendered with
niimath and assembled with Pillow.

Paths default to the repository root, so the script runs from any directory.

Needs nibabel, numpy and Pillow, plus a niimath binary on the PATH.

Volumes are chosen by DIFFUSION WEIGHTING, not by index: every volume whose
b-value reaches --minb is shown.  The b-values come from a .bval beside the
input, found automatically.  Without one, the --first/--last range is used as-is.

Compressed inputs are inflated ONCE into a temporary directory and every frame
is rendered from that. gzip is not seekable, so niimath must otherwise inflate a
whole series to reach one volume, once per volume per row -- which is what left
earlier runs on the 2 GB large series unfinished. Budget room for one
uncompressed copy of each row (--tmpdir moves it).

Example:
    python3 tools/make_anim.py                       # y 0.3, b >= 50
    python3 tools/make_anim.py --minb 1500           # the shells that matter
    python3 tools/make_anim.py -o z -s 0.5 -f 1 -l 30

    # the large series, every b > 1500 volume, rows with different filenames:
    python3 tools/make_anim.py -b benchmark-large --minb 1500 -s 0.45 -d 150 \\
        --bval benchmark-large/dwi/sub-01_ses-1p5mm_dir-AP_dwi.bval \\
        -R "mrtrix/real.nii.gz=input (complex, rotated to real)" \\
        -R "mrtrix/denoised.nii.gz=MRtrix3 MPPCA*" \\
        -R "svht_denoise/denoised.nii.gz=svht_denoise" \\
        /tmp/anim_EDDEN1p5mm.gif

    NOT images/anim_EDDEN1p5mm.gif: the copy README embeds lives there and is
    hand-compressed to ~2.4 MB, so writing to it replaces it with a ~7 MB
    regenerated one.  Compress first, then move it into place.
"""

import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import nibabel as nib
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent.parent

DEFAULT_ROWS = ["input=input", "dwidenoise=MRtrix3 dwidenoise", "svht_denoise=svht_denoise"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-b", "--benchmark", type=Path, default=ROOT / "benchmark",
                   help="folder holding input/, dwidenoise/, svht_denoise/ (default: <repo>/benchmark)")
    p.add_argument("-n", "--name", default="dwi.nii.gz",
                   help="NIfTI filename inside each method folder (default: dwi.nii.gz)")
    p.add_argument("-f", "--first", type=int, default=0,
                   help="first volume to consider, indexed from 0 (default: 0)")
    p.add_argument("-l", "--last", type=int, default=None,
                   help="last volume to consider, inclusive (default: the final volume)")
    p.add_argument("--minb", type=float, default=50.0,
                   help="show only volumes whose b-value reaches this (default: 50, "
                        "which drops the b=0 images)")
    p.add_argument("--bval", type=Path, default=None,
                   help="FSL .bval file (default: the only *.bval beside the input)")
    p.add_argument("-R", "--rows", action="append", metavar="FOLDER[/FILE]=LABEL",
                   help="one row per flag, top to bottom; give FOLDER/FILE when the rows "
                        "do not share a filename (default: %s)" % ", ".join(DEFAULT_ROWS))
    p.add_argument("-o", "--orient", choices=["x", "y", "z"], default="y",
                   help="slice orientation (default: y)")
    p.add_argument("-s", "--slice", type=float, default=0.3,
                   help="slice position, 0..1 fraction or negative for absolute (default: 0.3)")
    p.add_argument("-r", "--range", type=float, nargs=2, metavar=("MIN", "MAX"),
                   help="intensity window shared by all rows (default: robust range of input)")
    p.add_argument("-p", "--pct", type=float, default=99.5,
                   help="percentile of the input defining the window maximum (default: 99.5)")
    p.add_argument("-z", "--zoom", type=float, default=1.0,
                   help="niimath scale factor for each panel (default: 1)")
    p.add_argument("-d", "--duration", type=int, default=250,
                   help="milliseconds per frame (default: 250)")
    p.add_argument("--labels", action=argparse.BooleanOptionalAction, default=True,
                   help="draw a row label on each panel (default: on)")
    p.add_argument("--niimath", default="niimath", help="path to the niimath binary")
    p.add_argument("--tmpdir", type=Path, default=None,
                   help="where to inflate .nii.gz inputs (default: the system temp dir). "
                        "Needs room for one uncompressed copy of every row.")
    p.add_argument("out", nargs="?", type=Path, default=ROOT / "images" / "anim_compare.gif",
                   help="output GIF (default: <repo>/images/anim_compare.gif). Do NOT aim this "
                        "at images/anim_EDDEN1p5mm.gif, the copy README embeds: it is "
                        "hand-compressed and would be overwritten.")
    return p.parse_args()


def window(path, vols, pct):
    """Robust intensity window taken from the displayed volumes of the input.

    Indexes the proxy per volume rather than materialising the series: on a
    decompressed cache that reads only the volumes shown, where
    asanyarray(dataobj) pulled all 297 of them (2.1 GB) to look at 180.
    """
    proxy = nib.load(str(path)).dataobj
    data = np.stack([np.asanyarray(proxy[..., v]) for v in vols])
    return 0.0, float(np.percentile(data, pct))


def uncompressed_bytes(path):
    """What `path` will occupy once inflated, from its header alone."""
    img = nib.load(str(path))
    n = int(np.prod(img.shape)) * int(img.get_data_dtype().itemsize)
    return n + 352


# THE REASON THIS SCRIPT IS USABLE ON A LARGE SERIES.  gzip is not seekable, so
# niimath must inflate the WHOLE file to reach one volume, and the frame loop
# calls it once per volume per row: 540 full inflations of a 2 GB series for a
# 180-frame three-row animation, which is where earlier runs on the
# large series were left to die.  Inflating ONCE up front turns that into three.
#
# It is worth this much and no more.  Measured, niimath does NOT read an
# uncompressed file whole -- it seeks to the volume it was asked for, so three
# renders off a 173 MB .nii took 0.119 s TOTAL.  So there is nothing left for the
# concurrency an older comment here proposed; the remaining cost is the inflation
# itself, which is serial in gzip and already down to once per row.
def uncompress(path, tmpdir, tag, quiet=False):
    """Inflate `path` into `tmpdir` once; return the path to render from.

    `tag` makes the temporary name unique per ROW, which is not decoration: the
    default rows are input/dwi.nii.gz, dwidenoise/dwi.nii.gz and
    svht_denoise/dwi.nii.gz, so keying the cache on the basename alone had all
    three overwrite one file and every row render from whichever landed last.
    """
    if path.suffix != ".gz":
        return path
    out = Path(tmpdir) / f"row{tag}.nii"
    t0 = time.monotonic()
    with gzip.open(path, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst, length=8 << 20)
    if not quiet:
        mb = out.stat().st_size / (1 << 20)
        print(f"  inflated {path.name} -> {mb:.0f} MiB in {time.monotonic() - t0:.1f} s")
    return out


def find_bval(folder, given):
    """The .bval to read b-values from, or None."""
    if given is not None:
        return given if given.is_file() else sys.exit(f"missing --bval file: {given}")
    hits = sorted(folder.glob("*.bval"))
    if len(hits) == 1:
        return hits[0]
    # Zero is the ordinary case for data that never had one; more than one is
    # ambiguous and guessing would silently animate the wrong shells.
    if len(hits) > 1:
        print(f"several .bval files in {folder}; pass --bval to choose")
    return None


def select(volumes_path, folder, args):
    """Which volume indices to animate, and why."""
    nvol = nib.load(str(volumes_path)).shape[3]
    last = nvol - 1 if args.last is None else args.last
    if last >= nvol:
        sys.exit(f"--last {last} exceeds the final volume index {nvol - 1}")
    if last < args.first:
        sys.exit(f"--last ({last}) precedes --first ({args.first})")
    span = list(range(args.first, last + 1))

    bval = find_bval(folder, args.bval)
    if bval is None:
        print(f"no .bval found, showing volumes {args.first}..{last} unfiltered")
        return span
    b = np.loadtxt(bval, ndmin=1)
    if b.size != nvol:
        sys.exit(f"{bval} has {b.size} b-values for {nvol} volumes")
    keep = [v for v in span if b[v] >= args.minb]
    if not keep:
        sys.exit(f"no volume in {args.first}..{last} reaches b >= {args.minb:g} "
                 f"(the data runs {b.min():g}..{b.max():g})")
    print(f"b >= {args.minb:g} keeps {len(keep)} of {len(span)} volumes "
          f"(b {b[keep].min():g}..{b[keep].max():g})")
    return keep


def render(niimath, nii, vol, args, lo, hi, png):
    cmd = [niimath, str(nii), "-crop", str(vol), "1", "-bitmap", "-u", "0"]
    if args.zoom != 1.0:
        cmd += ["-s", f"{args.zoom:g}"]
    cmd += ["-t", f"{lo:g}", f"{hi:g}", f"-{args.orient}", f"{args.slice:g}", str(png)]
    png.unlink(missing_ok=True)
    # niimath exits non-zero after a successful -bitmap -s, so judge it by its output
    done = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if not png.is_file():
        sys.exit(f"niimath failed ({' '.join(cmd)}):\n{done.stderr}")
    img = Image.open(png).convert("L")
    img.load()
    return img


def label(img, text):
    """Write text in the top-left corner, with a dark halo for legibility."""
    img = img.copy()
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    for dx, dy in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        draw.text((3 + dx, 2 + dy), text, fill=0, font=font)
    draw.text((3, 2), text, fill=255, font=font)
    return img


def main():
    args = parse_args()
    niimath = shutil.which(args.niimath)
    if niimath is None:
        sys.exit(f"niimath not found: {args.niimath}")

    # A row is FOLDER=LABEL, and FOLDER may carry a filename when the rows do not
    # share one -- comparing a rotated input against two denoised series is the
    # ordinary case for complex data, and those three are not called the same.
    rows, volumes = [], []
    for spec in (args.rows or DEFAULT_ROWS):
        where, _, text = spec.partition("=")
        rows.append((where, text or where))
        volumes.append(args.benchmark / where if "/" in where
                       else args.benchmark / where / args.name)
    for nii in volumes:
        if not nii.is_file():
            sys.exit(f"missing volume: {nii}")

    vols = select(volumes[0], volumes[0].parent, args)

    frames = []
    with tempfile.TemporaryDirectory(dir=args.tmpdir) as tmp:
        # Refuse up front rather than 90% through a long run: three inflated
        # copies of the large series are 6.4 GB, and the failure when they do not
        # fit lands after the slowest step.
        need = sum(uncompressed_bytes(n) for n in volumes if n.suffix == ".gz")
        if need:
            free = shutil.disk_usage(tmp).free
            print(f"inflating {sum(1 for n in volumes if n.suffix == '.gz')} series "
                  f"({need / (1 << 30):.1f} GiB) into {tmp}")
            if need > free * 0.95:
                sys.exit(f"not enough room in {tmp}: need {need / (1 << 30):.1f} GiB, "
                         f"{free / (1 << 30):.1f} GiB free. Point --tmpdir somewhere roomier.")
            volumes = [uncompress(n, tmp, i) for i, n in enumerate(volumes)]

        lo, hi = args.range if args.range else window(volumes[0], vols, args.pct)
        print(f"intensity window {lo:g}..{hi:g}")

        png = Path(tmp) / "frame.png"
        for n, vol in enumerate(vols, 1):
            panels = [render(niimath, nii, vol, args, lo, hi, png) for nii in volumes]
            if args.labels:
                panels = [label(p, text) for p, (_, text) in zip(panels, rows)]
            width = max(p.width for p in panels)
            tile = Image.new("L", (width, sum(p.height for p in panels)))
            y = 0
            for p in panels:
                tile.paste(p, ((width - p.width) // 2, y))
                y += p.height
            frames.append(tile)
            print(f"\rvolume {vol} ({n}/{len(vols)})", end="", flush=True)
    print()

    frames[0].save(args.out, save_all=True, append_images=frames[1:],
                   duration=args.duration, loop=0, optimize=True)
    print(f"wrote {args.out} ({frames[0].width}x{frames[0].height}, {len(frames)} frames)")


if __name__ == "__main__":
    main()
