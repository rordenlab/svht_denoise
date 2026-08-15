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

Example:
    python3 tools/make_anim.py                       # y 0.3, b >= 50
    python3 tools/make_anim.py --minb 1500           # the shells that matter
    python3 tools/make_anim.py -o z -s 0.5 -f 1 -l 30
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
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
    p.add_argument("out", nargs="?", type=Path, default=ROOT / "images" / "anim_compare.gif",
                   help="output GIF (default: <repo>/images/anim_compare.gif, the one README.md embeds)")
    return p.parse_args()


def window(path, vols, pct):
    """Robust intensity window taken from the displayed volumes of the input."""
    data = np.asanyarray(nib.load(str(path)).dataobj)[..., vols]
    return 0.0, float(np.percentile(data, pct))


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


# SLOW ON A .nii.gz, AND BADLY: gzip is not seekable, so niimath inflates the
# WHOLE series to reach one volume, and this calls it once per volume per row --
# 540 full decompressions of a 2 GB series for a 180-frame three-row animation.
#
# Two fixes, in order of what they are worth. Decompressing ONCE to a temporary
# .nii and rendering every frame from that turns O(frames) inflations into one,
# which is the ~500x. Running the niimath calls concurrently is the easy ~14x on
# top. Neither is done here because the script has only ever been used on the
# small bundled series, where the whole thing takes seconds.
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

    lo, hi = args.range if args.range else window(volumes[0], vols, args.pct)
    print(f"intensity window {lo:g}..{hi:g}")

    frames = []
    with tempfile.TemporaryDirectory() as tmp:
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
