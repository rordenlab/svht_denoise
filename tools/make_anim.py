#!/usr/bin/env python3
"""Build an animated grayscale GIF comparing denoising methods.

One column, three rows (top to bottom): input data, MRtrix3 dwidenoise,
svht_denoise. Each frame is one diffusion volume; frames are rendered with
niimath and assembled with Pillow.

Paths default to the repository root, so the script runs from any directory.

Needs nibabel, numpy and Pillow, plus a niimath binary on the PATH.

Example:
    python3 tools/make_anim.py                       # y 0.3, volumes 1..30
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

ROWS = [
    ("input", "input"),
    ("dwidenoise", "MRtrix3 dwidenoise"),
    ("svht_denoise", "svht_denoise"),
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-b", "--benchmark", type=Path, default=ROOT / "benchmark",
                   help="folder holding input/, dwidenoise/, svht_denoise/ (default: <repo>/benchmark)")
    p.add_argument("-n", "--name", default="dwi.nii.gz",
                   help="NIfTI filename inside each method folder (default: dwi.nii.gz)")
    p.add_argument("-f", "--first", type=int, default=1,
                   help="first volume to display, indexed from 0 (default: 1)")
    p.add_argument("-l", "--last", type=int, default=30,
                   help="last volume to display, inclusive (default: 30)")
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


def window(path, first, last, pct):
    """Robust intensity window taken from the displayed volumes of the input."""
    data = np.asanyarray(nib.load(str(path)).dataobj)[..., first:last + 1]
    return 0.0, float(np.percentile(data, pct))


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
    if args.last < args.first:
        sys.exit(f"--last ({args.last}) precedes --first ({args.first})")

    volumes = [args.benchmark / folder / args.name for folder, _ in ROWS]
    for nii in volumes:
        if not nii.is_file():
            sys.exit(f"missing volume: {nii}")

    nvol = nib.load(str(volumes[0])).shape[3]
    if args.last >= nvol:
        sys.exit(f"--last {args.last} exceeds the final volume index {nvol - 1}")

    lo, hi = args.range if args.range else window(volumes[0], args.first, args.last, args.pct)
    print(f"intensity window {lo:g}..{hi:g}, volumes {args.first}..{args.last}")

    frames = []
    with tempfile.TemporaryDirectory() as tmp:
        png = Path(tmp) / "frame.png"
        for vol in range(args.first, args.last + 1):
            panels = [render(niimath, nii, vol, args, lo, hi, png) for nii in volumes]
            if args.labels:
                panels = [label(p, text) for p, (_, text) in zip(panels, ROWS)]
            width = max(p.width for p in panels)
            tile = Image.new("L", (width, sum(p.height for p in panels)))
            y = 0
            for p in panels:
                tile.paste(p, ((width - p.width) // 2, y))
                y += p.height
            frames.append(tile)
            print(f"\rvolume {vol}", end="", flush=True)
    print()

    frames[0].save(args.out, save_all=True, append_images=frames[1:],
                   duration=args.duration, loop=0, optimize=True)
    print(f"wrote {args.out} ({frames[0].width}x{frames[0].height}, {len(frames)} frames)")


if __name__ == "__main__":
    main()
