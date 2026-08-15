#!/usr/bin/env python3
"""Diff two NIfTI-1 images (.nii or .nii.gz), numpy-only, no nibabel.

usage: python3 nii_diff.py <a.nii[.gz]> <b.nii[.gz]> [--thresh F]
"""
import argparse
import sys
import gzip
import struct
import numpy as np

DT_NAMES = {2: "uint8", 4: "int16", 8: "int32", 16: "float32", 64: "float64", 512: "uint16"}
DT_FMT = {2: "u1", 4: "i2", 8: "i4", 16: "f4", 64: "f8", 512: "u2"}


def read_nii(path):
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    if len(raw) < 348:
        raise ValueError(f"{path}: too small to be a NIfTI-1 header")

    # sizeof_hdr tells us the byte order: 348 native, else byteswapped.
    if struct.unpack("<i", raw[0:4])[0] == 348:
        e = "<"
    elif struct.unpack(">i", raw[0:4])[0] == 348:
        e = ">"
    else:
        raise ValueError(f"{path}: sizeof_hdr is not 348 in either byte order")

    dim = struct.unpack(e + "8h", raw[40:56])
    ndim = dim[0]
    if not (1 <= ndim <= 7):
        raise ValueError(f"{path}: bad dim[0]={ndim}")
    shape = tuple(int(s) for s in dim[1 : 1 + ndim])

    datatype = struct.unpack(e + "h", raw[70:72])[0]
    vox_offset = struct.unpack(e + "f", raw[108:112])[0]
    scl_slope = struct.unpack(e + "f", raw[112:116])[0]
    scl_inter = struct.unpack(e + "f", raw[116:120])[0]
    if datatype not in DT_FMT:
        raise ValueError(f"{path}: unsupported datatype code {datatype}")

    np_dtype = np.dtype(e + DT_FMT[datatype])
    count = int(np.prod(shape))
    offset = int(round(vox_offset))
    if offset + count * np_dtype.itemsize > len(raw):
        raise ValueError(f"{path}: file truncated relative to header dims")

    data = np.frombuffer(raw, dtype=np_dtype, count=count, offset=offset).astype(np.float64)
    if scl_slope != 0:
        data = data * scl_slope + scl_inter
    return data, shape, DT_NAMES[datatype]


def main():
    p = argparse.ArgumentParser(description="Diff two NIfTI-1 images.")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--thresh", type=float, default=1e-4,
                   help="report voxels differing by more than this times RMS(a) (default: 1e-4)")
    opts = p.parse_args()
    args, thresh = [opts.a, opts.b], opts.thresh

    try:
        a, shape_a, dtype_a = read_nii(args[0])
        b, shape_b, dtype_b = read_nii(args[1])
    except (OSError, ValueError) as ex:
        print(f"error: {ex}", file=sys.stderr)
        sys.exit(2)

    print(f"a: {args[0]}  shape={shape_a}  dtype={dtype_a}")
    print(f"b: {args[1]}  shape={shape_b}  dtype={dtype_b}")
    if shape_a != shape_b:
        print(f"error: shape mismatch {shape_a} vs {shape_b}", file=sys.stderr)
        sys.exit(2)

    finite_a = np.isfinite(a)
    rms_a = float(np.sqrt(np.mean(a[finite_a] ** 2))) if finite_a.any() else float("nan")
    max_abs_a = float(np.max(np.abs(a[finite_a]))) if finite_a.any() else float("nan")
    print(f"reference: RMS(a)={rms_a:.6g}  max|a|={max_abs_a:.6g}")

    diff = a - b
    abs_diff = np.abs(diff)
    if np.isfinite(abs_diff).any():
        idx = int(np.nanargmax(abs_diff))
        max_diff = float(abs_diff[idx])
        coord = tuple(int(c) for c in np.unravel_index(idx, shape_a, order="F"))
        print(f"max abs diff: {max_diff:.6g}  at flat index {idx}  voxel {coord}")
    else:
        print("max abs diff: n/a (no finite differences)")

    norm_a = float(np.sqrt(np.nansum(a**2)))
    norm_diff = float(np.sqrt(np.nansum(diff**2)))
    rel_l2 = norm_diff / norm_a if norm_a != 0 else float("nan")
    print(f"relative L2: {rel_l2:.6g}")

    mask = abs_diff > thresh * rms_a
    n_diff = int(np.sum(mask))
    pct = 100.0 * n_diff / a.size
    print(f"voxels with |a-b| > {thresh:g}*RMS(a): {n_diff} ({pct:.4g}%)")

    if n_diff > 0:
        coords = np.unravel_index(np.flatnonzero(mask), shape_a, order="F")
        xyz = np.stack(coords[:3], axis=-1)
        n_xyz = len(np.unique(xyz, axis=0))
        if len(shape_a) > 3:
            n_vol = len(np.unique(np.ravel_multi_index(coords[3:], shape_a[3:])))
        else:
            n_vol = 1
        print(f"  clustering: {n_diff} voxels, {n_xyz} distinct (x,y,z) columns, {n_vol} distinct volume indices")


if __name__ == "__main__":
    main()
