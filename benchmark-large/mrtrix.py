#!/usr/bin/env python3
"""Real-valued MPPCA denoising of the EDDEN sample data.

Reproduces the "MPPCA*" arm of Manzano Patron et al. (2024): the complex data
(magnitude + phase) are rotated onto the real axis, and MRtrix3 `dwidenoise` --
the *magnitude* MPPCA implementation -- is then run on that real-valued series.

The rotation is a NumPy port of NIFTI_COMP_to_REAL.m (Steen Moeller, 10/24/2023,
shipped with the EDDEN code), which is NIFTI_NORDIC.m with the denoising removed:

  1. scale the phase image to radians and form I = mag * exp(i*phase)
  2. remove the static (background) phase, taken from the first high-signal volume
  3. remove the residual dynamic phase, estimated per slice/volume with a
     low-pass k-space filter (Tukey^3 window)
  4. keep the real part

Usage:  python3 mrtrix.py MAG PHASE OUTDIR [--rotate-only]

All three paths are required.
"""

import subprocess
import sys
from pathlib import Path

import nibabel as nib
import numpy as np


PHASE_FILTER_WIDTH = 3  # ARG.phase_filter_width default


def tukey(n, alpha=1.0):
    """scipy.signal.windows.tukey / MATLAB tukeywin, symmetric."""
    x = np.linspace(0, 1, n)
    w = np.ones(n)
    lo, hi = x < alpha / 2, x >= 1 - alpha / 2
    w[lo] = 0.5 * (1 + np.cos(2 * np.pi / alpha * (x[lo] - alpha / 2)))
    w[hi] = 0.5 * (1 + np.cos(2 * np.pi / alpha * (x[hi] - 1 + alpha / 2)))
    return w


def lowpass_2d(vol):
    """Per-slice/volume low-pass filter in k-space, over the first two axes.

    One VOLUME at a time.  The filter is independent per (slice, volume), so
    this is the same arithmetic as filtering the 4D block in one call -- but
    numpy's FFT promotes to complex128, and at 140x140x92x297 that intermediate
    is 8 GiB with several live at once.  Streaming bounds it at ~60 MB.
    """
    nx, ny = vol.shape[:2]
    wx = (tukey(nx) ** PHASE_FILTER_WIDTH)[:, None, None]
    wy = (tukey(ny) ** PHASE_FILTER_WIDTH)[None, :, None]
    out = np.empty(vol.shape, dtype=np.complex64)
    for t in range(vol.shape[3]):
        v = vol[..., t]
        k = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(v, axes=(0, 1)), axes=(0, 1)), axes=(0, 1))
        k *= wx
        k *= wy
        out[..., t] = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(k, axes=(0, 1)), axes=(0, 1)), axes=(0, 1))
    return out


def comp_to_real(mag_file, phase_file, out_file):
    mag_img = nib.load(mag_file)
    mag = np.abs(np.asarray(mag_img.dataobj, dtype=np.float32))
    phs = np.asarray(nib.load(phase_file).dataobj, dtype=np.float32)

    # Scale the raw phase (e.g. Siemens -4096..4094) to radians, centred on 0.
    span = phs.max() - phs.min()
    centre = (phs.max() + phs.min()) / span / 2
    phs = (phs / span - centre) * 2 * np.pi
    print(f"Phase range {phs.min():.2f} to {phs.max():.2f} rad (should be -pi to pi)")

    img = (mag * np.exp(1j * phs)).astype(np.complex64)

    # Normalise by the smallest non-zero magnitude of the first volume.
    v0 = np.abs(img[..., 0])
    scale = v0[v0 != 0].min()
    img /= scale

    # Static phase: reference is the first volume whose mean magnitude is
    # within 5% of the maximum (i.e. a b=0 volume).
    tt = np.abs(img).mean(axis=(0, 1, 2))
    ref = int(np.argmax(tt > 0.95 * tt.max()))
    img *= np.exp(-1j * np.angle(img[..., ref]))[..., None]

    # Residual dynamic phase.
    img *= np.exp(-1j * np.angle(lowpass_2d(img)))

    real = np.nan_to_num(img.real) * scale
    nib.save(nib.Nifti1Image(real.astype(np.float32), mag_img.affine, mag_img.header), out_file)
    return out_file


def main():
    # --rotate-only stops after writing real.nii.gz.  benchmark_large needs that:
    # if this script runs dwidenoise itself, the benchmark cannot time it, and
    # the run it did time would silently be the rotation alone.
    args = [a for a in sys.argv[1:] if a != "--rotate-only"]
    rotate_only = len(args) != len(sys.argv) - 1
    if len(args) != 3:
        sys.exit("usage: mrtrix.py MAG PHASE OUTDIR [--rotate-only]")
    mag, phase, out = Path(args[0]), Path(args[1]), Path(args[2])
    out.mkdir(parents=True, exist_ok=True)
    real = comp_to_real(mag, phase, out / "real.nii.gz")
    if rotate_only:
        return

    # Default MRtrix3 magnitude MPPCA: patch extent is chosen automatically as
    # the smallest odd P with P^3 >= number of volumes.
    subprocess.run(
        ["dwidenoise", str(real), str(out / "denoised.nii.gz"),
         "-noise", str(out / "noise.nii.gz"), "-force"],
        check=True,
    )


if __name__ == "__main__":
    main()
