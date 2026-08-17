#!/usr/bin/env python3
"""MP-PCA denoising with DIPY, for comparison against the root README benchmark.

    python3 benchmark/dipy_mppca.py <in.nii.gz> <out.nii.gz>

patch_radius=2 is a 5x5x5 window, matching what svht_denoise selects
automatically for a 36-volume series (smallest odd k with k^3 > N) and what
dwidenoise uses by default -- the three are compared on the same patch.

Timed as a whole process, so the figure includes interpreter startup and
imports. That is what a user waits for; the in-process split is printed on
stderr so the two can be told apart.
"""
import sys
import time

t_start = time.perf_counter()

import nibabel as nib
from dipy.denoise.localpca import mppca

t_import = time.perf_counter()

if len(sys.argv) != 3:
    sys.exit(f"usage: {sys.argv[0]} <in.nii.gz> <out.nii.gz>")
src, dst = sys.argv[1], sys.argv[2]

img = nib.load(src)
data = img.get_fdata()
t_load = time.perf_counter()

denoised = mppca(data, patch_radius=2)
t_denoise = time.perf_counter()

nib.save(nib.Nifti1Image(denoised, img.affine, img.header), dst)
t_end = time.perf_counter()

print(
    f"import {t_import - t_start:6.2f} s   "
    f"load {t_load - t_import:6.2f} s   "
    f"mppca {t_denoise - t_load:7.2f} s   "
    f"save {t_end - t_denoise:6.2f} s   "
    f"total {t_end - t_start:7.2f} s",
    file=sys.stderr,
)
