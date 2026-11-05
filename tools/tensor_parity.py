#!/usr/bin/env python3
"""Check svht_denoise's tensor mode against the ME-PCA Python package.

``-tensor gd`` estimates each echo's noise level from the matrix mode's
Gavish-Donoho median map, and puts each echo in units of it; ME-PCA estimates
one level with a Marchenko-Pastur fit, which svht_denoise does not implement.  So the two are compared at the same KNOWN noise variance,
where ME-PCA runs only its recursive pass and its window rules -- exactly what
svht_denoise shares with it -- and never its noise fit.  Three parts, all on
data generated here with a fixed seed:

1. Kernel: single tensors of assorted shapes go through
   ``mepca.tensor_pca.denoise_array_recursive_tensor`` (float64) and through
   ``dn_tensortest run``, which calls the C kernel directly.  Compared: every
   mode's rank and the denoised tensor.
2. Image: multi-echo series are written as one 4D NIfTI per echo and denoised by
   ``svht_denoise -tensor gd -tensor_sigma2`` and by
   ``mepca.tensor_pca.denoise_recursive_tensor`` on the same float32 values (as
   float64).  Compared: every denoised echo, the noise map, the per-mode rank
   maps and the SNR-gain map, including which voxels are uncovered.
3. Estimator: svht_denoise's own per-echo noise maps against the true SD, and
   its denoising error against the same run at the true variance (no reference).

SVD signs.  With ``subtract_mean`` the result depends on the sign of every
singular vector (later levels take means across the reshaped V*S), and LAPACK
defines no sign.  svht_denoise fixes one -- each V column's largest-magnitude
entry is positive -- and this script applies the same convention to the
reference, so the demeaned cases compare the algorithms rather than two LAPACK
accidents.  Without centring the sign cannot matter.  The demeaned kernel cases
also print how far the UNPATCHED reference moves.

Usage::

    python tools/tensor_parity.py --bin src/svht_denoise --kernel src/dn_tensortest

Needs numpy, nibabel and mepca (``pip install git+https://github.com/tsalo/ME-PCA@python-package``).
Exit status is 0 when every case is within tolerance.
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import mepca.tensor_pca as tensor_pca
import nibabel as nb
import numpy as np
from mepca.tensor_pca import denoise_array_recursive_tensor, denoise_recursive_tensor

LAPACK_SVD = tensor_pca._svd

# Tolerances.  The C kernel works in double through the smaller Gram of each
# unfolding while the reference calls LAPACK's SVD, so the two agree to roundoff;
# the image level is float32 on output and averages online rather than by
# sum-then-divide.
KERNEL_TOL = 1e-9
IMAGE_TOL = 1e-5
MODES = ("space", "time", "echo")
NOISE_SD = 20.0


def signed_svd(x):
    """Economy SVD with svht_denoise's sign convention.

    Each V column's largest-magnitude entry (the first, on a tie) is positive.
    """
    u, s, vh = LAPACK_SVD(x)
    biggest = vh[np.arange(vh.shape[0]), np.argmax(np.abs(vh), axis=1)]
    flip = biggest < 0
    u[:, flip] *= -1
    vh[flip] *= -1
    return u, s, vh


def tucker(shape, ranks, rng, sd):
    """Return a low Tucker rank tensor plus Gaussian noise, and the clean tensor."""
    core = rng.standard_normal(ranks) * 10
    factors = [rng.standard_normal((n, r)) for n, r in zip(shape, ranks)]
    clean = np.einsum("abc,ia,jb,kc->ijk", core, *factors) + 100.0
    return clean + sd * rng.standard_normal(shape), clean


def rel_l2(a, b):
    """Return the relative L2 difference of ``a`` from ``b``."""
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), np.finfo(float).tiny))


def run_kernel_case(kernel, x, shrink, demean, sigma2):
    """Denoise ``x`` with the C kernel; return the output and the ranks."""
    proc = subprocess.run(
        [
            str(kernel),
            "run",
            *(str(n) for n in x.shape),
            "optshrink" if shrink else "truncate",
            str(int(demean)),
            repr(float(sigma2)),
        ],
        input=x.astype(np.float64).tobytes(order="F"),
        capture_output=True,
        check=True,
    )
    out = np.frombuffer(proc.stdout, dtype=np.float64).reshape(x.shape, order="F")
    return out, np.array([int(f) for f in proc.stderr.decode().split()])


def reference_kernel(x, shrink, demean, sigma2):
    """Denoise ``x`` with the reference at a known variance."""
    return denoise_array_recursive_tensor(
        x,
        opt_shrink=shrink,
        subtract_mean=demean,
        sigma2=sigma2,
    )


def kernel_cases(kernel):
    """Compare the C kernel with the reference on single tensors."""
    rng = np.random.default_rng(20261104)
    shapes = [
        ((3, 7, 5), (2, 2, 2)),
        ((27, 40, 3), (3, 2, 2)),
        ((8, 60, 4), (2, 3, 2)),
        ((125, 80, 4), (4, 3, 2)),
        ((5, 125, 30), (2, 3, 2)),
        ((64, 30, 4), (3, 3, 2)),
    ]
    options = [
        dict(shrink=True, demean=False, sigma2=1.0),
        dict(shrink=False, demean=False, sigma2=1.0),
        dict(shrink=True, demean=False, sigma2=4.0),
        dict(shrink=True, demean=True, sigma2=1.0),
        dict(shrink=False, demean=True, sigma2=0.5),
        dict(shrink=True, demean=False, sigma2=0.0),
    ]
    failures = 0
    for shape, ranks in shapes:
        x, _ = tucker(shape, ranks, rng, 1.0)
        for opt in options:
            ref, _, ref_ranks = reference_kernel(x, **opt)
            out, c_ranks = run_kernel_case(kernel, x, **opt)
            d_out = rel_l2(out, ref)
            ok = d_out <= KERNEL_TOL and np.array_equal(c_ranks, ref_ranks)
            failures += not ok
            note = ""
            if opt["demean"]:
                tensor_pca._svd = LAPACK_SVD
                raw = reference_kernel(x, **opt)[0]
                tensor_pca._svd = signed_svd
                note = f" (LAPACK's own signs move the reference by {rel_l2(raw, ref):.2g})"

            print(
                f"  {'ok  ' if ok else 'FAIL'} kernel {shape} {opt}: rel L2 {d_out:.2g}, "
                f"ranks {c_ranks.tolist()} vs {ref_ranks.tolist()}{note}"
            )

    return failures


def multi_echo(shape, n_vols, tes, rng, sd):
    """Return synthetic multi-echo magnitude of shape (x, y, z, t, e), and its truth."""
    x, y, z = np.meshgrid(*(np.arange(n) for n in shape), indexing="ij")
    s0 = 800 + 300 * np.sin(0.4 * x) * np.cos(0.3 * y) + 40 * z
    r2s = 0.02 + 0.005 * np.cos(0.25 * (x + y))
    t = np.arange(n_vols)
    bold = (
        1
        + 0.05 * np.sin(0.6 * t[None, None, None, :] + 0.2 * x[..., None])
        + 0.03 * np.cos(0.9 * t)[None, None, None, :]
    )
    decay = np.exp(-np.asarray(tes)[None, None, None, :] * r2s[..., None])
    clean = s0[..., None, None] * bold[..., None] * decay[..., None, :]
    noisy = clean + sd * rng.standard_normal(clean.shape)
    return np.abs(noisy).astype(np.float32), clean


def write_echoes(data, directory, affine):
    """Write each echo of ``data`` as a 4D NIfTI; return the paths."""
    paths = []
    for i_echo in range(data.shape[-1]):
        path = directory / f"in_echo-{i_echo + 1}.nii"
        img = nb.Nifti1Image(data[..., i_echo], affine)
        img.header.set_xyzt_units("mm", "sec")
        img.header["pixdim"][4] = 2.0
        img.to_filename(path)
        paths.append(path)

    return paths


def run_cli(binary, case_dir, echoes, window, order, extra):
    """Run ``svht_denoise -tensor gd`` on ``echoes``, writing every map."""
    args = [
        str(binary),
        str(echoes[0]),
        str(case_dir / "out.nii"),
        "-quiet",
        "-tensor",
        "gd",
        "-tensor_window",
        ",".join(str(w) for w in window),
        "-tensor_order",
        ",".join(order),
        "-noise",
        str(case_dir / "noise.nii"),
        "-rank",
        str(case_dir / "rank.nii"),
        "-snr_gain",
        str(case_dir / "gain.nii"),
        *extra,
    ]
    for echo in echoes[1:]:
        args += ["-echo", str(echo)]

    subprocess.run(args, check=True)


def nonfinite_equal(a, b):
    """Return True where ``a`` and ``b`` are both finite or both not."""
    return np.array_equal(np.isfinite(a), np.isfinite(b))


def image_case(binary, name, data, window, options, mask, order, tmp):
    """Denoise ``data`` with the CLI and the reference; return 1 on a mismatch."""
    case_dir = Path(tempfile.mkdtemp(prefix=f"{name}_", dir=tmp))
    affine = np.diag([2.0, 2.0, 2.5, 1.0])
    echoes = write_echoes(data, case_dir, affine)
    options = {"sigma2": NOISE_SD**2, **options}
    c_flags = {
        "stride": ("-tensor_stride", lambda v: str(v)),
        "center_assign": ("-aggregator", lambda v: "exclusive" if v else "uniform"),
        "opt_shrink": ("-filter", lambda v: "optshrink" if v else "truncate"),
        "subtract_mean": ("-tensor_demean", lambda v: "y" if v else "n"),
        "sigma2": ("-tensor_sigma2", lambda v: repr(float(v))),
    }
    extra = []
    for key, value in options.items():
        flag, fmt = c_flags[key]
        extra += [flag, fmt(value)]

    if mask is not None:
        mask_path = case_dir / "mask.nii"
        nb.Nifti1Image(mask.astype(np.uint8), affine).to_filename(mask_path)
        extra += ["-mask", str(mask_path)]

    run_cli(binary, case_dir, echoes, window, order, extra)

    # The reference, on the same float32 values, in float64.
    groups = {"space": [0, 1, 2], "time": [3], "echo": [4]}
    ref_mask = None
    if mask is not None:
        ref_mask = mask.any(axis=2) if len(window) == 2 else mask

    denoised, sigma2_map, p_map, snr_gain = denoise_recursive_tensor(
        data.astype(np.float64),
        window,
        indices=[groups[m] for m in order],
        mask=ref_mask,
        **options,
    )
    if len(window) == 2:
        # The reference's maps for an in-plane window are per column; the C
        # tool writes them to every slice of that column.
        n_z = data.shape[2]
        sigma2_map = np.repeat(sigma2_map[:, :, None], n_z, axis=2)
        p_map = np.repeat(p_map[:, :, None, :], n_z, axis=2)
        snr_gain = np.repeat(snr_gain[:, :, None], n_z, axis=2)

    worst = 0.0
    for i_echo in range(data.shape[-1]):
        out = nb.load(case_dir / f"out_echo-{i_echo + 1}.nii").get_fdata()
        worst = max(worst, rel_l2(out, denoised[..., i_echo]))

    # -noise has one volume per echo; at a known variance every one is the same.
    noise_echoes = nb.load(case_dir / "noise.nii").get_fdata()
    noise = noise_echoes[..., 0]
    same_echoes = all(np.array_equal(noise_echoes[..., i], noise, equal_nan=True) for i in range(data.shape[-1]))
    rank = nb.load(case_dir / "rank.nii").get_fdata()
    gain = nb.load(case_dir / "gain.nii").get_fdata()
    ref_noise = np.sqrt(sigma2_map)
    covered = np.isfinite(ref_noise)
    d_noise = rel_l2(noise[covered], ref_noise[covered]) if covered.any() else 0.0
    d_rank = float(np.max(np.abs(rank[covered] - p_map[covered]))) if covered.any() else 0.0
    finite_gain = np.isfinite(snr_gain)
    d_gain = rel_l2(gain[finite_gain], snr_gain[finite_gain]) if finite_gain.any() else 0.0
    patterns = (
        same_echoes
        and nonfinite_equal(noise, ref_noise)
        and nonfinite_equal(rank, p_map)
        and nonfinite_equal(gain, snr_gain)
    )
    ok = (
        worst <= IMAGE_TOL
        and d_noise <= IMAGE_TOL
        and d_rank <= 1e-4
        and d_gain <= IMAGE_TOL
        and patterns
    )
    print(
        f"  {'ok  ' if ok else 'FAIL'} image {name}: echoes rel L2 {worst:.2g}, noise {d_noise:.2g}, "
        f"ranks max {d_rank:.2g}, gain {d_gain:.2g}, uncovered voxels {int((~covered).sum())}"
        f"{'' if patterns else ', NaN pattern differs'}"
    )
    return int(not ok)


def image_cases(binary, data, tmp):
    """Compare the CLI with the reference on whole multi-echo series."""
    mask = np.zeros(data.shape[:3], dtype=bool)
    mask[2:8, 3:9, 1:4] = True
    mask[4:6, 5:7, 2] = False
    cases = [
        ("local-3x3x3", (3, 3, 3), {}, None, MODES),
        ("in-plane-4x4", (4, 4), {}, None, MODES),
        ("larger-variance", (3, 3, 3), {"sigma2": 4 * NOISE_SD**2}, None, MODES),
        ("no-shrink", (3, 3, 3), {"opt_shrink": False}, None, MODES),
        ("demean", (3, 3, 3), {"subtract_mean": True}, None, MODES),
        ("stride-2-mask", (3, 3, 3), {"stride": 2}, mask, MODES),
        ("stride-4", (3, 3, 3), {"stride": 4}, None, MODES),
        ("center-mask", (3, 3, 3), {"center_assign": True}, mask, MODES),
        ("in-plane-center-mask", (4, 4), {"center_assign": True}, mask, MODES),
        ("order-time-space-echo", (3, 3, 3), {}, None, ("time", "space", "echo")),
        ("order-echo-time-space", (3, 3, 2), {}, None, ("echo", "time", "space")),
    ]
    failures = 0
    for name, window, options, case_mask, order in cases:
        failures += image_case(binary, name, data, window, options, case_mask, order, tmp)

    return failures


def estimator_case(binary, data, clean, tmp):
    """Report the Gavish-Donoho noise map against the truth; 1 if it costs > 5%."""
    case_dir = Path(tempfile.mkdtemp(prefix="estimator_", dir=tmp))
    echoes = write_echoes(data, case_dir, np.diag([2.0, 2.0, 2.5, 1.0]))
    errors = {}
    for label, extra in (("map", []), ("oracle", ["-tensor_sigma2", repr(NOISE_SD**2)])):
        run_cli(binary, case_dir, echoes, (3, 3, 3), MODES, extra)
        out = np.stack(
            [nb.load(case_dir / f"out_echo-{i + 1}.nii").get_fdata() for i in range(data.shape[-1])],
            axis=-1,
        )
        errors[label] = rel_l2(out, clean)
        if label == "map":
            noise = nb.load(case_dir / "noise.nii").get_fdata()

    ok = errors["map"] <= 1.05 * errors["oracle"]
    medians = " ".join(f"{np.median(noise[..., i]):.4g}" for i in range(noise.shape[-1]))
    print(
        f"  {'ok  ' if ok else 'FAIL'} estimator: per-echo noise map medians {medians} against a "
        f"true SD of {NOISE_SD:g}; error to truth {errors['map']:.4g}, at the true variance "
        f"{errors['oracle']:.4g}, input {rel_l2(data, clean):.4g}"
    )
    return int(not ok)


def main():
    """Run every case; exit non-zero on any mismatch."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--bin", type=Path, required=True, help="svht_denoise binary")
    parser.add_argument("--kernel", type=Path, required=True, help="dn_tensortest binary")
    args = parser.parse_args()

    tensor_pca._svd = signed_svd
    print("kernel parity, known variance")
    failures = kernel_cases(args.kernel)
    data, clean = multi_echo((12, 11, 6), 14, (15.0, 35.0, 55.0), np.random.default_rng(1104), NOISE_SD)
    with tempfile.TemporaryDirectory(prefix="tensor_parity_") as tmp:
        print("image parity, known variance")
        failures += image_cases(args.bin, data, tmp)
        print("estimated noise, against the truth")
        failures += estimator_case(args.bin, data, clean, tmp)

    print(f"\n{failures} failure{'' if failures == 1 else 's'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
