## Introduction

svht_denoise provides a principled, data-driven approach for noise reduction in multi-volume NIfTI datasets, such as diffusion-weighted MRI. The tool applies Singular Value Decomposition (SVD) over local image patches to identify and truncate noise components using the optimal Singular Value Hard Thresholding (SVHT) framework derived by Gavish and Donoho (2014).

The tool works on magnitude volumes, and on complex data when the phase is also available. Magnitude reconstruction shifts thermal Gaussian noise into a non-Gaussian (Rician / non-central chi) distribution whose mean is not zero, and that noise floor biases exactly the low-SNR voxels denoising is meant to help. Given the phase, `-phase` rotates the complex data onto the real axis first, so the noise stays zero-mean Gaussian and the floor is avoided rather than denoised.

<img src=images/anim_compare.gif width="300">

## Building

Requires a C99 compiler on a POSIX system (macOS, Linux, WSL), plus libm, pthreads and zlib — all standard on those platforms. On macOS the build supplies its own zlib (see below), so only libm and pthreads come from the system.

```
git clone --recursive https://github.com/rordenlab/svht_denoise
cd svht_denoise/src
make
make test     # optional: regression tests, no dependencies
```

`--recursive` fetches [zlib-ng](https://github.com/zlib-ng/zlib-ng), the one bundled dependency, pinned to an exact commit. On macOS `make` builds it and links it statically, which is roughly **3× faster at writing `.nii.gz`** — deflate is the largest serial block in a compressed run — and leaves the binary with no `libz` dependency at all. Reading is unchanged; inflate was never the bottleneck.

Nothing breaks without it. A plain `git clone` leaves the submodule empty, and a build path containing a space — or any character outside letters, digits, dot, underscore and hyphen — defeats zlib-ng's own build; either way the build falls back to the system zlib and prints which one it used, and why:

```
make ZLIBNG=0                 # force the system zlib
make ZLIBNG_ROOT=<dir>        # use a prebuilt ZLIB_COMPAT zlib-ng tree
git submodule update --init ../third_party/zlib-ng    # from src, after a plain clone
```

Elsewhere the system zlib is the default, since this was measured only on Apple Silicon; `ZLIBNG_ROOT` still works there if you want it.

This writes the `svht_denoise` executable into `src`. Run `./svht_denoise -help` for the full list of options.

```
./svht_denoise dwi.nii.gz denoised.nii.gz
./svht_denoise dwi.nii.gz denoised.nii.gz -mask brain.nii.gz -noise sigma.nii.gz
```

## Usage

The `-help` option provides full usage details. In general, this tool behaves similarly to MRtrix dwidenoise and therefore has similar options and considerations. For example with regards to [patch size](https://mrtrix.readthedocs.io/en/latest/dwi_preprocessing/denoising.html#patch-size).

## Complex data

Pass the magnitude series as the input and its phase with `-phase`, and the two are rotated onto the real axis before denoising:

```
./svht_denoise mag.nii.gz denoised.nii.gz -phase phase.nii.gz
```

The phase must match the magnitude in all four dimensions and share its world transform, which is checked. `-real` writes the rotated series before denoising, if you want to see what the denoiser was given.

Only the phase *scale* has to be right; any constant offset is removed along with the background phase. By default the units are taken from the observed range, and the convention used is reported on stderr (unless `-quiet`) so it can be checked:

| observed | read as |
| --- | --- |
| a span within 0.1 of 2π | radians, unchanged, wherever centred |
| any other span | whatever unit makes that span one turn |
| constant | a constant rotation, removed whatever its scale |

The middle row is the rule the reference implementation uses, and it is correct for **any** encoding that covers a whole turn — scanner integers (Siemens writes -4096..4094), degrees, cycles, [-1, 1] — without needing to know which of them it is. Every acquired phase image covers a turn, because its background is noise whose phase is uniform over the circle.

It is wrong for data that does *not* cover a turn, and no amount of cleverness fixes that: a range of [-0.5, 0.5] is either half a radian or one whole cycle, and nothing in the file says which. For that case, say so — `-phaseunits radians|degrees|turns` (`cycles` is a synonym for `turns`; `auto` is the default) skips the inference entirely:

```
./svht_denoise mag.nii.gz denoised.nii.gz -phase phase.nii.gz -phaseunits radians
```

The rotation removes the per-voxel static background phase, taken from the first volume whose mean magnitude reaches 95% of the largest (a b=0 image in a diffusion series), then removes the residual per-slice, per-volume phase left by bulk motion, estimated by low-pass filtering the complex image in plane. The real part is kept. This is the "MPPCA\*" preprocessing of Manzano Patron et al. (2024), and the rotation itself is the phase-removal front end of NORDIC (Moeller et al., 2021), specified by `NIFTI_COMP_to_REAL.m`. No phase unwrapping is involved anywhere: every rotation is computed as `conj(z)/|z|` straight from the complex value, so wrapped angles never arise.

**The output is then signed.** Background voxels scatter about zero instead of piling up on the noise floor, which is the entire point, but it means the result is no longer a magnitude image. On the sample data in [edden/](edden/) the rotation takes the background mean from 70.0 to 43.8, makes 22% of background voxels negative, and raises the background standard deviation from 78.4 to 81.7. Those figures are the *rotated series*, what `-real` writes; after denoising, the same background reads 43.4, 3.2% negative, standard deviation 71.1. The rise in spread is expected rather than a regression — magnitude reconstruction rectifies noise onto the positive half-line, and undoing that returns the signed channel's own scatter.

That residual +43.8 is not a bug and is not removable by this method: the phase estimate is derived from the data, so each noise sample is partly rotated onto the alignment of its own neighbourhood. The reference MRtrix3 pipeline leaves the same bias — 43.6 against 43.8 on the rotated series, and 43.37 against 43.41 after denoising.

The rotation is serial, so it adds about 0.37 s of wall clock to a 2.1 s run (~18%) while being only ~2.7% of total CPU time — the gap between those two figures is why a single wall-clock sample badly misjudges it. It costs rather more memory: a float32 copy of the phase series plus a fixed set of 3D scratch volumes, taking peak RSS on the sample data from 127 to 211 MiB.

## macOS release packages

A signed, notarized, stapled installer can be built from `src`. It installs an Apple Silicon (arm64) binary into `/usr/local/bin`, targets macOS 11 or newer, and links only against libSystem — zlib-ng is embedded statically, so there is no `libz` dependency. The packaging script refuses to build without the zlib-ng submodule, and separately proves the static link happened by checking that the finished binary does not name `libz` at all; the system-zlib fallback that is right for a working copy would otherwise ship silently, being ~3× slower with nothing about the binary to say so.

Releases are arm64 only. This is a compute-bound tool, and the Intel Macs a universal build would have served are the slowest machines that could run it — on a platform macOS 26 is the last release to support. Building from source on an Intel Mac is unaffected: a plain `make` produces a native binary there as it always did. Cutting a *release* does require an Apple Silicon host, because the packaging script runs the test suite against the arm64 binary it builds; it says so and stops before building if run elsewhere. The installer itself declares `hostArchitectures="arm64"`, so an Intel Mac refuses it rather than installing a binary that cannot run.

Store the notarization credentials once per machine — this prompts for an [app-specific password](https://appleid.apple.com), not your Apple ID password:

```
make macos-notary-profile APPLE_ID='you@example.com' TEAM_ID='ABCDE12345'
```

Then build a release. Two *different* certificates are needed, both from the same developer account: "Developer ID Application" signs the executable, "Developer ID Installer" signs the `.pkg`.

```
make macos-release \
  VERSION=0.1.20260808 \
  MACOS_SIGN_IDENTITY='Developer ID Application: Your Name (ABCDE12345)' \
  MACOS_INSTALLER_IDENTITY='Developer ID Installer: Your Name (ABCDE12345)'
```

That checks the credentials actually authenticate before doing any work, then builds, runs the full test suite, signs, packages, notarizes, staples and verifies — writing to `dist/`. The binary that ships is the one the tests ran against: the build and the suite happen in a single `make` invocation, so a stale or differently-configured binary cannot be substituted between them. With a single arm64 slice, the suite exercises exactly the code that ships — the universal build could only ever test the host slice and signed the other one unexecuted. `security find-identity -v` lists your identities.

For local testing without certificates, `make macos-pkg-adhoc` produces an unsigned package; Gatekeeper will refuse it on another machine. `scripts/verify_macos_pkg.sh <pkg>` checks a package without installing it — signature, that nothing outside `/usr/lib` and `/System` is linked, and that the binary runs; it reports the architectures rather than asserting them, since the arm64-only gate runs in `package_macos.sh` before anything is signed. That last check *runs* the packaged executable, so use it on packages you built, not on one from an untrusted source.

## Licensing

The underlying algorithm for optimal singular value hard thresholding (SVHT) is based on the work of Gavish and Donoho (2014).

 - Algorithm & Theory: The method computes the unknown-noise-level threshold of Equation 4, tau = omega(beta) * y_med, where omega(beta) = lambda_star(beta) / sqrt(mu_beta) is assembled from the closed form of Equation 11 and an exactly bisected Marchenko-Pastur median, rather than from the cubic approximation to omega given in Equation 5. While reference MATLAB scripts for the method were distributed under GPL-3.0, this project is an independent clean-room C implementation built directly from the paper's mathematical formulations rather than derived from the original reference code.
 - CLI Design: The command-line interface, flag names, and operational parameters were intentionally structured to follow the conventions of dwidenoise from MRtrix3 wherever the two tools share a concept (`-mask`, `-noise`, `-extent`, `-nthreads`, `-quiet`). No code from MRtrix3 was incorporated.
 - Complex Data: The `-phase` front end implements the real-axis rotation that Manzano Patron et al. (2024) [PMID: 40800437](https://pubmed.ncbi.nlm.nih.gov/40800437/) apply before MPPCA, and which originates as the phase-removal stage of NORDIC (Moeller et al., NeuroImage 226:117539, 2021). Steen Moeller's MATLAB `NIFTI_COMP_to_REAL.m` was read as a specification of the method; no code was translated from it, and the two deliberate divergences, along with the derivation that replaces its per-slice FFTs with a 7-tap convolution, are recorded at the head of `dn_phase.c`. The convention of reading the phase units off the observed data range follows [niimath](https://github.com/rordenlab/niimath) (BSD-2-Clause), whose `medic.c` does the same.
 - AI Assistance: Generative AI tools were used during translation and refactoring.
 - Prior Work: Applying random matrix theory to denoise diffusion MRI was established by Veraart, Novikov, Fieremans and colleagues (2016), whose [dwidenoise](https://mrtrix.readthedocs.io/en/dev/reference/commands/dwidenoise.html) remains the reference implementation and the benchmark this tool is measured against. svht_denoise is an independent implementation of a different estimator, and shares no code with it: dwidenoise fits the Marchenko-Pastur distribution to the patch eigenspectrum to estimate the noise level sigma, which then sets the cutoff, whereas svht_denoise applies the Gavish and Donoho hard threshold directly, at omega(beta) times the median singular value, so the noise level never enters the retain-or-discard decision. An estimate is still available on request: the optional `-noise` map reports sigma = y_med / sqrt(M * mu_beta), Equation 26 of the paper, but it is derived after the fact and does not influence the denoised output.

This executable and its source code are licensed under the Mozilla Public License 2.0 (MPL-2.0), matching the permissive license used by the core MRtrix3 codebase. Note that dwidenoise itself is not covered by that permissive license: it carries an additional notice from New York University and the University of Antwerp restricting it to non-commercial research and excluding clinical care. svht_denoise carries no such restriction.

The macOS release binary statically embeds [zlib-ng](https://github.com/zlib-ng/zlib-ng), © 1995-2024 Jean-loup Gailly and Mark Adler, under the zlib licence — permissive and compatible with MPL-2.0 §3.3. Its text is in [packaging/LICENSE-zlib-ng.txt](packaging/LICENSE-zlib-ng.txt). A `make ZLIBNG=0` build links the system zlib instead and embeds nothing.

## Benchmarking

A sample diffusion dataset is provided to evaluate the time and peak memory required to process a provided diffusion scan. This demonstrates that the performance on these metrics is similar to [dwidenoise](https://mrtrix.readthedocs.io/en/dev/reference/commands/dwidenoise.html). However, note these are different algorithms and the outputs are not expected to be identical. Further, the input image is intentionally small to minimize the size of this repository: 100x100x54 voxels at 2.2mm isotropic resolution, masked, with 36 volumes (a b=0 image, 30 half-sphere b=2000 directions, and five trailing b=0 images).

Run these from `src`, after `mkdir -p ../benchmark/svht_denoise ../benchmark/dwidenoise`:

```
/usr/bin/time -l ./svht_denoise ../benchmark/input/dwi.nii.gz ../benchmark/svht_denoise/dwi.nii.gz
/usr/bin/time -l ./svht_denoise ../benchmark/input/dwi.nii.gz ../benchmark/svht_denoise/dwi1.nii.gz -nthreads 1
/usr/bin/time -l dwidenoise ../benchmark/input/dwi.nii.gz ../benchmark/dwidenoise/dwi.nii.gz
/usr/bin/time -l dwidenoise ../benchmark/input/dwi.nii.gz ../benchmark/dwidenoise/dwi1.nii.gz -nthreads 1
```

And, from the repository root, with a Python environment that has DIPY (`mkdir -p benchmark/dipy` first):

```
/usr/bin/time -l python3 benchmark/dipy_mppca.py benchmark/input/dwi.nii.gz benchmark/dipy/dwi.nii.gz
```

Here are findings for an Apple MacBook with M4 Pro CPU (14 core, 10 performance). Time is elapsed wall clock, and peak RAM is the maximum resident set size, both as reported by `/usr/bin/time -l`.

The `svht_denoise` rows are the default macOS build, with the bundled static zlib-ng, and are the median of five runs (spread 1650-1680 ms and 14070-14260 ms respectively). Since which zlib went in changes the wall clock and nothing else, a timing is only comparable to another taken with the same one — `make` prints which it used. The `dwidenoise` rows are unchanged and were measured in an earlier session.

Note how differently the two rows moved: zlib-ng took 430 ms off the 14-thread run (2080 → 1650, a 21% saving) but only ~650 ms off the single-threaded one (14750 → 14100, about 4%). Compression is serial either way, so it is a large share of the wall clock only once the denoise itself is spread across cores — which is also why it was worth fixing.

The `dipy mppca` row is DIPY 1.11.0 on Python 3.12 (NumPy 2.1.3), median of three runs, at `patch_radius=2` — a 5x5x5 window, the same patch the other two use here. It has no thread option and uses one core (user time equals wall clock), so there is no 14-thread row for it; the figure is whole-process, of which interpreter startup and imports are 0.3 s. It is the closest algorithmic comparison to `dwidenoise` rather than to this tool: both are Marchenko-Pastur PCA, and svht_denoise is a different estimator. The memory is the more interesting column — DIPY works in float64 on the whole series, where the two C tools stream float32.

| Method       | Threads | Time (ms) | Peak RAM (MiB) |
| ------------ | ------- | --------- | ------------- |
| svht_denoise |      14 |      1650 |           152 |
| svht_denoise |       1 |     14100 |           151 |
| dwidenoise   |      14 |      1950 |           228 |
| dwidenoise   |       1 |     11620 |           226 |
| dipy mppca   |       1 |     59590 |           849 |

**That table predates a round of internal optimisation and has not been re-measured.** Measured before-and-after in one sitting on this same bundled series, at 14 threads and median of five alternating rounds, that work went from 1.79 s to 1.50 s of wall clock (-16%) and from 19.69 s to 16.01 s of CPU time (-19%), with peak RSS unchanged at 152 MiB. Those absolute times are not comparable with the table's 1650 ms — different sitting, different conditions — which is why the rows are left as they were measured rather than partially updated. Re-measuring all three tools in one sitting is still outstanding.

The saving grows with the series. On a 100x100x58 series of 138 volumes at 7x7x7 (not in this repository, same 14 threads, mains power, quiet machine, median of five alternating rounds) the same change reads 94.66 s to 58.67 s of wall clock (-38%) and 1220.96 s to 726.46 s of CPU (-40%), peak RSS 778.1 to 775.5 MiB. Nothing about the interface or the output format changed. One of those changes does move values — a local `hypot` in the eigensolver's QL inner loop — at 2.24e-10 relative L2 error on that series, with the rank map byte-identical.

`scripts/bench.sh <input> <output> <binary>...` produced those figures. It alternates candidates and reverses their order every other round so neither sits in the same thermal position twice, reports the median of the rounds, and refuses to time at all on a loaded machine, on battery, or in Low Power Mode — battery alone was measured at 2x, and nothing in the numbers afterwards says which state a timing was taken in.

## Links

 - [Gavish and Donoho](https://ieeexplore.ieee.org/document/6846297) describe the singular value hard thresholding method. The pre-print [arxiv.1305.5870v3](https://arxiv.org/abs/1305.5870) is included with this repository.
 - [dwidenoise](https://mrtrix.readthedocs.io/en/dev/reference/commands/dwidenoise.html) is a popular and robust open source tool for de-noising diffusion images. It is based on Veraart et al. (2016) works [PMID:27523449](https://pubmed.ncbi.nlm.nih.gov/27523449/) and [PMID:26599599](https://pubmed.ncbi.nlm.nih.gov/26599599/). Note that it is restricted to non-commercial research usage.
 - [dwidenoise2](https://github.com/Lestropie/dwidenoise2) uses a novel method for de-noising, albeit it is also restricted to non-commercial usage.
