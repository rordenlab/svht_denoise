## Introduction

svht_denoise removes thermal noise from multi-volume NIfTI datasets such as diffusion-weighted MRI. It decomposes local image patches by SVD and truncates the noise components using the optimal Singular Value Hard Thresholding of Gavish and Donoho (2014), so no noise level has to be supplied.

It works on magnitude volumes, and on complex data when the phase is also available. Magnitude reconstruction pushes Gaussian noise into a Rician distribution with a non-zero mean, and that noise floor biases exactly the low-SNR voxels denoising is meant to help. Given the phase, `-phase` rotates the complex data onto the real axis first, so the noise stays zero-mean Gaussian and the floor is avoided rather than denoised. It can also remove Gibbs ringing (`-degibbs`). The image below shows the b-value 2000 volumes for the 1.5mm [OpenNeuro ds004666](https://openneuro.org/datasets/ds004666/) dataset (the small voxels and high b-values both contribute to the noisy raw images).

<img src=images/anim_EDDEN1p5mm.gif width="420">

## Usage

```
svht_denoise dwi.nii.gz denoised.nii.gz
svht_denoise dwi.nii.gz denoised.nii.gz -mask brain.nii.gz -noise sigma.nii.gz
svht_denoise mag.nii.gz denoised.nii.gz -phase phase.nii.gz
svht_denoise dwi.nii.gz out.nii.gz -degibbs y
```

`-help` lists every option. The interface follows MRtrix3 `dwidenoise` wherever the two share a concept, including [patch size](https://mrtrix.readthedocs.io/en/latest/dwi_preprocessing/denoising.html#patch-size).

**Complex data.** Only the phase *scale* matters; any constant offset is removed along with the background phase. The unit is read off the observed range by default, which is correct for any encoding covering a full turn — scanner integers, degrees, cycles, radians — and the convention chosen is reported on stderr. Data that does *not* span a turn is ambiguous and must be declared with `-phaseunits radians|degrees|turns`. **The output is then signed**, because background voxels scatter about zero instead of piling onto the noise floor.

**Gibbs ringing.** `-degibbs y` denoises then removes ringing, which is the order these belong in; `-degibbs o` removes ringing only, and is the one mode that accepts a 3D image. Do not use either for partial Fourier acquisitions. Neither combines with `-mask`: a mask leaves a hard zero edge, and ringing removal would treat that edge as signal, throwing ringing back into the masked region. `-degibbs o` additionally refuses `-noise`, `-rank`, `-phase`, `-real` and `-extent`, which all describe a denoiser it is not running.

## Compiling

Needs a C99 compiler, libm, pthreads and zlib.

```
git clone --recursive https://github.com/rordenlab/svht_denoise
cd svht_denoise/src
make
make test
```

`--recursive` fetches [zlib-ng](https://github.com/zlib-ng/zlib-ng), the one bundled dependency, which is ~3× faster at writing `.nii.gz`. It is the default **on macOS only**; elsewhere the system zlib is used unless you ask, since it was measured only on Apple Silicon. Also on macOS, the eigensolver's `dsytrd`/`dormtr` and the patch Gram's `cblas_dsyrk` come from Accelerate, worth -36% CPU. Neither is a library to install, and `make` prints what went in — read that line before trusting a timing.

```
make ZLIBNG=1               # zlib-ng off macOS, or force it on
make ZLIBNG=0               # system zlib
make ZLIBNG_ROOT=<dir>      # a prebuilt ZLIB_COMPAT zlib-ng tree
make ACCELERATE=0           # portable kernels (the default off macOS)
make DEGIBBS=0              # omit -degibbs entirely
```

After a plain (non-recursive) clone, `git submodule update --init ../third_party/zlib-ng` from `src` populates it; a build path containing anything outside letters, digits, dot, underscore and hyphen also defeats zlib-ng's own build. Either way the build falls back and says why.

Accelerate changes the arithmetic slightly (~3e-10 relative L2, identical rank map). Byte-identical output **across thread counts** holds in either build, and `-degibbs` is byte-identical across both.

## Benchmark

Apple M4 Pro (14 cores), macOS 26.6, clang 21, mains power, all 14 threads, writing `.nii.gz`. Wall clock and peak resident set size from `/usr/bin/time -l`.

| Dataset                 | method               | Time (s) | Peak RAM (MiB) |
| ----------------------- | -------------------- | -------- | -------------- |
| small 100×100×54 × 36   | dwidenoise           | 2.0      | 227            |
| large 140×140×92 × 297  | dwidenoise           | 1128.4   | 6172           |
| small 100×100×54 × 36   | svht_denoise         | **1.4**  | **157**        |
| large 140×140×92 × 297  | svht_denoise         | **549.9**| **4163**       |
| small 100×100×54 × 36   | mrdegibbs            | 2.9      | 157            |
| large 140×140×92 × 297  | mrdegibbs            | 116.0    | 4100           |
| small 100×100×54 × 36   | svht_denoise degibbs | **1.8**  | 160            |
| large 140×140×92 × 297  | svht_denoise degibbs | **73.3** | 4102           |

On the large series that is **2.05× faster than `dwidenoise` using 2.0 GB less**, and **1.58× faster than `mrdegibbs`**. Both denoisers are given the same already-rotated input, so the rows measure the denoiser rather than the front end; the complex rotation is a shared, untimed NumPy step. `svht_denoise` can do that rotation itself in the same pass (`-phase`), which is how you would really run it.

`svht_denoise` and `dwidenoise` are different estimators, not two spellings of one, so the outputs are not expected to match — `dwidenoise` fits Marchenko-Pastur to estimate sigma, where svht_denoise thresholds directly and never uses a noise level to decide. `-noise` still reports one (sigma = y_med/sqrt(M*mu_beta), the paper's Equation 26), but it is derived after the fact and does not influence the output. `-degibbs` *is* meant to match `mrdegibbs`, and does: relative L2 of 1e-13 to 3e-12 with no voxel differing by more than one float32 ULP, and byte-identical on the test fixture.

The small dataset ships in [benchmark/](benchmark/): 100×100×54 at 2.2 mm, masked, 36 volumes. The large one is [OpenNeuro ds004666](https://openneuro.org/datasets/ds004666/) (EDDEN, 1.5 mm, magnitude and phase, b up to 3010); [benchmark-large/benchmark_large](benchmark-large/benchmark_large) fetches it, runs both pipelines and writes the animation above. Run it with no argument and it fetches then benchmarks; `benchmark_large fetch` only downloads. The benchmark half runs only on files that are present, so a checkout without the data reports what is missing instead of failing.

One caveat on the degibbs rows: this binary links zlib-ng while MRtrix3 links the system zlib, so a compressed write flatters us by ~1.3× that is not the kernel's doing — uncompressed and alternated, `-degibbs o` is 1.52× `mrdegibbs` on CPU at 104² and 1.36× at 140². `scripts/bench.sh` alternates candidates, reverses their order each round, takes a median, and refuses to time on a loaded machine, on battery, or in Low Power Mode.

## Links

**Method**

- Gavish M, Donoho DL. The optimal hard threshold for singular values is 4/√3. *IEEE Trans Inf Theory* 2014;60(8):5040-5053. [DOI](https://ieeexplore.ieee.org/document/6846297) — the threshold implemented here; the [preprint](https://arxiv.org/abs/1305.5870) is bundled with this repository.
- Manjón JV, Coupé P, Concha L, Buades A, Collins DL, Robles M.  Diffusion Weighted Image Denoising Using Overcomplete Local PCA [PMID 24019889](https://pubmed.ncbi.nlm.nih.gov/24019889/) PLoS One. 2013;8:e73021. doi: 10.1371/journal.pone.0073021.
- Veraart J, Novikov DS, Christiaens D, Ades-Aron B, Sijbers J, Fieremans E. Denoising of diffusion MRI using random matrix theory. *NeuroImage* 2016;142:394-406. [PMID 27523449](https://pubmed.ncbi.nlm.nih.gov/27523449/) — MPPCA, the dwidenoise approach.
- Kellner E, Dhital B, Kiselev VG, Reisert M. Gibbs-ringing artifact removal based on local subvoxel-shifts. *Magn Reson Med* 2016;76:1574-1581. — the `-degibbs` method.
- Moeller S, Pisharady PK, Ramanna S, et al. NOise reduction with DIstribution Corrected (NORDIC) PCA. *NeuroImage* 2021;226:117539. — origin of the phase rotation.
- Manzano Patron JP, Moeller S, Andersson JLR, et al. Denoising diffusion MRI: considerations and implications for analysis. *Imaging Neuroscience* 2024. [PMID 40800437](https://pubmed.ncbi.nlm.nih.gov/40800437/) — the "MPPCA\*" pipeline `-phase` reproduces. Datasets, scripts and methods to quantify performance.

**Related Patents**

Tools such as [dwidenoise](https://mrtrix.readthedocs.io/en/dev/reference/commands/dwidenoise.html) and [dwidenoise2](https://github.com/Lestropie/dwidenoise2) rely on algorithms covered by active patents (e.g., [US10698065B2](https://patents.google.com/patent/US10698065B2/en)) and are explicitly restricted to non-commercial research use.

In contrast, `svht_denoise` builds upon earlier, unencumbered frameworks by Manjón et al. (2013) and Gavish & Donoho (2014). The algorithm calculates noise thresholds directly without relying on MPPCA noise estimation, and is implemented as an independent clean-room build. Consequently, the authors believe `svht_denoise` is not subject to these patent restrictions. However, commercial users should independently evaluate their legal compliance, as the authors provide no formal legal warranties.

**Software**

- [MRtrix3](https://www.mrtrix.org/) — [dwidenoise](https://mrtrix.readthedocs.io/en/dev/reference/commands/dwidenoise.html) and [mrdegibbs](https://mrtrix.readthedocs.io/en/dev/reference/commands/mrdegibbs.html), the reference implementations benchmarked here. Note `dwidenoise` is restricted to non-commercial research use.
- [dwidenoise2](https://github.com/Lestropie/dwidenoise2) — a newer method, also non-commercial only.
- [DIPY](https://dipy.org/) — `benchmark/dipy_mppca.py` runs its MPPCA for comparison.
- [niimath](https://github.com/rordenlab/niimath) — renders the animation frames.
- [zlib-ng](https://github.com/zlib-ng/zlib-ng) — bundled, for faster `.nii.gz` writing.
- [OpenNeuro ds004666](https://openneuro.org/datasets/ds004666/) — the large benchmark dataset.

## Licensing

MPL-2.0, matching the core MRtrix3 codebase. Unlike `dwidenoise`, which carries an additional NYU/Antwerp notice restricting it to non-commercial research, svht_denoise has no such restriction.

- **Algorithm.** An independent clean-room C implementation from the Gavish and Donoho paper's formulations — Equation 4's threshold with Equation 11's closed form and an exactly bisected Marchenko-Pastur median — not derived from its GPL-3.0 reference MATLAB.
- **CLI.** Flag names follow `dwidenoise` conventions. No MRtrix3 code is in the denoiser.
- **Gibbs ringing.** `src/mrdegibbs/dg.c` **is adapted from** MRtrix3's `cmd/mrdegibbs.cpp`, by Ben Jeurissen and J-Donald Tournier, © 2008-2025 the MRtrix3 contributors, MPL-2.0. The `-axes` option is not carried over.
- **FFT.** `src/mrdegibbs/dg_fft.c` **is adapted from** Eigen's `unsupported/Eigen/src/FFT/ei_kissfft_impl.h`, © 2009 Mark Borgerding, and MPL-2.0 as distributed by Eigen. It derives in turn from [kissfft](https://github.com/mborgerding/kissfft), © 2003-2009 Mark Borgerding, which upstream is BSD-3-Clause. It is the same transform `mrdegibbs` reaches through `Eigen::FFT`; MRtrix3 links neither FFTW nor a BLAS.
- **Complex data.** The rotation follows NORDIC's phase-removal stage, with Steen Moeller's `NIFTI_COMP_to_REAL.m` read as a specification; no code was translated. Reading phase units off the data range follows [niimath](https://github.com/rordenlab/niimath) (BSD-2-Clause).
- **zlib-ng** © 1995-2024 Jean-loup Gailly and Mark Adler, zlib licence ([text](packaging/LICENSE-zlib-ng.txt)), statically linked into every macOS build (`ZLIBNG=0` links the system zlib and embeds nothing).
- **NIfTI I/O.** `src/nifti_io.c` is vendored public-domain code by Robert W Cox, Mark Jenkinson, Rick Reynolds and Chris Rorden, via [niimath](https://github.com/rordenlab/niimath). It is the largest third-party component here and is not modified.
- **AI assistance.** Generative AI tools were used during translation and refactoring.

The phase rotation's two deliberate divergences and its 7-tap convolution derivation are documented at the head of `src/dn_phase.c`; the degibbs design record is at the head of `src/mrdegibbs/dg.c`.

## macOS release packages

`make macos-release VERSION=<tag>` builds a signed, notarized, stapled arm64 installer targeting macOS 14 or newer, into `dist/`. It needs two certificates from one developer account — "Developer ID Application" for the executable and "Developer ID Installer" for the `.pkg` — plus credentials stored once with `make macos-notary-profile APPLE_ID=... TEAM_ID=...` (which prompts for an [app-specific password](https://appleid.apple.com), not your Apple ID password).

Releases are arm64 only, and `MACOSX_DEPLOYMENT_TARGET` must not go below 13.3, where Accelerate's new LAPACK symbols first exist. `make macos-pkg-adhoc` builds an unsigned package for local testing; `scripts/verify_macos_pkg.sh <pkg>` inspects one without installing it — note it *runs* the packaged executable, so use it on packages you built, not on one from an untrusted source. The Makefile's release section has the detail.

