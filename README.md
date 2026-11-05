## Introduction

svht_denoise removes thermal noise from multi-volume NIfTI datasets such as diffusion-weighted MRI. It decomposes local image patches by SVD and suppresses the noise components with the singular-value rules of Gavish and Donoho, so no noise level has to be supplied. By default:
- **Shrinkage:** each component is shrunk by their Frobenius-optimal rule (2017), against a noise map the tool estimates first with their median estimator. It is not truncated by their optimal hard threshold (2014).
- **Patches:** they are spherical in mm, centred on every second voxel, and averaged where they overlap (Manjón et al. 2013).
- **Magnitude input:** it is variance-stabilised (Foi 2011), so the Rician noise floor is removed rather than denoised.

It works on magnitude volumes, and on complex data when the phase is also available. Magnitude reconstruction pushes Gaussian noise into a Rician distribution with a non-zero mean, and that noise floor biases exactly the low-SNR voxels denoising is meant to help. Given the phase, `-phase` rotates the complex data onto the real axis first, so the noise stays zero-mean Gaussian and the floor is avoided rather than denoised. It can also remove Gibbs ringing (`-degibbs`). The image below shows the 180 volumes with b > 1500 (the b≈2000 and b≈3000 shells) from the 1.5 mm [OpenNeuro ds004666](https://openneuro.org/datasets/ds004666/) dataset (the small voxels and high b-values both contribute to the noisy raw images).

<img src=images/anim_EDDEN1p5mm.gif width="420">

## Usage

```
svht_denoise dwi.nii.gz denoised.nii.gz
svht_denoise dwi.nii.gz denoised.nii.gz -mask brain.nii.gz -noise sigma.nii.gz
svht_denoise mag.nii.gz denoised.nii.gz -phase phase.nii.gz
svht_denoise dwi.nii.gz out.nii.gz -degibbs y
svht_denoise dwi.nii.gz out.nii.gz -filter optthresh -aggregator exclusive -shape cube -vst n   # the original pipeline
```

**Magnitude or not.** The variance-stabilising transform (`-vst`) assumes magnitude data. It is off with `-phase` (which refuses `-vst y`), and by default for any input containing a negative value, since that cannot be a magnitude image. A series already rotated to real elsewhere is caught by its negative voxels; pass `-vst n` for any non-magnitude input that has none (a masked one may not). `-noise_dof` sets the number of receive channels for sum-of-squares reconstructions (default 1, Rician).

**Fewer than 3 noise columns** (volumes minus demeaned shells, e.g. a 2-volume series) cannot give a noise level, so the run falls back to `-filter optthresh` without `-vst`, and refuses `-noise`.

**Per-shell demeaning** (`-demean y`, reading the input's `.bval` or `-bval`; shells are clustered as in MRtrix3, sorted neighbours less than 80 s/mm² apart sharing a shell) is available, but off by default: measured, it is neutral at best and up to 6.7% worse RMSE under the defaults.

`-help` lists every option. The interface follows MRtrix3 `dwidenoise` wherever the two share a concept, including [patch size](https://mrtrix.readthedocs.io/en/latest/dwi_preprocessing/denoising.html#patch-size). Patch size starts from k, by default the smallest odd k with k³ > the number of volumes (`-extent k` overrides). A cube is k³ voxels. A sphere, measured in mm from the header spacing, is the smallest ball of whole distance shells holding at least k³ voxels, so about 1.24k across, and it can overshoot by up to one shell (147 voxels at k = 5, not 125): the voxel count of a radius is a lattice-point count with no closed form, roughly (4/3)πr³ / voxel volume. The run summary prints the radius and voxel count.

**Complex data.** Only the phase *scale* matters; any constant offset is removed along with the background phase. The unit is read off the observed range by default, which is correct for any encoding covering a full turn — scanner integers, degrees, cycles, radians — and the convention chosen is reported on stderr. Data that does *not* span a turn is ambiguous and must be declared with `-phaseunits radians|degrees|turns`. **The output is then signed**, because background voxels scatter about zero instead of piling onto the noise floor. One exception: `-degibbs` assumes magnitude data and truncates negative input to zero, so combining it with `-phase`, or with any other signed input, gives up that property; the run warns when it does. `-real` writes the rotated series before either stage and is always signed.

**Gibbs ringing.** `-degibbs y` denoises then removes ringing, which is the order these belong in; `-degibbs o` removes ringing only, and is the one mode that accepts a 3D image. For partial-Fourier acquisitions, add `-pF 0.875` (7/8) or `-pF 0.75` (6/8); it assumes y is the phase-encode direction, and 6/8 additionally needs an even y dimension. Only those two factors are implemented, and any value other than those two and `1.0` (full k-space, the default) is refused rather than rounded to the nearest — the pipeline differs per factor instead of varying continuously. The two are not equally exact: 6/8 reduces to the ordinary method applied unchanged to the odd and even y-columns, where 7/8 is a resampling construction whose interpolation choice is ours, and reference implementations differ from each other by more than we differ from any of them. Negative input is truncated to zero first, because the method assumes magnitude data and a flat zero baseline is what its total-variation search wants in background regions; truncating rather than rectifying keeps zero-mean noise from becoming a Rician-like positive floor. That is a no-op on a magnitude image. Ringing correction still undershoots, so the *output* can be negative — except at `-pF 0.875`, which clamps its own output as its reference does. In-plane dimensions must be at least 5, and at least 5 again after a `-pF` factor splits them. `-degibbs` does not combine with `-mask` in either mode, at any factor, including full k-space: a mask leaves a hard zero edge, and ringing removal would treat that edge as signal, throwing ringing back into the masked region. `-degibbs o` additionally refuses `-noise`, `-rank`, `-phase`, `-real` and `-extent`, which all describe a denoiser it is not running.

### Multi-echo (tensor) denoising

```
svht_denoise echo-1.nii.gz den.nii.gz -tensor gd -echo echo-2.nii.gz -echo echo-3.nii.gz
```

`-tensor gd` denoises the echoes of one multi-echo series jointly and writes `den_echo-1.nii.gz`, `den_echo-2.nii.gz`, and so on. Each window is a space × time × echo **tensor**, not a matrix with echoes concatenated onto time. It is decomposed one mode at a time in the recursive scheme of Olesen et al. (2022): unfold, keep the components above the noise bulk edge, pass the reduced tensor on to the next mode, and rebuild on the way back. By default the last mode's kept values get the Gavish–Donoho Frobenius shrinkage.

In tensor mode `-filter` takes the same three names as in matrix mode: `optshrink` (the default) as above, `truncate` (the same cuts with nothing shrunk), or `optthresh`. Tensor `optthresh` replaces each mode's edge cut with Gavish and Donoho's 2014 optimal hard threshold for a **known** noise level, s > λ*(β) σ √n, where n is the unfolding's longer side and β the ratio of its sides, and keeps what survives unshrunk. It cuts every mode higher than the edge, since λ*(β) > 1 + √β. That is a different rule from matrix `optthresh`, which uses their unknown-noise threshold, ω(β) times each patch's median singular value. The median needs many noise singular values, and the echo mode has only as many as there are echoes, most of them signal. So tensor `optthresh` reads the same per-echo noise as the other two filters. ME-PCA has no hard-threshold option, so `optthresh` is outside the parity check below; `make test` pins it with unit cases on either side of the threshold and an anchor.

The noise level comes from the **same Gavish–Donoho median-estimator map the matrix denoiser reads**, estimated from each echo's own series. The tensor kernel assumes one noise variance across the whole tensor, so each echo is first put in units of **its own** noise and mapped back afterwards. Magnitude echoes are variance-stabilized against their map (`-vst`, Rician by default or non-central χ with `-noise_dof`, with matrix mode's re-estimation round and exact-unbiased inverse); rotated complex echoes are divided by it. The run summary prints each echo's median σ. No Marchenko–Pastur fit is used. Olesen et al. pool one over every mode's spectrum, and that part is deliberately not implemented. Each mode is cut at the noise bulk edge σ(√m + √c), the same edge matrix-mode `optshrink` discards below. `-tensor_sigma2` gives one known noise variance for every echo instead. Tensor mode shares `-phase`, `-phaseunits`, `-filter`, `-vst`, `-noise_dof`, `-preserve_noise_bias` and `-aggregator` with matrix mode, with the meanings above; it refuses the matrix denoiser's other options, its own are all `-tensor_*`, and `-help` lists them.

- **Noise across echoes.** Thermal noise is the same at every echo in k-space, but it need not be the same in the reconstructed images that the denoiser sees. On ds006131 (Siemens, MB6 with in-plane GRAPPA, the same bandwidth on every echo), the noise-only scans put echo 5's noise at 0.51 of echo 1's. That ratio ranges from 0.11 to 0.77 across the brain, which points to echo-specific reconstruction rather than physics. A single level pooled over the echoes read 29% low for echo 1 and 39–43% high for echo 5. Each echo's own map read within 1–6% on every echo. `-tensor_noise pooled` keeps one level for all echoes, for reconstructions where they do share one; on that run it under-denoised the first echo, over-cut the last, and, with the VST, biased echo 5's mean by −8%.
- **Complex data.** `-phase`, given once per echo in echo order, takes each echo's phase and rotates that echo onto the real axis with it, exactly as `-phase` does for a single series in matrix mode. The output is then signed.
- **Shrinker dimensions.** Under `-filter optshrink`, the last mode's kept values are shrunk by the Gavish–Donoho 2017 rule, and `-tensor_shrinker` picks the dimensions it is given. The default, `olesen`, is the unfolding less the kept rank, (m − r) × (c − r), as Olesen et al. and ME-PCA use. `gd` is the observed m × c, as Gavish and Donoho derive it and as matrix mode uses. Algebraically, `olesen` is `gd` plus 4rσ² on every kept squared value: a component just above the edge keeps a weight of about 2√r/(√m + √c) where `gd`'s fades to zero. Neither form is derived for this setting. Gavish and Donoho assume noise iid over the matrix being shrunk, and the last level's matrix is the projection kept by the earlier levels, whose directions were chosen from the same noise. So the default is set by measurement rather than parity.
  - **Benchmark.** A 16-slice slab of ds006131 (110×110×16, 243 volumes, 5 echoes). The truth is matrix-mode output, Rician noise follows each echo's noRF σ map at 1× and 2×, and errors are relative RMSE over the brain mask.
  - **Default mode order:** the last mode is the echo mode, and its five components all sit far above the edge. So `gd`, `olesen` and `truncate` agree to within 0.03% at either noise level, magnitude or complex. On the real slab the two shrinkers' outputs differ by 4e-6–4e-5 relative L2, with identical ranks.
  - **Time as the last mode** (`-tensor_order space,echo,time`): there the shrunk components do sit near the edge. `olesen` beats `gd` at 1× (0.011730 against 0.011788, from 0.0487 in the input) and at 2× (0.017990 against 0.018012). `truncate` beats both at 1× (0.011519) and loses to both at 2× (0.018038). The shrinker that removes less does better here, which is the direction the iid assumption's failure predicts.
  - **Verdict:** `olesen` stays the default, and `gd` is there for comparison. `-filter optthresh` cuts harder and lost at 1× in both orders (0.013528 magnitude, 0.015610 complex).

- **Windows.** The default is a local 5×5×5 window. `-tensor_window 4,4` slides in x and y and spans every slice, as the [ME-PCA](https://github.com/tsalo/ME-PCA) package's `run_TensorMPPCA` does; the two are not interchangeable. The window rules are ME-PCA's, not the matrix denoiser's. A window sits wherever it fits inside the image, with no edge shifting, and overlapping estimates are averaged uniformly at `-tensor_stride` (default 1); `-aggregator exclusive` instead writes only each window's center voxel (its center column for an in-plane window), which is ME-PCA's center assignment, and `-aggregator gaussian` is matrix mode's only. A voxel that no window covers keeps its input value (to float rounding, through the noise normalization and back) and reads NaN in every map. `-mask` selects the windows that touch it, and each such window writes every voxel it covers, so voxels outside the mask are **not** zeroed. The noise map is not masked, as in matrix mode.
- **Input.** Every echo is a 4D magnitude series on the same grid, with the same number of volumes and the same TR, and all of this is checked. Estimating the noise needs at least 3 volumes per echo. Negative values are refused: ME-PCA takes `abs()` of its input, which would silently destroy a signed series. Complex data come in as magnitude plus one `-phase` per echo.
- **Outputs.** `-noise` is the noise standard deviation, one volume per echo, in that echo's units. `-rank` has one volume per mode, in `-tensor_order`, and its header says which. `-snr_gain` is Olesen et al.'s SNR-gain heuristic, computed from the mean ranks; it reads NaN where nothing was kept.
- **Validation.** At a known noise variance, `tools/tensor_parity.py` checks the kernel and the CLI against ME-PCA (float64) on synthetic data; that runs only ME-PCA's recursive pass and window rules, never its noise fit. Across 36 kernel cases and 11 image cases (local and in-plane windows, with and without shrinkage, demeaning, strides, masks, center assignment and mode orders), the two agree to ≤ 2e-11 relative L2 per tensor and ≤ 6e-8 per image, with identical ranks and identical uncovered voxels. On the same data with noise of SD 20, the per-echo maps' medians read 20.4–20.8, and the error to the truth (0.00926, from 0.0375 in the input) is no worse than at the true variance (0.00939). `make test` pins the geometry, the refusals and numerical anchors, and runs the same oracle comparison, without needing Python.
- **Signs under demeaning.** With `-tensor_demean y`, later levels take means across the reshaped V·S of the previous one, so the result depends on each singular vector's sign, which LAPACK leaves undefined: on the parity cases ME-PCA's own output moved by up to 0.4% with it. svht_denoise fixes the sign (each V column's largest-magnitude entry is positive), which makes a demeaned run reproducible across platforms; the parity script applies the same convention to ME-PCA. Without demeaning, the sign cannot matter.
- **Speed and memory.** A 110×110×72 run of 243 volumes and 5 echoes (ds006131, 764,000 windows of 5×5×5) took 13.5 min of wall time and 323 CPU-minutes on 28 threads for magnitude, and 12.1 min for complex, with a peak of 8.2 GiB, since every echo is held in and out at once. Matrix mode on the same echoes, one after another, took 11.7 min and 158 CPU-minutes, peaking at 1.7 GiB per echo.

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
| small 100×100×54 × 36   | dwidenoise2          | crash¹   | —              |
| large 140×140×92 × 297  | dwidenoise2          | 388.4    | 8471           |
| small 100×100×54 × 36   | svht_denoise         | **1.2**  | **173**        |
| large 140×140×92 × 297  | svht_denoise         | **155.1**| **4263**       |
| small 100×100×54 × 36   | mrdegibbs            | 2.9      | 157            |
| large 140×140×92 × 297  | mrdegibbs            | 116.0    | 4100           |
| small 100×100×54 × 36   | svht_denoise degibbs | **1.8**  | 160            |
| large 140×140×92 × 297  | svht_denoise degibbs | **73.3** | 4102           |

On the large series that is **7.3× faster than `dwidenoise` using 1.9 GB less**, **2.5× faster than `dwidenoise2` using half the memory**, and **1.58× faster than `mrdegibbs`**.

- **Configuration.** All rows use each tool's defaults, single runs. dwidenoise2 was built natively from its pinned MRtrix3 commit and given the gradient table in MRtrix format (`-grad`); its FSL reader rejected this repository's bvec. svht_denoise's VST stays off here because the input is already rotated to real, and so is signed. dwidenoise2 has no such check and applies its magnitude VST to it.
- **¹ dwidenoise2 segfaults on the small series.** That series is brain-masked, and dwidenoise2 crashes on a large constant background; filling the background with noise lets it run. On that noise-filled copy dwidenoise2 took 3.0 s and 258 MiB, against svht_denoise's 1.5 s and 171 MiB (median of 3).
- **Shared input.** Every denoiser is given the same already-rotated input, so the rows measure the denoiser rather than the front end; the complex rotation is a shared, untimed NumPy step. `svht_denoise` can do that rotation itself in the same pass (`-phase`), which is how you would really run it.

**Against dwidenoise2.** On the large series as magnitude, with defaults for both and the gradient table given to dwidenoise2: svht_denoise took **145 s and 4.2 GiB peak**, dwidenoise2 375 s and 8.3 GiB (single runs, `.nii` output). For accuracy, use a synthetic truth: the `dwidenoise` output on `benchmark/` plus Rician noise of known σ, scored inside the brain. Error (RMSE) and mean bias at b=2000:

| σ   | svht_denoise original | svht_denoise default | dwidenoise | dwidenoise2 |
| --- | --------------------- | -------------------- | ---------- | ----------- |
| 10  | 7.060 / +1.8          | 6.363 / −0.8         | 6.903 / +1.8 | **6.264** / −0.3 |
| 30  | 18.801 / +13.4        | **13.493** / −0.2    | 19.629 / +13.4 | 13.986 / −0.4 |
| 60  | 44.243 / +41.9        | **25.202** / +8.9    | 46.606 / +42.0 | 25.851 / +7.4 |

The truth image is itself a denoiser's output and keeps some residual noise, which flatters the low-σ row toward whichever method removes least; read that row with care.

`svht_denoise` and `dwidenoise` are different estimators, not two spellings of one, so the outputs are not expected to match. `dwidenoise` fits Marchenko-Pastur to estimate sigma. svht_denoise estimates it with Gavish and Donoho's median estimator (sigma = y_med/sqrt(M*mu_beta), the 2014 paper's Equation 26), from larger patches on a coarse grid; under `-filter optthresh` it thresholds without using sigma at all. `-degibbs` *is* meant to match `mrdegibbs`, and does: relative L2 of 1e-13 to 3e-12 with no voxel differing by more than one float32 ULP, and byte-identical on the test fixture. That is a full-k-space claim only — `mrdegibbs` has no `-pF`, so nothing in it extends to the partial-Fourier paths. It also assumes an input with no negative voxel, which every magnitude image satisfies and a *degibbsed* one does not: undershoot means `-degibbs` output typically contains negatives, and feeding that back in truncates them where `mrdegibbs` would not. Chaining the stage twice, or running it on an already-signed series, is therefore outside the parity claim by construction.

The small dataset ships in [benchmark/](benchmark/): 100×100×54 at 2.2 mm, masked, 36 volumes. The large one is [OpenNeuro ds004666](https://openneuro.org/datasets/ds004666/) (EDDEN, 1.5 mm, magnitude and phase, b up to 3010); [benchmark-large/benchmark_large](benchmark-large/benchmark_large) fetches it, runs both pipelines and writes the animation above. Run it with no argument and it fetches then benchmarks; `benchmark_large fetch` only downloads. The benchmark half runs only on files that are present, so a checkout without the data reports what is missing instead of failing.

One caveat on the degibbs rows: this binary links zlib-ng while MRtrix3 links the system zlib, so a compressed write flatters the *wall-clock* rows above by ~1.3× that is not the kernel's doing. Measured uncompressed and alternated, on **CPU time** rather than wall, `-degibbs o` is 1.52× `mrdegibbs` at 104² and 1.36× at 140² — a different column from the table, so the two are not meant to reconcile arithmetically. `scripts/bench.sh` alternates candidates, reverses their order each round, takes a median, and refuses to time on a loaded machine, on battery, or in Low Power Mode.

## Links

**Method**

- Gavish M, Donoho DL. The optimal hard threshold for singular values is 4/√3. *IEEE Trans Inf Theory* 2014;60(8):5040-5053. [DOI](https://ieeexplore.ieee.org/document/6846297) — the `-filter optthresh` threshold and the median noise estimator; the [preprint](https://arxiv.org/abs/1305.5870) is bundled with this repository.
- Gavish M, Donoho DL. Optimal shrinkage of singular values. *IEEE Trans Inf Theory* 2017;63(4):2137-2152. [preprint](https://arxiv.org/abs/1405.7511) — the default `-filter optshrink`.
- Foi A. Noise estimation and removal in MR imaging: the variance-stabilization approach. *IEEE ISBI* 2011:1809-1814. — the `-vst` transform and its exact-unbiased inverse.
- Manjón JV, Coupé P, Concha L, Buades A, Collins DL, Robles M.  Diffusion Weighted Image Denoising Using Overcomplete Local PCA [PMID 24019889](https://pubmed.ncbi.nlm.nih.gov/24019889/) PLoS One. 2013;8:e73021. doi: 10.1371/journal.pone.0073021. — the overlapping-patch averaging (`-aggregator`).
- Veraart J, Novikov DS, Christiaens D, Ades-Aron B, Sijbers J, Fieremans E. Denoising of diffusion MRI using random matrix theory. *NeuroImage* 2016;142:394-406. [PMID 27523449](https://pubmed.ncbi.nlm.nih.gov/27523449/) — MPPCA, the dwidenoise approach.
- Kellner E, Dhital B, Kiselev VG, Reisert M. Gibbs-ringing artifact removal based on local subvoxel-shifts. *Magn Reson Med* 2016;76:1574-1581. — the `-degibbs` method.
- Lee HH, Novikov DS, Fieremans E. Removal of Partial Fourier-Induced Gibbs (RPG) Ringing artifacts in MRI. *Magn Reson Med* 2021;86:2733-2750. [PMID 34227142](https://pubmed.ncbi.nlm.nih.gov/34227142/) — the RPG method used by `-degibbs -pF` for 7/8 and 6/8 partial Fourier data.
- Moeller S, Pisharady PK, Ramanna S, et al. NOise reduction with DIstribution Corrected (NORDIC) PCA. *NeuroImage* 2021;226:117539. — origin of the phase rotation.
- Manzano Patron JP, Moeller S, Andersson JLR, et al. Denoising diffusion MRI: considerations and implications for analysis. *Imaging Neuroscience* 2024. [PMID 40800437](https://pubmed.ncbi.nlm.nih.gov/40800437/) — the "MPPCA\*" pipeline `-phase` reproduces. Datasets, scripts and methods to quantify performance.
- Olesen JL, Ianus A, Østergaard L, Shemesh N, Jespersen SN. Tensor denoising of multidimensional MRI data. *Magn Reson Med* 2022. [doi:10.1002/mrm.29478](https://doi.org/10.1002/mrm.29478) — the recursive mode-wise scheme of `-tensor gd`, whose own Marchenko–Pastur noise fit is not used.

**Related Patents**

Tools such as [dwidenoise](https://mrtrix.readthedocs.io/en/dev/reference/commands/dwidenoise.html) and [dwidenoise2](https://github.com/Lestropie/dwidenoise2) rely on algorithms covered by active patents (e.g., [US10698065B2](https://patents.google.com/patent/US10698065B2/en)) and are explicitly restricted to non-commercial research use.

In contrast, `svht_denoise` builds upon earlier, unencumbered frameworks by Manjón et al. (2013) and Gavish & Donoho (2014). The algorithm estimates the noise level with the Gavish–Donoho median estimator rather than MPPCA's Marchenko–Pastur fit, and is implemented as an independent clean-room build. Consequently, the authors believe `svht_denoise` is not subject to these patent restrictions. However, commercial users should independently evaluate their legal compliance, as the authors provide no formal legal warranties. The tensor mode (`-tensor gd`) takes its noise level from the same Gavish–Donoho median estimator; it does not implement the Marchenko–Pastur fit of Olesen et al.'s tensor MP-PCA.

**Software**

- [MRtrix3](https://www.mrtrix.org/) — [dwidenoise](https://mrtrix.readthedocs.io/en/dev/reference/commands/dwidenoise.html) and [mrdegibbs](https://mrtrix.readthedocs.io/en/dev/reference/commands/mrdegibbs.html), the reference implementations benchmarked here. Note `dwidenoise` is restricted to non-commercial research use.
- [dwidenoise2](https://github.com/Lestropie/dwidenoise2) — a newer method, also non-commercial only.
- [DIPY](https://dipy.org/) — `benchmark/dipy_mppca.py` runs its MPPCA for comparison.
- [niimath](https://github.com/rordenlab/niimath) — renders the animation frames.
- [zlib-ng](https://github.com/zlib-ng/zlib-ng) — bundled, for faster `.nii.gz` writing.
- [OpenNeuro ds004666](https://openneuro.org/datasets/ds004666/) — the large benchmark dataset.

## Licensing

MPL-2.0, matching the core MRtrix3 codebase. Unlike `dwidenoise`, which carries an additional NYU/Antwerp notice restricting it to non-commercial research, svht_denoise has no such restriction.

- **Algorithm.** An independent clean-room C implementation from the Gavish and Donoho paper's formulations — Equation 4's threshold with Equation 11's closed form and an exactly bisected Marchenko-Pastur median — not derived from its GPL-3.0 reference MATLAB. The 2017 shrinker, the overlapping-patch averaging and the variance-stabilising transform are likewise written from their papers. [dwidenoise2](https://github.com/Lestropie/dwidenoise2), whose author co-develops this tool, was read as a design reference for which refinements to adopt; none of its code is used.
- **CLI.** Flag names follow `dwidenoise` conventions. No MRtrix3 code is in the denoiser.
- **Eigensolver.** `src/dn_eig.c`'s `tred2`/`tql2` **are adapted from** [niimath](https://github.com/rordenlab/niimath)'s `src/tensor.c` (BSD-2-Clause), where the same routines are specialised to n = 3; they are generalised here to a runtime *n*, with an iteration cap and a status return added. They are in turn the public-domain JAMA translation of the EISPACK Fortran and the Algol procedures of Bowdler, Martin, Reinsch and Wilkinson. `dn.h`'s checked size arithmetic, `dn_nii.c`'s datatype conversion and `dn_nii.h`'s grid comparison follow niimath too.
- **Gibbs ringing.** `src/mrdegibbs/dg.c` **is adapted from** MRtrix3's `cmd/mrdegibbs.cpp`, by Ben Jeurissen and J-Donald Tournier, © 2008-2025 the MRtrix3 contributors, MPL-2.0. The `-axes` option is not carried over.
- **Partial Fourier.** The `-pF` RPG correction in that same file is *not* MRtrix3's — `mrdegibbs` has no equivalent. It is an independent clean-room implementation from Lee et al. 2021 and a written procedural description of the method. Neither the authors' own toolbox, which is licensed for non-commercial use, nor TORTOISE's GPL implementation was read or translated; both were used only as black boxes to compare output against.
- **FFT.** `src/mrdegibbs/dg_fft.c` **is adapted from** Eigen's `unsupported/Eigen/src/FFT/ei_kissfft_impl.h`, © 2009 Mark Borgerding, and MPL-2.0 as distributed by Eigen. It derives in turn from [kissfft](https://github.com/mborgerding/kissfft), © 2003-2009 Mark Borgerding, which upstream is BSD-3-Clause. It is the same transform `mrdegibbs` reaches through `Eigen::FFT`; MRtrix3 links neither FFTW nor a BLAS.
- **Complex data.** `src/dn_phase.c`'s rotation follows NORDIC's phase-removal stage, with Steen Moeller's `NIFTI_COMP_to_REAL.m` read as a specification; no code was translated into the shipped tool. Reading phase units off the data range follows [niimath](https://github.com/rordenlab/niimath) (BSD-2-Clause).
- **Tensor mode.** `src/dn_tensor.c` and `src/dn_tensor_run.c` are written in C from Olesen et al. (2022)'s recursive scheme and the Gavish–Donoho 2017 shrinker; the noise map is this tool's own (`src/dn_sigma.c`). The [ME-PCA](https://github.com/tsalo/ME-PCA) Python package (`python-package` branch, LGPL-2.1) was read as the specification of the recursion's conventions and the window rules, and is run as a black box by `tools/tensor_parity.py`, which never reaches its noise fit. That package is itself a translation of Olesen's GPL-3.0 MATLAB [Tensor-MP-PCA](https://github.com/Neurophysics-CFIN/Tensor-MP-PCA). The C code has its own structure (SVD through the smaller Gram, a fixed sign convention), but **its provenance relative to those licences has not been formally reviewed.**
- **Benchmark harness.** `benchmark-large/mrtrix.py` **is a NumPy port of** `NIFTI_COMP_to_REAL.m` (Steen Moeller, 2023, shipped with the EDDEN code), itself `NIFTI_NORDIC.m` with the denoising removed. It exists only to build the comparison arm and is not part of the tool or of any release artefact.
- **zlib-ng** © 1995-2024 Jean-loup Gailly and Mark Adler, zlib licence ([text](packaging/LICENSE-zlib-ng.txt)), statically linked into every macOS build (`ZLIBNG=0` links the system zlib and embeds nothing).
- **NIfTI I/O.** `src/nifti_io.c` is vendored public-domain code by Robert W Cox, Mark Jenkinson, Rick Reynolds and Chris Rorden, via [niimath](https://github.com/rordenlab/niimath). It is the largest third-party component here and is not modified.
- **AI assistance.** Generative AI tools were used during translation and refactoring.

The phase rotation's two deliberate divergences and its 7-tap convolution derivation are documented at the head of `src/dn_phase.c`; the degibbs design record is at the head of `src/mrdegibbs/dg.c`.

## macOS release packages

`make macos-release VERSION=<tag>` builds a signed, notarized, stapled arm64 installer targeting macOS 14 or newer, into `dist/`. It needs two certificates from one developer account — "Developer ID Application" for the executable and "Developer ID Installer" for the `.pkg` — plus credentials stored once with `make macos-notary-profile APPLE_ID=... TEAM_ID=...` (which prompts for an [app-specific password](https://appleid.apple.com), not your Apple ID password).

Releases are arm64 only, and `MACOSX_DEPLOYMENT_TARGET` must not go below 13.3, where Accelerate's new LAPACK symbols first exist. `make macos-pkg-adhoc` builds an unsigned package for local testing; `scripts/verify_macos_pkg.sh <pkg>` inspects one without installing it — note it *runs* the packaged executable, so use it on packages you built, not on one from an untrusted source. The Makefile's release section has the detail.
