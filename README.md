## Introduction

svht_denoise provides a principled, data-driven approach for noise reduction in multi-volume NIfTI datasets, such as diffusion-weighted MRI. The tool applies Singular Value Decomposition (SVD) over local image patches to identify and truncate noise components using the optimal Singular Value Hard Thresholding (SVHT) framework derived by Gavish and Donoho (2014).

The current version operates on magnitude-reconstructed volumes. Because magnitude conversion shifts thermal Gaussian noise into a non-Gaussian (Rician/Non-Central Chi) distribution, processing raw complex data (real/imaginary or magnitude/phase) remains a planned enhancement to eliminate low-SNR background bias.

<img src=images/anim_compare.gif width="300">

## Building

Requires a C99 compiler on a POSIX system (macOS, Linux, WSL) and zlib. The only other libraries are libm and pthreads, both standard on those platforms; there are no third-party dependencies.

```
git clone https://github.com/rordenlab/svht_denoise
cd svht_denoise/src
make
```

This writes the `svht_denoise` executable into `src`. Run `./svht_denoise -help` for the full list of options.

```
./svht_denoise dwi.nii.gz denoised.nii.gz
./svht_denoise dwi.nii.gz denoised.nii.gz -mask brain.nii.gz -noise sigma.nii.gz
```

## Licensing

The underlying algorithm for optimal singular value hard thresholding (SVHT) is based on the work of Gavish and Donoho (2014).

 - Algorithm & Theory: The method computes the unknown-noise-level threshold of Equation 4, tau = omega(beta) * y_med, where omega(beta) = lambda_star(beta) / sqrt(mu_beta) is assembled from the closed form of Equation 11 and an exactly bisected Marchenko-Pastur median, rather than from the cubic approximation to omega given in Equation 5. While reference MATLAB scripts for the method were distributed under GPL-3.0, this project is an independent clean-room C implementation built directly from the paper's mathematical formulations rather than derived from the original reference code.
 - CLI Design: The command-line interface, flag names, and operational parameters were intentionally structured to follow the conventions of dwidenoise from MRtrix3 wherever the two tools share a concept (`-mask`, `-noise`, `-extent`, `-nthreads`, `-quiet`). No code from MRtrix3 was incorporated.
 - AI Assistance: Generative AI tools were used during translation and refactoring.
 - Prior Work: Applying random matrix theory to denoise diffusion MRI was established by Veraart, Novikov, Fieremans and colleagues (2016), whose [dwidenoise](https://mrtrix.readthedocs.io/en/dev/reference/commands/dwidenoise.html) remains the reference implementation and the benchmark this tool is measured against. svht_denoise is an independent implementation of a different estimator, and shares no code with it: dwidenoise fits the Marchenko-Pastur distribution to the patch eigenspectrum to estimate the noise level sigma, which then sets the cutoff, whereas svht_denoise applies the Gavish and Donoho hard threshold directly, at omega(beta) times the median singular value, so the noise level never enters the retain-or-discard decision. An estimate is still available on request: the optional `-noise` map reports sigma = y_med / sqrt(M * mu_beta), Equation 26 of the paper, but it is derived after the fact and does not influence the denoised output.

This executable and its source code are licensed under the Mozilla Public License 2.0 (MPL-2.0), matching the permissive license used by the core MRtrix3 codebase. Note that dwidenoise itself is not covered by that permissive license: it carries an additional notice from New York University and the University of Antwerp restricting it to non-commercial research and excluding clinical care. svht_denoise carries no such restriction.

## Benchmarking

A sample diffusion dataset is provided to evaluate the time and peak memory required to process a provided diffusion scan. This demonstrates that the performance on these metrics is similar to [dwidenoise](https://mrtrix.readthedocs.io/en/dev/reference/commands/dwidenoise.html). However, note these are different algorithms and the outputs are not expected to be identical. Further, the input image is intentionally small to minimize the size of this repository: 100x100x54 voxels at 2.2mm isotropic resolution, masked, with 36 volumes (a b=0 image, 30 half-sphere b=2000 directions, and five trailing b=0 images).

Run these from `src`, after `mkdir -p ../benchmark/svht_denoise ../benchmark/dwidenoise`:

```
/usr/bin/time -l ./svht_denoise ../benchmark/input/dwi.nii.gz ../benchmark/svht_denoise/dwi.nii.gz
/usr/bin/time -l ./svht_denoise ../benchmark/input/dwi.nii.gz ../benchmark/svht_denoise/dwi1.nii.gz -nthreads 1
/usr/bin/time -l dwidenoise ../benchmark/input/dwi.nii.gz ../benchmark/dwidenoise/dwi.nii.gz
/usr/bin/time -l dwidenoise ../benchmark/input/dwi.nii.gz ../benchmark/dwidenoise/dwi1.nii.gz -nthreads 1
```

Here are findings for an Apple MacBook with M4 Pro CPU (14 core, 10 performance). Time is elapsed wall clock, and peak RAM is the maximum resident set size, both as reported by `/usr/bin/time -l`.

| Method       | Threads | Time (ms) | Peak RAM (MiB) |
| ------------ | ------- | --------- | ------------- |
| svht_denoise |      14 |      2080 |           152 |
| svht_denoise |       1 |     14750 |           150 |
| dwidenoise   |      14 |      1950 |           228 |
| dwidenoise   |       1 |     11620 |           226 |

## Links

 - [Gavish and Donoho](https://ieeexplore.ieee.org/document/6846297) describe the singular value hard thresholding method. The pre-print [arxiv.1305.5870v3](https://arxiv.org/abs/1305.5870) is included with this repository.
 - [dwidenoise](https://mrtrix.readthedocs.io/en/dev/reference/commands/dwidenoise.html) is a popular and robust open source tool for de-noising diffusion images. It is based on Veraart et al. (2016) works [PMID:27523449](https://pubmed.ncbi.nlm.nih.gov/27523449/) and [PMID:26599599](https://pubmed.ncbi.nlm.nih.gov/26599599/). Note that it is restricted to non-commercial research usage.
 - [dwidenoise2](https://github.com/Lestropie/dwidenoise2) uses a novel method for de-noising, albeit it is also restricted to non-commercial usage.