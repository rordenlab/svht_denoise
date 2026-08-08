svht_denoise for macOS
======================

This installer places a single executable in /usr/local/bin, which is on the
default PATH, so it is ready to use immediately:

  svht_denoise

The binary is universal (Apple Silicon and Intel). This build requires macOS 11
or newer; the minimum is set at build time by MACOSX_DEPLOYMENT_TARGET. It is
entirely self-contained: it links only against libSystem and libz, both of which
ship with macOS. There is no folder to keep together and nothing
else to download. You may copy it to any other Mac and it will run.

What it does
------------

svht_denoise reduces thermal noise in multi-volume NIfTI datasets, such as
diffusion-weighted MRI, by singular value hard thresholding over local image
patches (Gavish and Donoho, 2014).

  svht_denoise dwi.nii.gz denoised.nii.gz

With a mask and a noise-level map:

  svht_denoise dwi.nii.gz denoised.nii.gz -mask brain.nii.gz -noise sigma.nii.gz

If you also have the phase, pass it and the complex data is rotated onto the
real axis first. Noise then stays zero-mean Gaussian rather than Rician, so the
low-SNR noise floor is avoided rather than denoised:

  svht_denoise mag.nii.gz denoised.nii.gz -phase phase.nii.gz

Note that the output of a -phase run is SIGNED: background voxels scatter about
zero instead of piling up on the noise floor, which is the point, but the result
is no longer a magnitude image.

Run `svht_denoise -help` for the full option list.

Verifying this package
----------------------

The package is signed with a Developer ID Installer certificate and notarized by
Apple. The notarization ticket is stapled to the package, so it validates with no
network connection. To check for yourself:

  pkgutil --check-signature svht_denoise-<version>-macos-universal.pkg
  spctl --assess --type install --verbose=4 svht_denoise-<version>-macos-universal.pkg
  xcrun stapler validate svht_denoise-<version>-macos-universal.pkg

The installed executable is signed with a Developer ID Application certificate.
A notarization ticket cannot be stapled to a bare executable, so this checks the
signature and its issuer rather than a ticket:

  codesign --verify --strict --verbose=2 /usr/local/bin/svht_denoise
  codesign -dvv /usr/local/bin/svht_denoise 2>&1 | grep Authority

Uninstalling
------------

  sudo rm /usr/local/bin/svht_denoise

Licence
-------

Mozilla Public License 2.0. See LICENSE in the source distribution at
https://github.com/rordenlab/svht_denoise
