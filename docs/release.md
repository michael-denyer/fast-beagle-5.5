# Releases

A release of this repository is a Beagle 5.5 edition of fast-beagle, versioned 5.5.x. [fast-beagle-5.4](https://github.com/michael-denyer/fast-beagle-5.4) releases the 5.4.x versions of the same `fast-beagle` package. Each release publishes a GitHub Release with:

- `fast-beagle-<version>.tar.gz`, the source tarball
- `fast-beagle-<version>-linux-x86_64.tar.gz`, `fast-beagle-<version>-linux-arm64.tar.gz` and `fast-beagle-<version>-macos-arm64.tar.gz`, each holding a `fast-beagle` binary with htslib linked statically, and the licenses
- `fast-beagle-<version>-macos-arm64.pkg`, a macOS installer that puts the macOS binary in `/usr/local/bin` and the licenses in `/usr/local/share/doc/fast-beagle`. The binary and the installer are signed with Developer ID certificates and Apple has notarised the installer.
- `SHA256SUMS`, the SHA-256 of every file above

The bioconda package `fast-beagle` builds from the source tarball with the recipe in `release/conda/`.

## Cut a release

1. Merge every change for the release into `main`, with the gate passing.
2. As a repository administrator, tag the `main` commit `v5.5.<n>` and push the protected tag. The release workflow fails on a malformed or unprotected release tag, on a commit that is not on `main`, and on a commit without a successful `gate` check run.

   ```bash
   git tag v5.5.1 origin/main
   git push origin v5.5.1
   ```

3. Wait for the release workflow run on the tag. It publishes the GitHub Release once every build and check has passed. If a job fails, nothing is published. Rerun a transient failure on the original tag. If the source needs a fix, merge it to `main` and use a new version tag.
4. Submit the recipe to [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) as `recipes/fast-beagle/`. The run's `conda-recipe` artifact holds `recipe.yaml` with the release's version and the tarball's sha256 filled in, with `build.sh` and `small.vcf`. After the first submission, bioconda's autobump bot opens the pull request for each new release.

## What the release workflow proves

`.github/workflows/release.yml` runs on a pushed `v5.5.*` tag, on pull requests that change the workflow, `release/` or the `Makefile`, and by hand. Pull requests, including forks, and manual runs build an unsigned installer and upload it as `preview-pkg-macos-arm64`. Branch and PR builds use the version in `release/conda/recipe.yaml`.

- On a tag or a manual run, the `source` job first checks the commit against the GitHub API. The commit must be on `main`, so the compare of `main` with the commit reports `identical` or `behind`, and it must have a `gate` check run that succeeded. Otherwise the job fails and nothing builds. A manual run takes the check so that it can be tried without pushing a tag.
- The `source` job makes the tarball with `git archive`, with every path under `fast-beagle-<version>/`. It fills in the recipe's version and sha256 from the tarball. Every other job builds from this tarball.
- The `conda` jobs build the recipe with rattler-build on `linux-64`, `linux-aarch64` and `osx-arm64`, from conda-forge and bioconda, as bioconda does. The build uses conda's compilers and its `CFLAGS` (`-O3` on `linux-aarch64`), and the Makefile appends its fixed floating-point flags after them. The Makefile's `-Werror` applies, so a warning from conda's compilers fails the build. rattler-build runs the recipe's test, a phasing run on `small.vcf`. The job then installs the package into a new environment and runs `tests/check-oracle.sh` on the installed `fast-beagle`, so every oracle case must match its hash at 1, 2 and 18 threads.
- The `binaries` jobs run `release/build-static.sh`. It builds htslib 1.24, the version that bioconda packages, from its checksum-verified release tarball, without libcurl, GCS and S3 support, plugins, bzip2 and lzma. It builds htslib against libdeflate 1.25, also from its checksum-verified release tarball, so that BGZF compression and decompression use libdeflate rather than zlib. It links htslib and libdeflate statically and fails if the binary loads a shared htslib or libdeflate. The package holds libdeflate's license as `LICENSE.libdeflate`. The Linux binaries build in `manylinux_2_28` containers pinned by digest, so they need glibc 2.28 or later, and run on the runner itself for the check. Each job runs `tests/check-oracle.sh` on the binary unpacked from its package.
- On a protected release-tag push, the `macos-pkg` job runs `release/build-pkg.sh` on the unpacked `macos-arm64` package, so the installer holds the binary that the `binaries` job checked. The script signs the binary with the Developer ID Application certificate under the hardened runtime, builds the installer with `pkgbuild`, signs it with the Developer ID Installer certificate, submits it to Apple's notary service and staples the ticket. It fails unless Apple accepts the submission and `spctl` accepts the stapled installer. The job then installs the `.pkg` on the runner and runs the recipe's phasing test with `/usr/local/bin/fast-beagle`. The preview build uses the same package layout and install test.
- The `publish` job runs only on a pushed tag, after every other job has passed. It writes `SHA256SUMS`, attests every tarball, the installer and `SHA256SUMS` with `actions/attest-build-provenance`, and creates the GitHub Release with the default `GITHUB_TOKEN`, which has `contents: write`, `id-token: write` and `attestations: write` in this job only. A failed attestation stops the job before the release is created. Anyone can check a downloaded file's provenance with `gh attestation verify <file> --repo michael-denyer/fast-beagle-5.5`.

rattler-build 0.76.1, micromamba 2.9.0-0, htslib 1.24 and libdeflate 1.25 are pinned by version and SHA-256 in the workflow and in `release/build-static.sh`.

## Signing secrets

Configure the `apple-signing` GitHub environment with **Selected branches and tags** and a single **Tag** rule matching `v5.5.*`. Keep an active tag ruleset on `refs/tags/v5.5.*` restricting creation, updates and deletion to release administrators. The signing job requires both a tag push and GitHub's protected-ref flag. Manual dispatches, including those targeting a tag, build previews.

Store the five secrets below only in `apple-signing`. To migrate, provision the environment from the original credentials, verify all five names, then remove their repository-scoped copies. Repository secrets remain accessible to same-repository workflows while those copies exist.

| Secret | Value |
| --- | --- |
| `MACOS_CERT_P12` | The Developer ID Application and Developer ID Installer identities, exported together from Keychain Access as one `.p12` and base64-encoded |
| `MACOS_CERT_PASSWORD` | The export password of that `.p12` |
| `NOTARY_KEY_P8` | An App Store Connect API team key with the Developer role, the contents of `AuthKey_<key id>.p8` |
| `NOTARY_KEY_ID` | The key's ID |
| `NOTARY_ISSUER_ID` | The issuer ID shown above the team keys |

Developer ID certificates expire after five years. Replace `MACOS_CERT_P12` and `MACOS_CERT_PASSWORD` when they are renewed.
