# Bioconda Release Checklist

This repo keeps a mirror of the Bioconda recipe in `recipes/somatem/`. The
canonical recipe that Bioconda builds lives in `bioconda/bioconda-recipes`, but
keeping this copy current makes release prep easier to review before opening a
Bioconda PR.

## Release Conventions

- Use the lowercase upstream repository URL: `https://github.com/treangenlab/somatem`.
- Tag upstream releases as `vX.Y.Z`, for example `v0.8.0`.
- Keep the conda package version numeric, for example `0.8.0`.
- For a new upstream version, keep `build:number` at `0`.
- Only bump `build:number` for recipe-only rebuilds of the same upstream version.

## 1. Prepare the Somatem Release

Update the release version in `nextflow.config` and commit the release changes.
The public tag must exist before the recipe checksum can be finalized, because
Bioconda verifies the checksum of GitHub's generated source tarball.

```bash
git tag -a v0.8.0 -m "somatem v0.8.0"
git push origin main
git push origin v0.8.0
```

## 2. Update This Repo's Recipe Copy

After the tag is available on GitHub, update the local recipe copy and checksum:

```bash
assets/scripts/update_bioconda_recipe.sh 0.8.0
```

Review the recipe diff:

```bash
git diff recipes/somatem/meta.yaml recipes/somatem/build.sh
```

The recipe should point at:

```yaml
source:
  url: https://github.com/treangenlab/somatem/archive/refs/tags/v{{ version }}.tar.gz
```

## 3. Test the Recipe Locally

For a quick conda-build smoke test from this repo:

```bash
conda build recipes/somatem/
```

For the closer Bioconda test, work from a fork of `bioconda-recipes`:

```bash
git checkout master
git pull upstream master
git checkout -b update-somatem-0.8.0
cp /path/to/somatem/recipes/somatem/meta.yaml recipes/somatem/meta.yaml
cp /path/to/somatem/recipes/somatem/build.sh recipes/somatem/build.sh
bioconda-utils lint --packages somatem
bioconda-utils build --packages somatem
```

## 4. Open the Bioconda PR

In your `bioconda-recipes` fork:

```bash
git status
git add recipes/somatem/meta.yaml recipes/somatem/build.sh
git commit -m "Update somatem to 0.8.0"
git push origin update-somatem-0.8.0
```

Open a pull request from your fork to `bioconda:master`. In the PR description,
include:

- Upstream release tag: `https://github.com/treangenlab/somatem/releases/tag/v0.8.0`
- Local tests run, for example `bioconda-utils lint --packages somatem`
- Any intentional recipe changes, especially dependency or wrapper-script changes

After CI passes, ask for the review label if needed:

```text
@BiocondaBot please add label
```

## Local Installed-Package Debugging

To test local wrapper changes inside an existing installed somatem conda
environment:

```bash
assets/scripts/prep_bioconda_dev_test.sh --prefix "$CONDA_PREFIX"
```

To point that installed package at this checkout while debugging:

```bash
assets/scripts/prep_bioconda_dev_test.sh --prefix "$CONDA_PREFIX" --link-source
```
