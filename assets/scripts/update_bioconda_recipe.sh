#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  update_bioconda_recipe.sh VERSION [options]

Updates recipes/somatem/meta.yaml for a tagged GitHub release and fills in the
sha256 from the GitHub-generated source tarball.

Options:
  --recipe PATH       Recipe meta.yaml to update (default: recipes/somatem/meta.yaml)
  --repo OWNER/REPO   GitHub repository (default: treangenlab/somatem)
  --tag-prefix TEXT   Tag prefix (default: v)
  -h, --help          Show this help message

Example:
  assets/scripts/update_bioconda_recipe.sh 0.8.0
USAGE
}

if [[ $# -eq 0 ]]; then
    usage >&2
    exit 1
fi

case "$1" in
    -h|--help)
        usage
        exit 0
        ;;
esac

version="$1"
shift

recipe="recipes/somatem/meta.yaml"
repo="treangenlab/somatem"
tag_prefix="v"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --recipe)
            recipe="$2"
            shift 2
            ;;
        --repo)
            repo="$2"
            shift 2
            ;;
        --tag-prefix)
            tag_prefix="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ ! "${version}" =~ ^[0-9]+([.][0-9]+)*([._-]?(a|b|rc|dev|post)[0-9]*)?$ ]]; then
    echo "Version '${version}' does not look like a conda package version." >&2
    exit 1
fi

if [[ ! -f "${recipe}" ]]; then
    echo "Recipe not found: ${recipe}" >&2
    exit 1
fi

for exe in awk curl perl sha256sum; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
        echo "Required command not found: ${exe}" >&2
        exit 1
    fi
done

tag="${tag_prefix}${version}"
url="https://github.com/${repo}/archive/refs/tags/${tag}.tar.gz"
recipe_url="https://github.com/${repo}/archive/refs/tags/${tag_prefix}{{ version }}.tar.gz"
tmpdir="$(mktemp -d)"
archive="${tmpdir}/somatem-${tag}.tar.gz"
trap 'rm -rf "${tmpdir}"' EXIT

echo "Downloading ${url}"
curl -fsSL "${url}" -o "${archive}"
sha256="$(sha256sum "${archive}" | awk '{print $1}')"

perl -0pi -e 's/{% set version = "[^"]+" %}/{% set version = "'"${version}"'" %}/' "${recipe}"
perl -0pi -e 's#url: https://github\.com/[^/]+/[^/]+/archive/refs/tags/[^\s]+\.tar\.gz#url: '"${recipe_url}"'#' "${recipe}"
perl -0pi -e 's/sha256: [A-Za-z0-9_]+/sha256: '"${sha256}"'/' "${recipe}"

echo "Updated ${recipe}"
echo "  version: ${version}"
echo "  tag:     ${tag}"
echo "  sha256:  ${sha256}"
