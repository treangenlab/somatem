#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  prep_bioconda_dev_test.sh --prefix CONDA_PREFIX [options]

Refreshes an installed somatem conda package with the local wrapper script.
Optionally replaces the installed pipeline directory with a symlink to this
checkout for quick recipe/debug testing.

Options:
  --prefix PATH      Conda environment prefix that contains somatem
  --version VERSION  Installed somatem package version (default: recipe version)
  --repo PATH        Local somatem checkout (default: git repo root)
  --link-source      Symlink share/somatem-VERSION to the local checkout
  -h, --help         Show this help message

Example:
  assets/scripts/prep_bioconda_dev_test.sh --prefix "$CONDA_PREFIX" --link-source
USAGE
}

repo_root="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
version=""
prefix="${CONDA_PREFIX:-}"
link_source=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)
            prefix="$2"
            shift 2
            ;;
        --version)
            version="$2"
            shift 2
            ;;
        --repo)
            repo_root="$2"
            shift 2
            ;;
        --link-source)
            link_source=1
            shift
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

if [[ -z "${version}" ]]; then
    version="$(sed -n 's/{% set version = "\(.*\)" %}/\1/p' "${repo_root}/recipes/somatem/meta.yaml" 2>/dev/null || true)"
fi

if [[ -z "${prefix}" ]]; then
    echo "No conda prefix provided. Use --prefix or activate the target environment." >&2
    exit 1
fi

if [[ -z "${version}" ]]; then
    echo "Could not determine somatem version. Use --version." >&2
    exit 1
fi

bin_dir="${prefix}/bin"
pipeline_dir="${prefix}/share/somatem-${version}"
wrapper="${bin_dir}/somatem"

if [[ ! -f "${repo_root}/bin/somatem" ]]; then
    echo "Local wrapper not found: ${repo_root}/bin/somatem" >&2
    exit 1
fi

mkdir -p "${bin_dir}"
cp "${repo_root}/bin/somatem" "${wrapper}"
sed -i "s|@@PIPELINE_DIR@@|${pipeline_dir}|g" "${wrapper}"
chmod +x "${wrapper}"

if [[ "${link_source}" -eq 1 ]]; then
    mkdir -p "$(dirname "${pipeline_dir}")"
    if [[ -e "${pipeline_dir}" && ! -L "${pipeline_dir}" ]]; then
        backup="${pipeline_dir}.original_backup"
        if [[ ! -e "${backup}" ]]; then
            mv "${pipeline_dir}" "${backup}"
        else
            echo "Backup already exists: ${backup}" >&2
            exit 1
        fi
    elif [[ -L "${pipeline_dir}" ]]; then
        rm "${pipeline_dir}"
    fi
    ln -s "${repo_root}" "${pipeline_dir}"
fi

echo "Updated ${wrapper}"
echo "Pipeline directory: ${pipeline_dir}"
