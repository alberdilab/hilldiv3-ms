#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_dir="$root/sources"
mkdir -p "$source_dir"

archive="$source_dir/Nitrifiers.tar.gz"
tree="$source_dir/gtdb-shallow.tree"
archive_url="https://zenodo.org/api/records/17162544/files/Nitrifiers.tar.gz/content"
tree_url="https://raw.githubusercontent.com/cmc-aau/mfd_sr_mags/957591d4104462ec6e52ab58d85a69798057e582/analysis/datasets/gtdb-shallow.tree"
archive_md5="9b65d9984cb41a929613758919ed2fcf"
tree_sha256="693976bea56093394704861c8229a26305a6aec21d8f60fde9debedc89acc112"

if [[ ! -f "$archive" ]]; then
  curl -fL --retry 3 "$archive_url" -o "$archive"
fi
if [[ "$(md5 -q "$archive" 2>/dev/null || md5sum "$archive" | cut -d' ' -f1)" != "$archive_md5" ]]; then
  echo "Nitrifiers archive checksum mismatch" >&2
  exit 1
fi

if [[ ! -f "$tree" ]]; then
  curl -fL --retry 3 "$tree_url" -o "$tree"
fi
if [[ "$(shasum -a 256 "$tree" | cut -d' ' -f1)" != "$tree_sha256" ]]; then
  echo "Source tree checksum mismatch" >&2
  exit 1
fi

mag_metadata="$source_dir/mags_shallow_all.tsv"
abundance="$source_dir/MFD_SRnodrep_tax_relative_abundance.tsv.xz"
if [[ ! -f "$mag_metadata" ]]; then
  tar -xOzf "$archive" ./Nitrifiers/data/mags_shallow_all.tsv > "$mag_metadata"
fi
if [[ ! -f "$abundance" ]]; then
  tar -xOzf "$archive" ./Nitrifiers/data/MFD_SRnodrep_tax_relative_abundance.tsv | xz -T2 -6 > "$abundance"
fi

echo "Verified source files in $source_dir"
