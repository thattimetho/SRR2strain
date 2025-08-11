#!/bin/bash

# Check supplied params
if [ $# -lt 2 ]; then
	echo "Usage: $0 <destination_dir> <original_dir> <dataset_name>"
	exit 1
fi

og_dir="${2}/*"

# Check if destination folder exists, and otherwise create dir
dest_subdir=$1
if [ ! -d "$dest_subdir" ]; then
  mkdir -p "$dest_subdir"
fi

# Get current date
  date=$(date +%Y_%m_%d)

# Set new run data output dir
dataset_name=$3
output_dir="${dest_subdir}/run_${date}_${dataset_name}"
run_config="${output_dir}/metadata/config_run_${date}_${dataset_name}.yml"

# Check if new output_dir exists, and otherwise create dir
if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

# Create necessary subfolders in new run data output dir
mkdir -p "$output_dir/instrain" "$output_dir/mappings" "$output_dir/sra_temp" "$output_dir/raw_reads" "$output_dir/metadata" "$output_dir/benchmarks" "$output_dir/logs"

# Copy template yml to configs dir if no existing yml matches new run id
if [ ! -f "$run_config" ]; then
  cp "resources/template.yml" "$run_config"
fi

# Loop over all files
for file in ${og_dir}; do
  # Check file extension
  ext="${file##*.}"

  # Move file to different dir depending on file extension
  case $ext in
    txt)
      ln -s "$(realpath "$file")" "$output_dir/metadata/$(basename "$file")"
      ;;
    # fa)
    #   ln -s "$(realpath "$file")" "$output_dir/metadata/$(basename "$file")"
    #   ;;
    # gff3)
    #   ln -s "$(realpath "$file")" "$output_dir/annotation/$(basename "$file")"
    #   ;;
    # bed4)
    #   ln -s "$(realpath "$file")" "$output_dir/annotation/$(basename "$file")"
    #   ;;
  esac

done

