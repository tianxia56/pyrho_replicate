#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status
set -o pipefail # Causes a pipeline to return the exit status of the last command that failed

# --- Configuration (Copied and adapted from your extraction workflow script) ---
# Input VCF directory (source of chromosome-wide VCFs)
RAW_VCF_DIR="/home/tx56/palmer_scratch/100kga/phased"

# Sample List directory (Absolute Path)
SAMPLE_LIST_DIR="/home/tx56/palmer_scratch/100kga/100kga/sample_lists"

# Input VCF filename pattern (Placeholder {chr} for chromosome number)
RAW_VCF_PATTERN="Pilot.adm.all.ch{chr}.phased.ft.vcf.gz"

# Sample List filename pattern (Placeholder {pop} for population)
SAMPLE_LIST_PATTERN="{pop}.IID.list" # Assumes files like IND.IID.list, KOR.IID.list

# Temporary file directory (can be adjusted)
TEMP_DIR="./temp_sample_check"
# =================================================

# --- Argument Check ---
if [ -z "$1" ]; then
    echo "Usage: $0 <chromosome_number>"
    echo "Example: $0 22"
    exit 1
fi
CHR_NUM_TO_CHECK="$1"

echo "--- Population Sample Count Check for Chromosome ${CHR_NUM_TO_CHECK} ---"
echo "Source Raw VCF Directory: $RAW_VCF_DIR"
echo "Sample List Directory: $SAMPLE_LIST_DIR"
echo "--------------------------------------------------"

# --- Validate Input Paths ---
if [ ! -d "$RAW_VCF_DIR" ]; then echo "Error: Raw VCF directory not found: $RAW_VCF_DIR" >&2; exit 1; fi
if [ ! -d "$SAMPLE_LIST_DIR" ]; then echo "Error: Sample list directory not found: $SAMPLE_LIST_DIR" >&2; exit 1; fi

# --- Create and clean temporary directory ---
mkdir -p "$TEMP_DIR"
# trap 'rm -rf "$TEMP_DIR"' EXIT # Cleanup temp dir on script exit

# --- Construct path to the specific chromosome VCF ---
raw_vcf_filename_to_check=$(echo "$RAW_VCF_PATTERN" | sed "s/{chr}/${CHR_NUM_TO_CHECK}/")
raw_vcf_gz_path_to_check="${RAW_VCF_DIR}/${raw_vcf_filename_to_check}"

if [ ! -f "$raw_vcf_gz_path_to_check" ]; then
    echo "Error: VCF file for chromosome ${CHR_NUM_TO_CHECK} not found at: ${raw_vcf_gz_path_to_check}" >&2
    exit 1
fi
echo "Analyzing VCF: ${raw_vcf_gz_path_to_check}"

# --- Step 1: Extract sample IDs from the VCF header and sort them ---
vcf_samples_sorted_temp="${TEMP_DIR}/vcf_chr${CHR_NUM_TO_CHECK}_samples_sorted.txt"
echo "Extracting and sorting sample IDs from VCF header..."
zcat "$raw_vcf_gz_path_to_check" | grep "^#CHROM" | head -n 1 | cut -f 10- | tr '\t' '\n' | sort -u > "$vcf_samples_sorted_temp"

if [ ! -s "$vcf_samples_sorted_temp" ]; then
    echo "Error: No sample IDs extracted from VCF header, or VCF is malformed: $raw_vcf_gz_path_to_check" >&2
    exit 1
fi
vcf_sample_count=$(wc -l < "$vcf_samples_sorted_temp")
echo "VCF for Chr ${CHR_NUM_TO_CHECK} contains ${vcf_sample_count} unique sample IDs."

# --- Step 2: Discover Populations from Sample Lists ---
echo "Discovering populations from sample lists in $SAMPLE_LIST_DIR..."
populations=()
while IFS= read -r -d $'\0' file; do
    pop_name=$(basename "$file")
    pop_name=${pop_name%.IID.list}
    if [[ -n "$pop_name" && "$pop_name" != "*" ]]; then
        populations+=("$pop_name")
    fi
done < <(find "$SAMPLE_LIST_DIR" -maxdepth 1 -name '*.IID.list' -type f -print0)

IFS=$'\n' populations=($(sort <<<"${populations[*]}"))
unset IFS

if [ ${#populations[@]} -eq 0 ]; then
    echo "Error: No population sample lists (*.IID.list) found in $SAMPLE_LIST_DIR." >&2
    exit 1
fi
echo "Found sample lists for populations: ${populations[*]}"
echo "--------------------------------------------------"
echo "Matching VCF samples against population lists..."
echo "--------------------------------------------------"

# --- Step 3: For each population, find common samples ---
for pop in "${populations[@]}"; do
    pop_sample_list_filename_to_check=$(echo "$SAMPLE_LIST_PATTERN" | sed "s/{pop}/$pop/")
    pop_sample_list_path_to_check="${SAMPLE_LIST_DIR}/${pop_sample_list_filename_to_check}"

    if [ ! -f "$pop_sample_list_path_to_check" ]; then
        echo "  Population ${pop}: Sample list file not found at ${pop_sample_list_path_to_check}. Skipping." >&2
        continue
    fi

    # Sort the population sample list for comm command
    pop_samples_sorted_temp="${TEMP_DIR}/pop_${pop}_samples_sorted.txt"
    sort -u "${pop_sample_list_path_to_check}" > "$pop_samples_sorted_temp"

    pop_list_sample_count=$(wc -l < "$pop_samples_sorted_temp")
    if [ ! -s "$pop_samples_sorted_temp" ]; then
        echo "  Population ${pop}: Sample list ${pop_sample_list_path_to_check} is empty. Skipping." >&2
        continue
    fi

    # Find common samples using comm
    common_samples_temp="${TEMP_DIR}/common_${pop}_chr${CHR_NUM_TO_CHECK}.txt"
    comm -12 "$vcf_samples_sorted_temp" "$pop_samples_sorted_temp" > "$common_samples_temp"

    common_sample_count=$(wc -l < "$common_samples_temp")

    echo "  Population ${pop}:"
    echo "    - Sample list (${pop_sample_list_filename_to_check}) contains: ${pop_list_sample_count} IDs"
    echo "    - Found in VCF Chr ${CHR_NUM_TO_CHECK}: ${common_sample_count} IDs"

    # Optional: List the common samples if count is low or for debugging
    # if [ "$common_sample_count" -lt 5 ] && [ "$common_sample_count" -gt 0 ]; then
    #     echo "      Common IDs: $(paste -sd, "$common_samples_temp")"
    # fi
done

echo "--------------------------------------------------"
echo "Sample check finished."
echo "Temporary files are in: $TEMP_DIR"
# Consider adding: rm -rf "$TEMP_DIR" # if you want automatic cleanup
# Or instruct user to remove it manually.
