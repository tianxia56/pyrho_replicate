#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status
set -o pipefail # Causes a pipeline to return the exit status of the last command that failed

# --- Configuration ---
RAW_VCF_DIR="/home/tx56/palmer_scratch/100kga/phased"
MASTER_INFO_FILE="/home/tx56/ycga_work/100kga/ind.info.txt"
RAW_VCF_PATTERN="Pilot.adm.all.ch{chr}.phased.ft.vcf.gz"
TEMP_DIR="./temp_sample_composition_check"
# =================================================

# --- Argument Check ---
if [ -z "$1" ]; then
    echo "Usage: $0 <chromosome_number>"
    echo "Example: $0 22"
    exit 1
fi
CHR_NUM_TO_CHECK="$1"

echo "--- VCF Sample Composition Check for Chromosome ${CHR_NUM_TO_CHECK} ---"
echo "Source Raw VCF Directory: $RAW_VCF_DIR"
echo "Master Sample Info File: $MASTER_INFO_FILE"
echo "--------------------------------------------------"

# --- Validate Input Paths ---
if [ ! -d "$RAW_VCF_DIR" ]; then echo "Error: Raw VCF directory not found: $RAW_VCF_DIR" >&2; exit 1; fi
if [ ! -f "$MASTER_INFO_FILE" ]; then echo "Error: Master sample info file not found: $MASTER_INFO_FILE" >&2; exit 1; fi

# --- Create and clean temporary directory ---
mkdir -p "$TEMP_DIR"
trap 'rm -rf "$TEMP_DIR"' EXIT

# --- Construct path to the specific chromosome VCF ---
raw_vcf_filename_to_check=$(echo "$RAW_VCF_PATTERN" | sed "s/{chr}/${CHR_NUM_TO_CHECK}/")
raw_vcf_gz_path_to_check="${RAW_VCF_DIR}/${raw_vcf_filename_to_check}"

if [ ! -f "$raw_vcf_gz_path_to_check" ]; then
    echo "Error: VCF file for chromosome ${CHR_NUM_TO_CHECK} not found at: ${raw_vcf_gz_path_to_check}" >&2
    exit 1
fi
echo "Analyzing VCF: ${raw_vcf_gz_path_to_check}"

# --- Step 1: Extract sample IDs from the VCF header ---
vcf_samples_list_temp="${TEMP_DIR}/vcf_chr${CHR_NUM_TO_CHECK}_samples.txt"
echo "Extracting sample IDs from VCF header..."
zcat "$raw_vcf_gz_path_to_check" | grep "^#CHROM" | head -n 1 | cut -f 10- | tr '\t' '\n' > "$vcf_samples_list_temp"

if [ ! -s "$vcf_samples_list_temp" ]; then
    echo "Error: No sample IDs extracted from VCF header or VCF malformed: $raw_vcf_gz_path_to_check" >&2
    exit 1
fi
vcf_sample_count=$(wc -l < "$vcf_samples_list_temp")
echo "VCF for Chr ${CHR_NUM_TO_CHECK} contains ${vcf_sample_count} sample IDs (columns)."

# --- Step 2: Process VCF samples against the master info file using awk ---
echo "Processing VCF samples and looking up metadata..."
echo "--------------------------------------------------"
echo "Counts of unique (POP | SUPOP | Countryoforigin) combinations found in VCF Chr ${CHR_NUM_TO_CHECK}:"

# Awk script (same as before)
awk_script='
    BEGIN { FS="\t"; OFS="\t"; separator="|"; header_skipped_info=0; }
    NR==FNR { 
        if (NR == 1 && $1 == "IID") { 
            header_skipped_info=1;
            next;
        }
        iid = $1; pop = $2; supop = $3; country = $4;
        if (pop == "") pop = "NA"; if (supop == "") supop = "NA"; if (country == "") country = "NA";
        info[iid] = pop separator supop separator country;
        next;
    }
    {
        vcf_iid = $1;
        if (vcf_iid in info) {
            composition_key = info[vcf_iid];
            counts[composition_key]++;
        }
    }
    END {
        # Just print data lines: Count POP SUPOP Country
        for (key in counts) {
            split(key, parts, separator);
            print counts[key] OFS parts[1] OFS parts[2] OFS parts[3];
        }
    }
'

# Execute awk, then sort, then add header, then format with column
(
  echo -e "Count\tPOP\tSUPOP\tCountryoforigin" # Print header first
  echo -e "----\t---\t-----\t-----------------" # Print separator line
  # awk script outputs data, then it's piped to sort
  awk "$awk_script" "$MASTER_INFO_FILE" "$vcf_samples_list_temp" | \
    sort -t $'\t' -k3,3 -k4,4 # Sort by SUPOP (field 3), then Countryoforigin (field 4)
) | column -t -s $'\t'


echo "--------------------------------------------------"
echo "Composition check finished."
echo "Temporary files are in: $TEMP_DIR"
# rm -rf "$TEMP_DIR" # Uncomment for automatic cleanup
