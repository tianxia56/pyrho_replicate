#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status
set -o pipefail # Causes a pipeline to return the exit status of the last command that exited with a non-zero status

# --- Configuration ---
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Input VCF directory (source of chromosome-wide VCFs)
RAW_VCF_DIR="/home/tx56/palmer_scratch/100kga/phased" # Location of chr*.vcf.gz files

# Sample List directory (Absolute Path)
SAMPLE_LIST_DIR="/home/tx56/palmer_scratch/100kga/100kga/sample_lists" # Location of POP.IID.list files

# Output base directory for the new per-population VCFs (Absolute Path)
OUTPUT_RMAP_DIR="/home/tx56/palmer_scratch/ga100k_rmap"

# Input VCF filename pattern (Placeholder {chr} for chromosome number)
RAW_VCF_PATTERN="Pilot.adm.all.ch{chr}.phased.ft.vcf.gz"

# Sample List filename pattern (Placeholder {pop} for population)
SAMPLE_LIST_PATTERN="{pop}.IID.list" # Assumes files like IND.IID.list, KOR.IID.list

# --- SLURM Configuration (can be adjusted) ---
# Log directory for this task (Absolute Path)
LOG_DIR_BASE="/home/tx56/palmer_scratch/ga100k_rmap"
SBATCH_SCRIPT_DIR="${LOG_DIR_BASE}/sbatch_scripts"
SBATCH_PARTITION="week"
SBATCH_TIME="2-00:00:00" # Adjusted time, VCFtools extraction can be quick or slow depending on VCF size and sample count
SBATCH_NODES="1"
SBATCH_NTASKS="1"
SBATCH_CPUS="1" # VCFtools is mostly single-threaded for this operation
SBATCH_MEM_PER_CPU="8000" # Per CPU (in MB), adjust if needed. 8GB should be plenty for most vcftools recode.

# --- Tools (Ensure they are in PATH within SLURM jobs) ---
VCFTOOLS_PATH="vcftools"
BGZIP_PATH="bgzip"
TABIX_PATH="tabix"
# =================================================

echo "--- VCF per-Population Extraction Workflow ---"
echo "Source Raw VCF Directory: $RAW_VCF_DIR"
echo "Source Raw VCF Pattern: $RAW_VCF_PATTERN"
echo "Sample List Directory: $SAMPLE_LIST_DIR"
echo "Sample List Pattern: $SAMPLE_LIST_PATTERN"
echo "Output Base Directory: $OUTPUT_RMAP_DIR"
echo "Log Directory: $LOG_DIR_BASE"
echo "--------------------------------------------------"

# --- Validate Input Paths ---
if [ ! -d "$RAW_VCF_DIR" ]; then echo "Error: Raw VCF directory not found: $RAW_VCF_DIR" >&2; exit 1; fi
if [ ! -d "$SAMPLE_LIST_DIR" ]; then echo "Error: Sample list directory not found: $SAMPLE_LIST_DIR" >&2; exit 1; fi

# --- Create necessary output and log directories ---
echo "Ensuring all output and log directories exist..."
mkdir -p "$OUTPUT_RMAP_DIR"
mkdir -p "$LOG_DIR_BASE"
mkdir -p "$SBATCH_SCRIPT_DIR"
echo "Output and log directories checked/created."

# --- Discover Populations from Sample Lists ---
echo "Discovering populations from sample lists in $SAMPLE_LIST_DIR..."
populations=()
while IFS= read -r -d $'\0' file; do
    pop_name=$(basename "$file")
    pop_name=${pop_name%.IID.list} # Use parameter expansion for suffix removal
    if [[ -n "$pop_name" && "$pop_name" != "*" ]]; then
        populations+=("$pop_name")
    fi
done < <(find "$SAMPLE_LIST_DIR" -maxdepth 1 -name '*.IID.list' -type f -print0)

# Sort the discovered populations
IFS=$'\n' populations=($(sort <<<"${populations[*]}"))
unset IFS

if [ ${#populations[@]} -eq 0 ]; then
    echo "Error: No population sample lists (*.IID.list) found in $SAMPLE_LIST_DIR." >&2
    exit 1
fi
echo "Found sample lists for populations: ${populations[*]}"
echo "IMPORTANT: Ensure these sample lists (e.g., IND.IID.list) contain one sample ID per line, with NO duplicates or _ID suffixes."

# --- Generate and Submit one SLURM job per Population ---
echo "Generating and submitting SLURM jobs for each population..."
job_ids=() # Array to store submitted job IDs

for pop in "${populations[@]}"; do
    echo "  Generating sbatch script for population ${pop}..."

    # Construct path to the sample list for this population
    pop_sample_list_filename=$(echo "$SAMPLE_LIST_PATTERN" | sed "s/{pop}/$pop/")
    pop_sample_list_path="${SAMPLE_LIST_DIR}/${pop_sample_list_filename}"

    if [ ! -f "$pop_sample_list_path" ]; then
        echo "  Warning: Sample list file not found for population ${pop}: ${pop_sample_list_path}. Skipping." >&2
        continue
    fi

    # Define SLURM job name and log paths
    job_name="${pop}_extract_vcf"
    sbatch_script_path="${SBATCH_SCRIPT_DIR}/submit_${job_name}.sbatch"
    slurm_output_log="${LOG_DIR_BASE}/${job_name}.%j.out"
    slurm_error_log="${LOG_DIR_BASE}/${job_name}.%j.err"

    # Per-population output directory (will be created by the sbatch script)
    pop_output_dir="${OUTPUT_RMAP_DIR}/${pop}"

    # --- Generate the sbatch script ---
    cat << EOF > "$sbatch_script_path"
#!/bin/bash
# --- SLURM Directives ---
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --time=${SBATCH_TIME}
#SBATCH --nodes=${SBATCH_NODES}
#SBATCH --ntasks=${SBATCH_NTASKS}
#SBATCH --cpus-per-task=${SBATCH_CPUS}
#SBATCH --mem-per-cpu=${SBATCH_MEM_PER_CPU}
#SBATCH --job-name=${job_name}
#SBATCH --output=${slurm_output_log}
#SBATCH --error=${slurm_error_log}

# --- Job Script Content ---
set -e # Exit immediately if a command exits with a non-zero status.
set -o pipefail # Causes a pipeline to return the exit status of the last command that exited with a non-zero status

echo "========================================================"
echo "SLURM Job ID: \$SLURM_JOB_ID"
echo "Job Name: \$SLURM_JOB_NAME"
echo "Running on host: \$(hostname)"
echo "Processing Population ID: ${pop}"
echo "Sample List: ${pop_sample_list_path}"
echo "Raw VCF Dir (source): ${RAW_VCF_DIR}"
echo "Raw VCF Pattern (source): ${RAW_VCF_PATTERN}"
echo "Output Directory for this Pop: ${pop_output_dir}"
echo "========================================================"

# --- Environment Setup & Tool Check ---
echo "Checking for required tools..."
command -v "${VCFTOOLS_PATH}" >/dev/null 2>&1 || { echo "Error: VCFtools command ('${VCFTOOLS_PATH}') not found." >&2; exit 1; }
command -v "${BGZIP_PATH}" >/dev/null 2>&1 || { echo "Error: bgzip command ('${BGZIP_PATH}') not found." >&2; exit 1; }
command -v "${TABIX_PATH}" >/dev/null 2>&1 || { echo "Error: tabix command ('${TABIX_PATH}') not found." >&2; exit 1; }
echo "All tools found."

# --- Create Population-Specific Output Directory ---
echo "Creating output directory: ${pop_output_dir}"
mkdir -p "${pop_output_dir}"

# --- Loop through Chromosomes 1-22 ---
job_failed_chr_count=0
for chr_num in {1..22}; do
    echo -e "\\n--- Processing Chromosome \${chr_num} for Population ${pop} ---"

    # Construct raw VCF filename from pattern
    raw_vcf_filename=\$(echo "${RAW_VCF_PATTERN}" | sed "s/{chr}/\${chr_num}/")
    raw_vcf_gz_path="${RAW_VCF_DIR}/\${raw_vcf_filename}"

    # Define output VCF.gz path for this pop & chr
    output_vcf_gz_path="${pop_output_dir}/${pop}.chr\${chr_num}.vcf.gz"
    vcftools_log_prefix="${pop_output_dir}/${pop}.chr\${chr_num}.vcftools_extract" # For VCFtools log

    echo "  Source VCF: \${raw_vcf_gz_path}"
    echo "  Output VCF: \${output_vcf_gz_path}"
    echo "  Sample List: ${pop_sample_list_path}"

    if [ ! -f "\${raw_vcf_gz_path}" ]; then
        echo "  Warning: Source VCF not found: \${raw_vcf_gz_path}. Skipping chromosome \${chr_num}." >&2
        ((job_failed_chr_count++))
        continue
    fi

    # --- Step 1: Extract samples using VCFtools and pipe to bgzip ---
    echo "  Extracting samples with VCFtools and compressing with bgzip..."
    # Using --stdout with vcftools and piping to bgzip is generally preferred for VCFs
    # It ensures block-gzipped output suitable for tabix.
    # Using a temporary log file for vcftools output within the loop
    if ! ("${VCFTOOLS_PATH}" --gzvcf "\${raw_vcf_gz_path}" \\
        --keep "${pop_sample_list_path}" \\
        --recode --recode-INFO-all \\
        --stdout 2> "\${vcftools_log_prefix}.log" | \\
        "${BGZIP_PATH}" -c > "\${output_vcf_gz_path}"); then
        
        echo "  Error: VCFtools or bgzip failed for Chr \${chr_num}." >&2
        echo "  Check VCFtools log: \${vcftools_log_prefix}.log"
        # Check specifically for keep errors if possible (VCFtools logs might vary)
        if grep -q -e "Error: Could not find allele frequency" -e "Keeping 0 individuals" -e "No data left" "\${vcftools_log_prefix}.log" 2>/dev/null ; then
             echo "  VCFtools Error Detail: Problems with sample IDs in '--keep' list (likely did not match VCF Sample IDs) or no variants left after selection." >&2
             echo "  DEBUG: First few sample IDs from VCF header:" >&2
             (gunzip -c "\${raw_vcf_gz_path}" 2>/dev/null || cat "\${raw_vcf_gz_path}") | grep -m1 '^#CHROM' | tr '\t' '\n' | tail -n 5 || echo "  DEBUG: Failed to read VCF header." >&2
             echo "  DEBUG: First few sample IDs from keep list:" >&2
             head -n 5 "${pop_sample_list_path}" || echo "  DEBUG: Failed to read keep list." >&2
        fi
        rm -f "\${output_vcf_gz_path}" # Clean up potentially incomplete output
        ((job_failed_chr_count++))
        continue # Move to next chromosome
    else
        echo "  VCFtools and bgzip successful for Chr \${chr_num}."
        # Clean up VCFtools log on success
        rm -f "\${vcftools_log_prefix}.log"
    fi

    # --- Step 2: Index the output VCF.gz file ---
    echo "  Indexing output VCF with tabix..."
    if ! "${TABIX_PATH}" -p vcf "\${output_vcf_gz_path}"; then
        echo "  Error: tabix failed for Chr \${chr_num} on \${output_vcf_gz_path}." >&2
        rm -f "\${output_vcf_gz_path}.tbi" # Clean up potentially incomplete index
        ((job_failed_chr_count++))
        continue # Move to next chromosome
    fi
    echo "  Tabix indexing successful for Chr \${chr_num}."

done # --- End of chromosome loop ---

# --- Final Job Status ---
echo "========================================================"
if [ \${job_failed_chr_count} -eq 0 ]; then
    echo "SLURM Job for Population ${pop} finished successfully."
    exit 0
else
    echo "SLURM Job for Population ${pop} finished with \${job_failed_chr_count} chromosome(s) failing." >&2
    exit 1
fi
echo "========================================================"

EOF
    # --- End of Here Document ---

    chmod +x "$sbatch_script_path"
    echo "    Generated sbatch script: $sbatch_script_path"
    sbatch_output=$(sbatch "$sbatch_script_path")
    job_id=$(echo "$sbatch_output" | awk '{print $NF}')
    if [[ "$job_id" =~ ^[0-9]+$ ]]; then
        echo "    Submitted job ID for ${pop}: ${job_id}"
        job_ids+=("${job_id}")
    else
        echo "    Error submitting job for population ${pop}: ${sbatch_output}" >&2
    fi
done # --- End of population loop ---

echo "--------------------------------------------------"
if [ ${#job_ids[@]} -gt 0 ]; then
    echo "All ${#job_ids[@]} population extraction jobs submitted."
    echo "Monitor jobs using: squeue -u $USER"
    echo "Check SLURM output/error logs in: $LOG_DIR_BASE (e.g., POP_extract_vcf.%j.out)"
    echo "Final VCFs will be in: ${OUTPUT_RMAP_DIR}/{pop}/{pop}.chr{chr}.vcf.gz"
else
    echo "No jobs were submitted. Check for errors in population discovery or sample list availability."
fi
echo "--------------------------------------------------"

exit 0
