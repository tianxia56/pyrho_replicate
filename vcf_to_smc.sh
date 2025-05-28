#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status
set -o pipefail # Causes a pipeline to return the exit status of the last command that exited with a non-zero status

# --- Configuration ---
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Base directory containing the per-population, per-chromosome VCFs
# (e.g., /home/tx56/palmer_scratch/ga100k_rmap/IND/IND.chr1.vcf.gz)
INPUT_VCF_BASE_DIR="/home/tx56/palmer_scratch/ga100k_rmap"

# Base output directory for SMC++ results
SMC_OUTPUT_BASE_DIR="/home/tx56/palmer_scratch/ga100k_smcpp_out" # CHANGED: Specific output for smc++

# Log directory for this SMC++ workflow
LOG_DIR_SMC="${SMC_OUTPUT_BASE_DIR}/logs" # Place logs within the SMC output directory
SBATCH_SCRIPT_DIR_SMC="${LOG_DIR_SMC}/sbatch_scripts"

# --- SMC++ Specific Parameters ---
SMC_PLUS_PLUS_PATH="smc++" # Ensure smc++ is in PATH or provide full path
BCFTOOLS_PATH="bcftools"   # Ensure bcftools is in PATH or provide full path
SMC_CUTOFF="1000"          # The -c parameter for vcf2smc

# --- SLURM Configuration (can be adjusted from your example) ---
SBATCH_PARTITION="day"     # From your example
SBATCH_TIME="23:00:00"     # From your example
SBATCH_NODES="1"
SBATCH_NTASKS="1"
SBATCH_CPUS_PER_TASK="1"   # smc++ vcf2smc is generally single-threaded
SBATCH_MEM_PER_CPU="50000" # From your example (50GB per CPU) - this is quite high, adjust if needed

# =================================================

echo "--- SMC++ vcf2smc Workflow ---"
echo "Input VCF Base Directory: $INPUT_VCF_BASE_DIR"
echo "SMC Output Base Directory: $SMC_OUTPUT_BASE_DIR"
echo "Log Directory: $LOG_DIR_SMC"
echo "SMC++ Cutoff (-c): $SMC_CUTOFF"
echo "--------------------------------------------------"

# --- Validate Input Path ---
if [ ! -d "$INPUT_VCF_BASE_DIR" ]; then
    echo "Error: Input VCF base directory not found: $INPUT_VCF_BASE_DIR" >&2
    exit 1
fi

# --- Create necessary output and log directories ---
echo "Ensuring all output and log directories exist..."
mkdir -p "$SMC_OUTPUT_BASE_DIR"
mkdir -p "$LOG_DIR_SMC"
mkdir -p "$SBATCH_SCRIPT_DIR_SMC"
echo "Output and log directories checked/created."

# --- Discover Populations from subdirectories in INPUT_VCF_BASE_DIR ---
echo "Discovering populations from subdirectories in $INPUT_VCF_BASE_DIR..."
populations=()
# Find directories (maxdepth 1) that are not hidden and are actual directories
while IFS= read -r -d $'\0' dir_path; do
    pop_name=$(basename "$dir_path")
    # Basic check to avoid hidden dirs or files, ensure it's a directory
    if [[ -d "$dir_path" && ! "$pop_name" =~ ^\..* ]]; then
        populations+=("$pop_name")
    fi
done < <(find "$INPUT_VCF_BASE_DIR" -maxdepth 1 -mindepth 1 -type d -print0)


if [ ${#populations[@]} -eq 0 ]; then
    echo "Error: No population subdirectories found in $INPUT_VCF_BASE_DIR." >&2
    echo "Expected structure: $INPUT_VCF_BASE_DIR/POP_NAME/POP_NAME.chrX.vcf.gz" >&2
    exit 1
fi

# Sort the discovered populations
IFS=$'\n' populations=($(sort <<<"${populations[*]}"))
unset IFS

echo "Found populations: ${populations[*]}"

# --- Generate and Submit one SLURM job per Population ---
echo "Generating and submitting SLURM jobs for each population..."
job_ids=() # Array to store submitted job IDs

for pop in "${populations[@]}"; do
    echo "  Processing population: ${pop}"

    # Define population-specific input VCF directory
    pop_input_vcf_dir="${INPUT_VCF_BASE_DIR}/${pop}"
    # Define population-specific output directory for SMC files
    pop_smc_output_dir="${SMC_OUTPUT_BASE_DIR}/${pop}" # e.g., .../ga100k_smcpp_out/IND

    # Check if the population VCF directory exists
    if [ ! -d "$pop_input_vcf_dir" ]; then
        echo "  Warning: VCF directory for population ${pop} not found: ${pop_input_vcf_dir}. Skipping." >&2
        continue
    fi

    # Define SLURM job name and log paths
    job_name="${pop}_smcpp_vcf2smc"
    sbatch_script_path="${SBATCH_SCRIPT_DIR_SMC}/submit_${job_name}.sbatch"
    slurm_output_log="${LOG_DIR_SMC}/${job_name}.%j.out"
    slurm_error_log="${LOG_DIR_SMC}/${job_name}.%j.err"

    # --- Generate the sbatch script ---
    cat << EOF > "$sbatch_script_path"
#!/bin/bash
# --- SLURM Directives (copied from main script config) ---
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --time=${SBATCH_TIME}
#SBATCH --nodes=${SBATCH_NODES}
#SBATCH --ntasks=${SBATCH_NTASKS}
#SBATCH --cpus-per-task=${SBATCH_CPUS_PER_TASK}
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
echo "Input VCF Directory for Pop: ${pop_input_vcf_dir}"
echo "Output SMC Directory for Pop: ${pop_smc_output_dir}"
echo "SMC++ Cutoff (-c): ${SMC_CUTOFF}"
echo "========================================================"

# --- Environment Setup & Tool Check ---
echo "Checking for required tools..."
command -v "${SMC_PLUS_PLUS_PATH}" >/dev/null 2>&1 || { echo "Error: smc++ command ('${SMC_PLUS_PLUS_PATH}') not found." >&2; exit 1; }
command -v "${BCFTOOLS_PATH}" >/dev/null 2>&1 || { echo "Error: bcftools command ('${BCFTOOLS_PATH}') not found." >&2; exit 1; }
echo "All tools found."

# --- Create Population-Specific Output Directory for SMC files ---
echo "Creating output directory for SMC files: ${pop_smc_output_dir}"
mkdir -p "${pop_smc_output_dir}"

# --- Loop through Chromosomes 1-22 ---
job_failed_chr_count=0
processed_chr_count=0
for chr_num in {1..22}; do
    echo -e "\\n--- Processing Chromosome \${chr_num} for Population ${pop} ---"

    # Define input VCF path for this pop & chr
    # Assumes naming convention like POP.chrNUM.vcf.gz (e.g., IND.chr1.vcf.gz)
    input_vcf_gz_path="${pop_input_vcf_dir}/${pop}.chr\${chr_num}.vcf.gz"

    # Define output SMC.gz path for this pop & chr
    output_smc_gz_path="${pop_smc_output_dir}/${pop}.chr\${chr_num}.smc.gz"

    echo "  Input VCF: \${input_vcf_gz_path}"
    echo "  Output SMC: \${output_smc_gz_path}"

    if [ ! -f "\${input_vcf_gz_path}" ]; then
        echo "  Warning: Input VCF not found: \${input_vcf_gz_path}. Skipping chromosome \${chr_num}." >&2
        ((job_failed_chr_count++))
        continue
    fi
    if [ ! -f "\${input_vcf_gz_path}.tbi" ]; then
        echo "  Warning: Index file (.tbi) for VCF not found: \${input_vcf_gz_path}.tbi. Skipping chromosome \${chr_num}." >&2
        echo "  Ensure VCFs are tabix indexed. This script does not re-index." >&2
        ((job_failed_chr_count++))
        continue
    fi

    # --- Extract individual IDs from the VCF header for the current chromosome ---
    echo "  Extracting individual list for ${pop}, Chr \${chr_num}..."
    IND_LIST=\$("${BCFTOOLS_PATH}" query -l "\${input_vcf_gz_path}" | paste -sd, -)
    if [ -z "\$IND_LIST" ]; then
        echo "  Error: Failed to extract individual list from VCF \${input_vcf_gz_path} or no samples found." >&2
        ((job_failed_chr_count++))
        continue
    fi
    # echo "    Individuals: \$IND_LIST" # Uncomment for debugging

    # --- Run smc++ vcf2smc ---
    echo "  Running smc++ vcf2smc for ${pop}, Chr \${chr_num}..."
    if ! "${SMC_PLUS_PLUS_PATH}" vcf2smc -c "${SMC_CUTOFF}" \\
        "\${input_vcf_gz_path}" \\
        "\${output_smc_gz_path}" \\
        "\${chr_num}" \\
        "${pop}:\${IND_LIST}"; then
        
        echo "  Error: smc++ vcf2smc failed for Chr \${chr_num}." >&2
        rm -f "\${output_smc_gz_path}" # Clean up potentially incomplete output
        ((job_failed_chr_count++))
        continue # Move to next chromosome
    fi
    echo "  smc++ vcf2smc successful for Chr \${chr_num}. Output: \$(basename \${output_smc_gz_path})"
    ((processed_chr_count++))

done # --- End of chromosome loop ---

# --- Final Job Status ---
echo "========================================================"
if [ \${job_failed_chr_count} -eq 0 ] && [ \${processed_chr_count} -gt 0 ]; then
    echo "SLURM Job for Population ${pop} (vcf2smc) finished successfully for \${processed_chr_count} chromosome(s)."
    exit 0
elif [ \${processed_chr_count} -eq 0 ]; then
    echo "SLURM Job for Population ${pop} (vcf2smc) processed 0 chromosomes. Check input VCFs and logs." >&2
    exit 1
else
    echo "SLURM Job for Population ${pop} (vcf2smc) finished with \${job_failed_chr_count} chromosome(s) failing out of $((processed_chr_count + job_failed_chr_count)) attempted." >&2
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
        echo "    Submitted SMC++ vcf2smc job for ${pop}. Job ID: ${job_id}"
        job_ids+=("${job_id}")
    else
        echo "    Error submitting SMC++ vcf2smc job for population ${pop}: ${sbatch_output}" >&2
    fi
done # --- End of population loop ---

echo "--------------------------------------------------"
if [ ${#job_ids[@]} -gt 0 ]; then
    echo "All ${#job_ids[@]} population SMC++ vcf2smc jobs submitted."
    echo "Monitor jobs using: squeue -u $USER"
    echo "Check SLURM output/error logs in: $LOG_DIR_SMC (e.g., POP_smcpp_vcf2smc.%j.out)"
    echo "Final SMC files will be in: ${SMC_OUTPUT_BASE_DIR}/{pop}/{pop}.chr{chr}.smc.gz"
else
    echo "No SMC++ vcf2smc jobs were submitted. Check for errors in population discovery or VCF directory availability."
fi
echo "--------------------------------------------------"

exit 0
