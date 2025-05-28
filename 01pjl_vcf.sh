#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status
set -o pipefail # Causes a pipeline to return the exit status of the last command that exited with a non-zero status

# --- Configuration ---
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# --- Target Population and 1000 Genomes Specifics ---
TARGET_POP="PJL"
PANEL_FILE="/home/tx56/1KGP/integrated_call_samples_v3.20130502.ALL.panel"
ONEKG_VCF_TEMPLATE="/home/tx56/1KGP/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
EXCLUDE_RSID_FILE="/home/tx56/1KGP/exclude_rsid.txt" 

# --- Output Configuration ---
OUTPUT_POP_DIR="${PWD}/${TARGET_POP}"

# --- SLURM Configuration ---
LOG_DIR_BASE="${PWD}/logs_extract_${TARGET_POP}"
SBATCH_SCRIPT_DIR="${LOG_DIR_BASE}/sbatch_scripts"
SBATCH_PARTITION="week"
SBATCH_TIME="2-00:00:00" 
SBATCH_NODES="1"
SBATCH_NTASKS="1"
SBATCH_CPUS="1" 
SBATCH_MEM_PER_CPU="80000" 

# --- Tools (Ensure they are in PATH within SLURM jobs) ---
VCFTOOLS_PATH="vcftools"
BCFTOOLS_PATH="bcftools" 
BGZIP_PATH="bgzip"
TABIX_PATH="tabix"
AWK_PATH="awk"
# =================================================

echo "--- VCF Extraction for Population ${TARGET_POP} from 1000 Genomes ---"
echo "Target Population: $TARGET_POP"
echo "Panel File: $PANEL_FILE"
echo "1000G VCF Template: $ONEKG_VCF_TEMPLATE"
echo "Exclude RSID File: $EXCLUDE_RSID_FILE"
echo "Output Directory: $OUTPUT_POP_DIR"
echo "Log Directory: $LOG_DIR_BASE"
echo "--------------------------------------------------"

if [ ! -f "$PANEL_FILE" ]; then
    echo "Error: Panel file not found: $PANEL_FILE" >&2
    exit 1
fi
if [ ! -f "$EXCLUDE_RSID_FILE" ]; then
    echo "Error: Exclude RSID file not found: $EXCLUDE_RSID_FILE" >&2
    exit 1
fi
echo "Panel file and Exclude RSID file found."
echo "Ensuring output and log directories exist..."
mkdir -p "$OUTPUT_POP_DIR"
mkdir -p "$LOG_DIR_BASE"
mkdir -p "$SBATCH_SCRIPT_DIR"
echo "Output and log directories checked/created."
pop_sample_list_path="${PWD}/${TARGET_POP}.samples.list"
echo "Generating sample list for ${TARGET_POP} from ${PANEL_FILE}..."
"${AWK_PATH}" -F'\t' -v pop="$TARGET_POP" 'NR > 1 && $2 == pop {print $1}' "$PANEL_FILE" > "$pop_sample_list_path"
if [ ! -s "$pop_sample_list_path" ]; then
    echo "Error: No samples found for population ${TARGET_POP} in ${PANEL_FILE}, or sample list is empty." >&2
    echo "Generated sample list path: ${pop_sample_list_path}"
    rm -f "$pop_sample_list_path"
    exit 1
fi
echo "Sample list created: $pop_sample_list_path (with $(wc -l < "$pop_sample_list_path") samples)"
echo "Generating sbatch script for population ${TARGET_POP}..."
job_name="${TARGET_POP}_extract_1kg_vcf_filtered" 
sbatch_script_path="${SBATCH_SCRIPT_DIR}/submit_${job_name}.sbatch"
slurm_output_log="${LOG_DIR_BASE}/${job_name}.%j.out"
slurm_error_log="${LOG_DIR_BASE}/${job_name}.%j.err"

cat << EOF > "$sbatch_script_path"
#!/bin/bash
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --time=${SBATCH_TIME} 
#SBATCH --nodes=${SBATCH_NODES}
#SBATCH --ntasks=${SBATCH_NTASKS}
#SBATCH --cpus-per-task=${SBATCH_CPUS}
#SBATCH --mem-per-cpu=${SBATCH_MEM_PER_CPU}
#SBATCH --job-name=${job_name} 
#SBATCH --output=${slurm_output_log} 
#SBATCH --error=${slurm_error_log} 

set -o pipefail 
set -x 

log_exit() {
  local exit_code=\$? 
  echo "--- SCRIPT EXITING (log_exit trap) ---" >&2 
  echo "Timestamp: \$(date)" >&2
  echo "Exit Code: \${exit_code}" >&2
  echo "Last command attempted (BASH_COMMAND): \${BASH_COMMAND}" >&2
  echo "Line number of BASH_COMMAND: \${BASH_LINENO}" >&2
  echo "Current function (FUNCNAME): \${FUNCNAME[*]}" >&2
  echo "Call stack (BASH_SOURCE BASH_LINENO FUNCNAME):" >&2
  for i in "\${!BASH_SOURCE[@]}"; do
    echo "  \${BASH_SOURCE[i]}:\${BASH_LINENO[i]} \${FUNCNAME[i]}" >&2
  done
  echo "Current working directory: \$(pwd)" >&2
  echo "--- END SCRIPT EXITING INFO ---" >&2
}
trap log_exit EXIT
set -e

echo "========================================================"
echo "SLURM Job ID: \$SLURM_JOB_ID; Node: \$(hostname); Shell: \$SHELL; Bash Version: \$BASH_VERSION"
echo "Current Date/Time: \$(date)"
echo "Processing Population ID: ${TARGET_POP}"
echo "Using Sample List: ${pop_sample_list_path}" 
echo "Using Exclude RSID list: ${EXCLUDE_RSID_FILE}"
echo "1000G VCF Template: ${ONEKG_VCF_TEMPLATE}"
echo "Output Directory for this Pop: ${OUTPUT_POP_DIR}"
echo "========================================================"

echo "Checking for required tools..."
command -v "${VCFTOOLS_PATH}" >/dev/null 2>&1 || { echo "Error: VCFtools command ('${VCFTOOLS_PATH}') not found." >&2; exit 1; }
command -v "${BCFTOOLS_PATH}" >/dev/null 2>&1 || { echo "Error: bcftools command ('${BCFTOOLS_PATH}') not found." >&2; exit 1; }
command -v "${BGZIP_PATH}" >/dev/null 2>&1 || { echo "Error: bgzip command ('${BGZIP_PATH}') not found." >&2; exit 1; }
command -v "${TABIX_PATH}" >/dev/null 2>&1 || { echo "Error: tabix command ('${TABIX_PATH}') not found." >&2; exit 1; }
echo "All tools found."

job_failed_chr_count=0
processed_chr_count=0
echo "DEBUG: Initializing job_failed_chr_count=\${job_failed_chr_count}, processed_chr_count=\${processed_chr_count}"
current_script_stage="Loop Initialization" 

for chr_num in \$(seq 1 22); do 
    current_script_stage="Top of loop for chr\${chr_num}"
    echo "DEBUG: Top of loop for chr_num=\${chr_num}"
    echo -e "\\n--- Processing Chromosome \${chr_num} for Population ${TARGET_POP} ---"
    onekg_vcf_gz_path=\$(echo "${ONEKG_VCF_TEMPLATE}" | sed "s/{chr}/\${chr_num}/g")
    current_script_stage="Constructed VCF path for chr\${chr_num}"
    output_vcf_gz_path="${OUTPUT_POP_DIR}/${TARGET_POP}.chr\${chr_num}.vcf.gz"
    pipeline_log_prefix="${OUTPUT_POP_DIR}/${TARGET_POP}.chr\${chr_num}.pipeline"
    echo "  Source 1000G VCF: \${onekg_vcf_gz_path}"
    echo "  Output VCF: \${output_vcf_gz_path}"
    echo "  Sample List: ${pop_sample_list_path}" 
    echo "  Exclude RSID List: ${EXCLUDE_RSID_FILE}" 
    current_script_stage="File checks for chr\${chr_num}"
    if [ ! -f "\${onekg_vcf_gz_path}" ]; then
        echo "  CRITICAL WARNING (stderr): Source VCF file NOT FOUND: \${onekg_vcf_gz_path}. Skipping chromosome \${chr_num}." >&2
        job_failed_chr_count=\$((job_failed_chr_count + 1))
        echo "DEBUG: job_failed_chr_count is now \${job_failed_chr_count} (VCF missing for chr\${chr_num})"
        continue
    fi
    if [ ! -f "\${onekg_vcf_gz_path}.tbi" ]; then
        echo "  CRITICAL WARNING (stderr): Index for source VCF file NOT FOUND: \${onekg_vcf_gz_path}.tbi. Skipping chromosome \${chr_num}." >&2
        job_failed_chr_count=\$((job_failed_chr_count + 1))
        echo "DEBUG: job_failed_chr_count is now \${job_failed_chr_count} (TBI missing for chr\${chr_num})"
        continue
    fi
    echo "DEBUG: Source VCF and index appear to exist for chr_num=\${chr_num}."
    echo "  Extracting samples (vcftools), filtering (bcftools view), sorting (bcftools sort), and compressing (bgzip)..."
    current_script_stage="VCFtools | BCFtools view | BCFtools sort | BGZIP pipeline for chr\${chr_num}"
    
    ( "${VCFTOOLS_PATH}" --gzvcf "\${onekg_vcf_gz_path}" \
        --keep "${pop_sample_list_path}" \
        --recode --recode-INFO-all \
        --stdout | \
      "${BCFTOOLS_PATH}" view -m2 -M2 --exclude "ID=@${EXCLUDE_RSID_FILE}" -O u - | \
      "${BCFTOOLS_PATH}" sort -O v -T "\${SLURM_TMPDIR:-\${PWD}}/bcftools_sort_tmp_\${SLURM_JOB_ID}_\${chr_num}" | \
      "${BGZIP_PATH}" -c > "\${output_vcf_gz_path}" \
    ) 2> "\${pipeline_log_prefix}.log"
    pipeline_exit_code=\$? 

    if [ \${pipeline_exit_code} -ne 0 ]; then
        echo "  ERROR (stderr): VCF processing pipeline FAILED for Chr \${chr_num} with exit code \${pipeline_exit_code}." >&2
        echo "  Check pipeline log: \${pipeline_log_prefix}.log"
        if [ -f "\${pipeline_log_prefix}.log" ]; then
             echo "--- Content of \${pipeline_log_prefix}.log ---" >&2
             cat "\${pipeline_log_prefix}.log" >&2
             echo "--- End of \${pipeline_log_prefix}.log ---" >&2
        fi
        rm -f "\${output_vcf_gz_path}" 
        job_failed_chr_count=\$((job_failed_chr_count + 1))
        echo "DEBUG: job_failed_chr_count is now \${job_failed_chr_count} after pipeline failure for chr\${chr_num}"
        exit 1 
    else
        echo "  VCF processing pipeline successful for Chr \${chr_num}."
        if [ -s "\${pipeline_log_prefix}.log" ]; then
            if ! grep -q -E -i "Error|Fail|Unable|Cannot|Warning" "\${pipeline_log_prefix}.log"; then
                 rm -f "\${pipeline_log_prefix}.log"
            else
                echo "  Pipeline log \${pipeline_log_prefix}.log contains potential issues, keeping it."
            fi
        else 
            rm -f "\${pipeline_log_prefix}.log"
        fi
    fi
    current_script_stage="Tabix for chr\${chr_num}"
    echo "  Checking if output VCF has data before indexing..."
    data_line_count=0 
    if [ -f "\${output_vcf_gz_path}" ] && [ -s "\${output_vcf_gz_path}" ]; then
        set +e
        data_line_count=\$(zcat "\${output_vcf_gz_path}" 2>/dev/null | grep -cv "^#")
        zcat_exit_status=\$?
        set -e
        if [ \${zcat_exit_status} -ne 0 ] && [ \${zcat_exit_status} -ne 141 ]; then 
            echo "  WARNING: zcat on \${output_vcf_gz_path} encountered an issue (exit code \${zcat_exit_status}), assuming 0 data lines." >&2
            data_line_count=0
        elif [ -z "\${data_line_count}" ]; then 
             data_line_count=0
        fi
    else
        echo "  Output VCF \${output_vcf_gz_path} does not exist or is empty (0 bytes)."
        data_line_count=0
    fi
    echo "  Number of data lines found: \${data_line_count}"
    if [ \${data_line_count} -gt 0 ]; then
        echo "  Output VCF \${output_vcf_gz_path} has \${data_line_count} data lines. Indexing with tabix..."
        "${TABIX_PATH}" -p vcf "\${output_vcf_gz_path}"
        tabix_exit_code=\$?
        if [ \${tabix_exit_code} -ne 0 ] || [ ! -f "\${output_vcf_gz_path}.tbi" ]; then
            echo "  ERROR (stderr): tabix FAILED for Chr \${chr_num} on \${output_vcf_gz_path}." >&2
            echo "  Tabix exit code: \${tabix_exit_code}. Index file exists: $( [ -f "\${output_vcf_gz_path}.tbi" ] && echo "Yes" || echo "No" )" >&2
            rm -f "\${output_vcf_gz_path}.tbi" 
            job_failed_chr_count=\$((job_failed_chr_count + 1))
            echo "DEBUG: job_failed_chr_count is now \${job_failed_chr_count} after tabix failure for chr\${chr_num}"
            if [ \${tabix_exit_code} -ne 0 ]; then 
                exit 1 
            fi
        else
            echo "  Tabix indexing successful for Chr \${chr_num}."
        fi
    else
        echo "  WARNING: Output VCF \${output_vcf_gz_path} is empty or header-only (0 data lines). Skipping tabix indexing for Chr \${chr_num}."
    fi
    current_script_stage="Incrementing processed_chr_count for chr\${chr_num}"
    processed_chr_count=\$((processed_chr_count + 1))
    job_failed_chr_count=\$((job_failed_chr_count + 0)) 
    echo "DEBUG: processed_chr_count is now \${processed_chr_count}"
    echo "DEBUG: End of loop body for chr_num=\${chr_num}."
done
current_script_stage="After chromosome loop, preparing final status"
echo "DEBUG: Chromosome loop finished. Processed: \${processed_chr_count}, Failed: \${job_failed_chr_count}"
echo "========================================================"
if [ \${job_failed_chr_count} -eq 0 ] && [ \${processed_chr_count} -eq 22 ]; then 
    echo "SLURM Job for Population ${TARGET_POP} (1000G extraction & filtering) finished successfully for all \${processed_chr_count} chromosome(s)."
    current_script_stage="Exiting successfully"
    exit 0
elif [ \${processed_chr_count} -eq 0 ] && [ \${job_failed_chr_count} -gt 0 ]; then
    echo "SLURM Job for Population ${TARGET_POP} (1000G extraction & filtering) processed 0 chromosomes and \${job_failed_chr_count} chromosomes failed to start (e.g. missing input)." >&2
    current_script_stage="Exiting - no chromosomes processed, initial failures"
    exit 1
elif [ \${job_failed_chr_count} -gt 0 ]; then
    echo "SLURM Job for Population ${TARGET_POP} (1000G extraction & filtering) finished with issues. Processed: \${processed_chr_count}, Failed during processing: \${job_failed_chr_count}." >&2
    current_script_stage="Exiting - some chromosomes failed during processing"
    exit 1
elif [ \${processed_chr_count} -lt 22 ] && [ \${job_failed_chr_count} -eq 0 ]; then 
    echo "SLURM Job for Population ${TARGET_POP} (1000G extraction & filtering) processed only \${processed_chr_count} out of 22 expected chromosomes, though no explicit processing errors were reported for these. This might be due to skipping chromosomes with missing inputs." >&2
    current_script_stage="Exiting - incomplete processing without explicit processing errors"
    exit 1   
else 
    echo "SLURM Job for Population ${TARGET_POP} (1000G extraction & filtering) finished with an undetermined status. Processed: \${processed_chr_count}, Failed: \${job_failed_chr_count}." >&2
    current_script_stage="Exiting - undetermined status"
    exit 1
fi
echo "========================================================"
current_script_stage="Very end of script - should not see this"
EOF
# --- End of Here Document ---

chmod +x "$sbatch_script_path"
echo "  Generated sbatch script: $sbatch_script_path"
sbatch_output=$(sbatch "$sbatch_script_path")
job_id=$(echo "$sbatch_output" | awk '{print $NF}')
if [[ "$job_id" =~ ^[0-9]+$ ]]; then
    echo "  Submitted job ID for ${TARGET_POP} extraction: ${job_id}"
    echo "  Temporary sample list ${pop_sample_list_path} will be used by the job."
    echo "  Consider removing it manually after the job completes if not needed: rm -f \"${pop_sample_list_path}\""
else
    echo "  Error submitting job for population ${TARGET_POP} extraction: ${sbatch_output}" >&2
    rm -f "$pop_sample_list_path"
    echo "  Cleaned up temporary sample list: $pop_sample_list_path"
fi
echo "--------------------------------------------------"
echo "Job submission process for ${TARGET_POP} complete."
echo "Monitor job using: squeue -u $USER -n ${job_name}"
echo "Check SLURM output/error logs in: $LOG_DIR_BASE"
echo "Final VCFs will be in: ${OUTPUT_POP_DIR}/${TARGET_POP}.chr{chr}.vcf.gz"
echo "--------------------------------------------------"

exit 0
