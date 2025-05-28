#!/bin/bash

# --- SLURM Configuration ---
#SBATCH --partition=week
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=smc_process
#SBATCH --output=slurm_logs_placeholder/%x_%j.out # Placeholder - see sbatch command
#SBATCH --error=slurm_logs_placeholder/%x_%j.err  # Placeholder - see sbatch command

set -e
set -o pipefail

# --- Argument Parsing ---
if [ -z "$1" ]; then
    echo "ERROR (${BASH_SOURCE[0]}): Population label (arg 1) is missing." >&2
    exit 1
fi
POP_LABEL_TO_USE="$1"
SPECIFIC_CHR_NUM_ARG="${2:-}"

# --- Base Directory for All Operations: SLURM_SUBMIT_DIR ---
# If running under Slurm, use SLURM_SUBMIT_DIR. Otherwise, fall back to PWD for local testing.
if [ -n "$SLURM_JOB_ID" ]; then
    BASE_OPERATING_DIR="$SLURM_SUBMIT_DIR"
    echo "INFO (${BASH_SOURCE[0]}): Running under Slurm. BASE_OPERATING_DIR set to SLURM_SUBMIT_DIR: ${BASE_OPERATING_DIR}"
    # Set SMC_CORES from Slurm allocation
    export SMC_CORES=${SLURM_CPUS_PER_TASK}
    echo "INFO (${BASH_SOURCE[0]}): SMC_CORES set to ${SMC_CORES} from SLURM_CPUS_PER_TASK"
    # Create Slurm log directory based on BASE_OPERATING_DIR (SLURM_SUBMIT_DIR)
    # This is where the #SBATCH --output/error directives should point.
    # The actual log file names are controlled by #SBATCH or sbatch command line.
    # This mkdir is just to ensure the directory exists if using relative paths in #SBATCH.
    mkdir -p "${BASE_OPERATING_DIR}/slurm_logs" # For example, if #SBATCH --output=slurm_logs/...
else
    # Not under Slurm (e.g., local bash execution for testing)
    BASE_OPERATING_DIR=$(pwd)
    echo "INFO (${BASH_SOURCE[0]}): Not running under Slurm. BASE_OPERATING_DIR set to PWD: ${BASE_OPERATING_DIR}"
    SMC_CORES=${SMC_CORES:-4} # Default for local run
    echo "INFO (${BASH_SOURCE[0]}): SMC_CORES set to ${SMC_CORES} (default or pre-existing env var)"
fi

# --- Path and File Name Configuration (Inputs) ---
# All relative paths for inputs/outputs are now relative to BASE_OPERATING_DIR
INPUT_VCF_DIR_RELPATH="${POP_LABEL_TO_USE}" # e.g. PJL/ (relative to BASE_OPERATING_DIR)
INPUT_VCF_PREFIX="${POP_LABEL_TO_USE}.chr"
INPUT_VCF_SUFFIX=".vcf.gz"
AA_TSV_DIR_ABS_PATH="/vast/palmer/pi/reilly/jfa38/datasets_for_annotation/ErinG_ARG_hg19_AncestralAlleles/tidied_AA_tsvs" # Keep absolute if central
AA_TSV_BASENAME_PREFIX="chr"
AA_TSV_BASENAME_SUFFIX="_ErinG_AA_hg19.tsv.gz"
FULL_MASK_FILE_RELPATH="20141020.strict_mask.whole_genome.bed" # Relative to BASE_OPERATING_DIR
SAMPLES_FILE_RELPATH="${POP_LABEL_TO_USE}.samples.list"       # Relative to BASE_OPERATING_DIR
REFERENCE_FASTA=""
UNIFORM_RECOMBINATION_RATE="1.25e-8"
PYTHON_ANNOTATE_SCRIPT_NAME="annotate_aa.py" # Assumed to be in BASE_OPERATING_DIR
PYTHON_ANNOTATE_SCRIPT_PATH="${BASE_OPERATING_DIR}/${PYTHON_ANNOTATE_SCRIPT_NAME}"

# --- Output Directory Structure Definition (relative to BASE_OPERATING_DIR) ---
POP_OUTPUT_DIR_ROOT="${BASE_OPERATING_DIR}/${POP_LABEL_TO_USE}_smc_processing_results"
VCF_SMC_INPUT_DIR="${POP_OUTPUT_DIR_ROOT}/vcf_to_smc_input"
SMC_ESTIMATE_INPUT_DIR="${POP_OUTPUT_DIR_ROOT}/smc_estimate_input"
SMC_ESTIMATED_MODEL_PER_CHR_BASE_DIR="${POP_OUTPUT_DIR_ROOT}/smc_estimated_initial_model_per_chr"
INTERMEDIATES_BASE_DIR="${POP_OUTPUT_DIR_ROOT}/all_intermediate_steps"

mkdir -p "${VCF_SMC_INPUT_DIR}"
mkdir -p "${SMC_ESTIMATE_INPUT_DIR}"
mkdir -p "${SMC_ESTIMATED_MODEL_PER_CHR_BASE_DIR}"
mkdir -p "${INTERMEDIATES_BASE_DIR}"

# --- Global definitions ---
SAMPLES_FILE_ABS_PATH="${BASE_OPERATING_DIR}/${SAMPLES_FILE_RELPATH}"
FULL_MASK_FILE_ABS_PATH="${BASE_OPERATING_DIR}/${FULL_MASK_FILE_RELPATH}"
if [ ! -f "${SAMPLES_FILE_ABS_PATH}" ]; then echo "ERROR: Samples list ${SAMPLES_FILE_ABS_PATH} not found." >&2; exit 1; fi
mapfile -t SAMPLES_POP < <(cat "${SAMPLES_FILE_ABS_PATH}")
if [ ${#SAMPLES_POP[@]} -lt 2 ]; then echo "ERROR: Not enough samples in ${SAMPLES_FILE_ABS_PATH}." >&2; exit 1; fi
DISTINGUISHED_INDIVIDUALS=("${SAMPLES_POP[0]}" "${SAMPLES_POP[1]}")
if [ ! -f "${FULL_MASK_FILE_ABS_PATH}" ]; then echo "ERROR: Full genome mask not found: ${FULL_MASK_FILE_ABS_PATH}" >&2; exit 1; fi
if [ ! -f "${PYTHON_ANNOTATE_SCRIPT_PATH}" ]; then echo "ERROR: Python script ${PYTHON_ANNOTATE_SCRIPT_PATH} not found." >&2; exit 1; fi


# Function to perform all processing for a single chromosome
process_single_chromosome() {
    local current_chr_num="$1"
    local current_target_chr_id="chr${current_chr_num}"
    local current_log_prefix_msg="POP: ${POP_LABEL_TO_USE}, ${current_target_chr_id}"

    log_chr() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] (${current_log_prefix_msg}) $1"; }
    log_chr "--- Starting processing for chromosome ${current_chr_num} ---"

    local local_smc_estimated_model_this_chr_dir="${SMC_ESTIMATED_MODEL_PER_CHR_BASE_DIR}/chr${current_chr_num}"
    local local_intermediates_this_chr_dir="${INTERMEDIATES_BASE_DIR}/chr${current_chr_num}"
    mkdir -p "${local_smc_estimated_model_this_chr_dir}" "${local_intermediates_this_chr_dir}"

    local local_original_input_vcf_abs_path="${BASE_OPERATING_DIR}/${INPUT_VCF_DIR_RELPATH}/${INPUT_VCF_PREFIX}${current_chr_num}${INPUT_VCF_SUFFIX}"
    local local_aa_tsv_original_abs_path="${AA_TSV_DIR_ABS_PATH}/${AA_TSV_BASENAME_PREFIX}${current_chr_num}${AA_TSV_BASENAME_SUFFIX}" # AA TSV path can remain absolute

    local local_vcf_body_renamed="${local_intermediates_this_chr_dir}/${POP_LABEL_TO_USE}.chr${current_chr_num}.phased.bodyRenamed.vcf.gz"
    local local_temp_renamed_header_file="${local_intermediates_this_chr_dir}/temp_renamed_header_chr${current_chr_num}.txt"
    local local_temp_chrom_map_file="${local_intermediates_this_chr_dir}/temp_chromosome_map_chr${current_chr_num}.txt"
    local local_vcf_chr_normalized_raw="${local_intermediates_this_chr_dir}/${POP_LABEL_TO_USE}.chr${current_chr_num}.phased.chrNormalized.raw.vcf.gz"
    local local_vcf_chr_normalized_cleaned="${local_intermediates_this_chr_dir}/${POP_LABEL_TO_USE}.chr${current_chr_num}.phased.chrNormalized.cleaned.vcf.gz"
    local local_aa_tsv_prepared="${local_intermediates_this_chr_dir}/chr${current_chr_num}_AA_prepared.tsv.gz"
    local local_vcf_annotated_py="${local_intermediates_this_chr_dir}/${POP_LABEL_TO_USE}.chr${current_chr_num}.phased.annotated.py.vcf.gz"
    local local_mask_file_final="${local_intermediates_this_chr_dir}/mask.${current_target_chr_id}.bed.gz"
    local local_python_stdout_log="${local_intermediates_this_chr_dir}/python_annotate_stdout_chr${current_chr_num}.log"
    local local_python_stderr_log="${local_intermediates_this_chr_dir}/python_annotate_stderr_chr${current_chr_num}.log"
    local local_final_vcf_for_smc_tool_input="${VCF_SMC_INPUT_DIR}/${POP_LABEL_TO_USE}.chr${current_chr_num}.phased.annotated.py.AAfiltered.vcf.gz"
    local local_smc_format_file_for_estimate="${SMC_ESTIMATE_INPUT_DIR}/${POP_LABEL_TO_USE}.chr${current_chr_num}.smc.gz"

    # Step 1: Normalize VCF
    log_chr "--- Starting Step 1: Normalize Input VCF ---"
    if [ ! -f "${local_original_input_vcf_abs_path}" ]; then log_chr "ERROR: VCF not found: ${local_original_input_vcf_abs_path}" >&2; return 1; fi
    if [ ! -f "${local_original_input_vcf_abs_path}.tbi" ] && [ ! -f "${local_original_input_vcf_abs_path}.csi" ]; then log_chr "ERROR: Index for VCF ${local_original_input_vcf_abs_path} not found." >&2; return 1; fi
    local original_vcf_chrom_name_sample=""
    set +e # Temporarily disable exit on error
    _command_output=$(zcat "${local_original_input_vcf_abs_path}" 2>/dev/null | grep -v "^#" | head -n 1 | awk '{print $1}')
    _pipe_statuses=("${PIPESTATUS[@]}")
    set -e # Re-enable exit on error
    if [ -n "$_command_output" ]; then original_vcf_chrom_name_sample="$_command_output"; else log_chr "ERROR: Could not get chr name from ${local_original_input_vcf_abs_path}. PIPESTATUS: ${_pipe_statuses[*]}" >&2; zcat "${local_original_input_vcf_abs_path}" | head -n 20 >&2; return 1; fi
    log_chr "Original VCF chr name sample: '${original_vcf_chrom_name_sample}'"
    if [ "${original_vcf_chrom_name_sample}" != "${current_target_chr_id}" ]; then
        log_chr "Normalizing chr name from '${original_vcf_chrom_name_sample}' to '${current_target_chr_id}'"
        echo -e "${original_vcf_chrom_name_sample}\t${current_target_chr_id}" > "${local_temp_chrom_map_file}"
        bcftools annotate --rename-chrs "${local_temp_chrom_map_file}" "${local_original_input_vcf_abs_path}" -O z -o "${local_vcf_body_renamed}" --threads "${SMC_CORES}"
        bcftools view -h "${local_vcf_body_renamed}" > "${local_temp_renamed_header_file}"
        sed -E -e "s/##contig=<ID=${original_vcf_chrom_name_sample}([,>])/##contig=<ID=${current_target_chr_id}\\1/g" "${local_temp_renamed_header_file}" > "${local_temp_renamed_header_file}.final"
        mv "${local_temp_renamed_header_file}.final" "${local_temp_renamed_header_file}"
        if ! grep -q "##contig=<ID=${current_target_chr_id}" "${local_temp_renamed_header_file}"; then echo "##contig=<ID=${current_target_chr_id}>" >> "${local_temp_renamed_header_file}"; fi
        bcftools reheader -h "${local_temp_renamed_header_file}" "${local_vcf_body_renamed}" -o "${local_vcf_chr_normalized_raw}"
    else
        cp "${local_original_input_vcf_abs_path}" "${local_vcf_chr_normalized_raw}"
    fi
    bcftools norm -m -any "${local_vcf_chr_normalized_raw}" -O z -o "${local_vcf_chr_normalized_cleaned}" --threads "${SMC_CORES}"
    tabix -f -p vcf "${local_vcf_chr_normalized_cleaned}"
    local input_vcf_current_step="${local_vcf_chr_normalized_cleaned}"
    log_chr "--- Step 1 finished. ---"

    # Step 2: Prepare AA TSV
    log_chr "--- Starting Step 2: Prepare AA TSV ---"
    if [ ! -f "${local_aa_tsv_original_abs_path}" ]; then log_chr "ERROR: AA TSV not found: ${local_aa_tsv_original_abs_path}" >&2; return 1; fi
    zcat "${local_aa_tsv_original_abs_path}" | awk -v tc="${current_target_chr_id}" 'BEGIN{OFS="\t"} (NR==1 && $0 !~ /^#/){print "#CHROM_AA","POS_AA","AA_VALUE"} ($0 !~ /^#/ && NF>=3){print tc,$2,$3} ($0 ~ /^#/){print $0}' | sort -k1,1V -k2,2n | bgzip -c > "${local_aa_tsv_prepared}"
    ps2=("${PIPESTATUS[@]}"); if [ ${ps2[1]} -ne 0 ] || [ ${ps2[2]} -ne 0 ]; then log_chr "ERROR: AA TSV prep (awk/sort) failed. PIPESTATUS: ${ps2[*]}" >&2; return 1; fi
    tabix -f -S 1 -s 1 -b 2 -e 2 "${local_aa_tsv_prepared}"
    log_chr "--- Step 2 finished. ---"

    # Step 3: Python Annotate
    log_chr "--- Starting Step 3: Python Annotate ---"
    python3 "${PYTHON_ANNOTATE_SCRIPT_PATH}" "${input_vcf_current_step}" "${local_aa_tsv_prepared}" "${local_vcf_annotated_py}" 1> "${local_python_stdout_log}" 2> "${local_python_stderr_log}"
    py_exit=$?; if [ $py_exit -ne 0 ]; then log_chr "ERROR: Python script failed (status $py_exit)." >&2; cat "${local_python_stderr_log}" >&2 ; return 1; fi
    tabix -f -p vcf "${local_vcf_annotated_py}"
    input_vcf_current_step="${local_vcf_annotated_py}"
    log_chr "--- Step 3 finished. ---"

    # Step 4: Filter VCF by AA
    log_chr "--- Starting Step 4: Filter by AA ---"
    bcftools view -i 'INFO/AA != "." && INFO/AA != "N"' -O z -o "${local_final_vcf_for_smc_tool_input}" "${input_vcf_current_step}" --threads "${SMC_CORES}"
    tabix -f -p vcf "${local_final_vcf_for_smc_tool_input}"
    input_vcf_current_step="${local_final_vcf_for_smc_tool_input}"
    log_chr "--- Step 4 finished. ---"

    # Step 5: Prepare Mask
    log_chr "--- Starting Step 5: Prepare Mask ---"
    local local_mask_file_plain_temp="${local_intermediates_this_chr_dir}/mask.${current_target_chr_id}.temp.plain.bed"
    if [ ! -f "${FULL_MASK_FILE_ABS_PATH}" ]; then log_chr "ERROR: Mask file not found: ${FULL_MASK_FILE_ABS_PATH}" >&2; return 1; fi # FULL_MASK_FILE_ABS_PATH is global
    grep "^${current_target_chr_id}[[:space:]]" "${FULL_MASK_FILE_ABS_PATH}" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' | sort -k1,1V -k2,2n > "${local_mask_file_plain_temp}"
    ps5=("${PIPESTATUS[@]}"); if [ ${ps5[0]} -gt 1 ]; then log_chr "ERROR: grep for mask failed (status ${ps5[0]})." >&2; return 1; fi
    if [ ! -s "${local_mask_file_plain_temp}" ] && [ ${ps5[0]} -eq 1 ]; then log_chr "WARNING: No mask regions for ${current_target_chr_id}."; fi
    bgzip -f -c "${local_mask_file_plain_temp}" > "${local_mask_file_final}"; rm -f "${local_mask_file_plain_temp}"
    tabix -f -p bed "${local_mask_file_final}"
    log_chr "--- Step 5 finished. ---"

    # Step 6: smc++ vcf2smc
    log_chr "--- Starting Step 6: smc++ vcf2smc ---"
    local local_smc_sample_list_str="${POP_LABEL_TO_USE}:$(IFS=,; echo "${SAMPLES_POP[*]}")"
    smc++ vcf2smc --cores "${SMC_CORES}" -m "${local_mask_file_final}" -d "${DISTINGUISHED_INDIVIDUALS[0]}" "${DISTINGUISHED_INDIVIDUALS[1]}" "${input_vcf_current_step}" "${local_smc_format_file_for_estimate}" "${current_target_chr_id}" "${local_smc_sample_list_str}"
    log_chr "--- Step 6 finished. ---"

    # Step 7: smc++ estimate
    log_chr "--- Starting Step 7: smc++ estimate ---"
    if [ ! -f "${local_smc_format_file_for_estimate}" ] || [ ! -s "${local_smc_format_file_for_estimate}" ]; then log_chr "ERROR: SMC input not found: ${local_smc_format_file_for_estimate}" >&2; return 1; fi
    smc++ estimate --cores "${SMC_CORES}" -o "${local_smc_estimated_model_this_chr_dir}" "${UNIFORM_RECOMBINATION_RATE}" "${local_smc_format_file_for_estimate}"
    local local_smc_model_json_this_chr="${local_smc_estimated_model_this_chr_dir}/model.final.json"
    if [ ! -f "${local_smc_model_json_this_chr}" ]; then log_chr "ERROR: Model JSON not created: ${local_smc_model_json_this_chr}" >&2; return 1; fi
    log_chr "--- Step 7 finished. ---"

    # Cleanup
    log_chr "--- Cleaning up temporary files ---"
    rm -f "${local_temp_renamed_header_file}" "${local_temp_renamed_header_file}.final" \
          "${local_temp_chrom_map_file}" \
          "${local_vcf_body_renamed}" "${local_vcf_body_renamed}.tbi" "${local_vcf_body_renamed}.csi" \
          "${local_vcf_chr_normalized_raw}" "${local_vcf_chr_normalized_raw}.tbi" "${local_vcf_chr_normalized_raw}.csi" \
          2>/dev/null || log_chr "Note: Some temp files already cleaned or not created."
    log_chr "--- SUCCESS for chromosome ${current_chr_num}. ---"
    return 0
}

# --- Main Script Logic ---
log_main() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] (POP: ${POP_LABEL_TO_USE}, MainScript@${HOSTNAME}) $1"; }

log_main "--- Initializing script for POP: ${POP_LABEL_TO_USE} ---"
log_main "BASE_OPERATING_DIR is ${BASE_OPERATING_DIR}"
log_main "Python script: ${PYTHON_ANNOTATE_SCRIPT_PATH}"
log_main "SMC_CORES: ${SMC_CORES}"
log_main "Tool checks passed."
log_main "Total samples: ${#SAMPLES_POP[@]}, DIs: ${DISTINGUISHED_INDIVIDUALS[*]}"

if [ -n "${SPECIFIC_CHR_NUM_ARG}" ]; then
    log_main "Processing single specified chromosome: ${SPECIFIC_CHR_NUM_ARG}"
    process_single_chromosome "${SPECIFIC_CHR_NUM_ARG}"
else
    log_main "Processing all autosomes (1-22) sequentially for population ${POP_LABEL_TO_USE}."
    for chr_num_loop in $(seq 1 22); do
        echo # Newline for readability between chromosome logs
        process_single_chromosome "${chr_num_loop}"
        # If process_single_chromosome fails, 'set -e' will cause script to exit.
    done
fi

log_main "--- Script ${BASH_SOURCE[0]} finished successfully for POP: ${POP_LABEL_TO_USE}. ---"
exit 0
