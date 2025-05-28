#!/bin/bash

# --- SLURM Configuration ---
#SBATCH --partition=week
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=4   
#SBATCH --mem=32G
#SBATCH --job-name=smc_pwd_FIXED_AA
#SBATCH --output=slurm_logs/%x_%j.out 
#SBATCH --error=slurm_logs/%x_%j.err  

set -e
set -o pipefail

# --- Argument Parsing and Mode Detection ---
BASE_OPERATING_DIR=$(pwd) 
SCRIPT_ACTUAL_LOCATION="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)" 

POP_LABEL_TO_USE=""
SPECIFIC_CHR_NUM_ARG_FROM_CMD="" 
RUN_MODE=""

if [ -n "$SLURM_JOB_ID" ]; then
    RUN_MODE="Slurm"
    if [ -z "$1" ]; then echo "ERROR (Slurm): Population label (arg 1) is missing." >&2; exit 1; fi
    POP_LABEL_TO_USE="$1"
    SPECIFIC_CHR_NUM_ARG_FROM_CMD="${2:-}" 
    SMC_CORES=${SLURM_CPUS_PER_TASK:-4}
    echo "INFO (${BASH_SOURCE[0]}): Running under Slurm. POP: ${POP_LABEL_TO_USE}. Optional CHR: ${SPECIFIC_CHR_NUM_ARG_FROM_CMD}."
else
    RUN_MODE="Direct"
    if [ -z "$1" ]; then echo "ERROR (Direct): Population label (arg 1) is missing." >&2; exit 1; fi
    POP_LABEL_TO_USE="$1"
    SPECIFIC_CHR_NUM_ARG_FROM_CMD="${2:-}"
    SMC_CORES=${SMC_CORES:-4} 
    echo "INFO (${BASH_SOURCE[0]}): Direct Bash execution. POP: ${POP_LABEL_TO_USE}. Optional CHR: ${SPECIFIC_CHR_NUM_ARG_FROM_CMD}."
fi
export SMC_CORES=$(( SMC_CORES > 0 ? SMC_CORES : 1 ))
echo "INFO (${BASH_SOURCE[0]}): BASE_OPERATING_DIR (PWD of launch): ${BASE_OPERATING_DIR}"
echo "INFO (${BASH_SOURCE[0]}): Script location: ${SCRIPT_ACTUAL_LOCATION}"
echo "INFO (${BASH_SOURCE[0]}): SMC_CORES set to ${SMC_CORES}."

# --- Path and File Name Configuration (Inputs) ---
PYTHON_ANNOTATE_SCRIPT_NAME="annotate_aa.py"
PYTHON_FILTER_SCRIPT_NAME="filter_vcf_for_smc.py" 
PYTHON_ANNOTATE_SCRIPT_PATH="${BASE_OPERATING_DIR}/${PYTHON_ANNOTATE_SCRIPT_NAME}"
PYTHON_FILTER_SCRIPT_PATH="${BASE_OPERATING_DIR}/${PYTHON_FILTER_SCRIPT_NAME}"

INPUT_VCF_DIR_RELPATH="${POP_LABEL_TO_USE}" 
INPUT_VCF_PREFIX="${POP_LABEL_TO_USE}.chr"
INPUT_VCF_SUFFIX=".vcf.gz"
AA_TSV_DIR_ABS_PATH="/vast/palmer/pi/reilly/jfa38/datasets_for_annotation/ErinG_ARG_hg19_AncestralAlleles/tidied_AA_tsvs" 
AA_TSV_BASENAME_PREFIX="chr"
AA_TSV_BASENAME_SUFFIX="_ErinG_AA_hg19.tsv.gz"
FULL_MASK_FILE_RELPATH="20141020.strict_mask.whole_genome.bed" 
REFERENCE_FASTA="" 
UNIFORM_RECOMBINATION_RATE="1.25e-8"
EXCLUDE_RSID_FILE_ABS_PATH="/home/tx56/1KGP/exclude_rsid.txt" 

# --- Output Directory Structure Definition (All relative to BASE_OPERATING_DIR) ---
POP_OUTPUT_DIR_ROOT="${BASE_OPERATING_DIR}/${POP_LABEL_TO_USE}_smc_processing_results_pwd_fix_aa" 
VCF_FINAL_FOR_SMC_DIR="${POP_OUTPUT_DIR_ROOT}/vcf_final_for_smc" 
SMC_FORMAT_FILES_DIR="${POP_OUTPUT_DIR_ROOT}/smc_format_files" 
SMC_MODEL_PER_CHR_DIR="${POP_OUTPUT_DIR_ROOT}/smc_estimated_model_per_chr"
INTERMEDIATES_PER_CHR_BASE_DIR="${POP_OUTPUT_DIR_ROOT}/intermediates_per_chr"

mkdir -p "${BASE_OPERATING_DIR}/slurm_logs" 
mkdir -p "${VCF_FINAL_FOR_SMC_DIR}" "${SMC_FORMAT_FILES_DIR}" \
           "${SMC_MODEL_PER_CHR_DIR}" "${INTERMEDIATES_PER_CHR_BASE_DIR}"

FULL_MASK_FILE_ABS_PATH="${BASE_OPERATING_DIR}/${FULL_MASK_FILE_RELPATH}"

if [ ! -f "${PYTHON_ANNOTATE_SCRIPT_PATH}" ]; then echo "ERROR: Annotate script ${PYTHON_ANNOTATE_SCRIPT_PATH} not found." >&2; exit 1; fi
if [ ! -f "${PYTHON_FILTER_SCRIPT_PATH}" ]; then echo "ERROR: Filter script ${PYTHON_FILTER_SCRIPT_PATH} not found." >&2; exit 1; fi
if [ ! -f "${FULL_MASK_FILE_ABS_PATH}" ]; then echo "ERROR: Full genome mask ${FULL_MASK_FILE_ABS_PATH} not found." >&2; exit 1; fi


process_single_chromosome() {
    local current_chr_num="$1"
    local current_target_chr_id="chr${current_chr_num}"
    local current_log_prefix_msg="POP: ${POP_LABEL_TO_USE}, ${current_target_chr_id}"

    log_chr() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] (${current_log_prefix_msg}) $1"; }
    log_chr "--- Starting processing ---"

    local chr_intermediates_dir="${INTERMEDIATES_PER_CHR_BASE_DIR}/chr${current_chr_num}"
    local chr_smc_model_dir="${SMC_MODEL_PER_CHR_DIR}/chr${current_chr_num}"
    mkdir -p "${chr_intermediates_dir}" "${chr_smc_model_dir}"

    local original_input_vcf="${BASE_OPERATING_DIR}/${INPUT_VCF_DIR_RELPATH}/${INPUT_VCF_PREFIX}${current_chr_num}${INPUT_VCF_SUFFIX}"
    local aa_tsv_file_original="${AA_TSV_DIR_ABS_PATH}/${AA_TSV_BASENAME_PREFIX}${current_chr_num}${AA_TSV_BASENAME_SUFFIX}"
    local current_file_prefix="${POP_LABEL_TO_USE}.chr${current_chr_num}"

    local vcf_step1_renamed_body="${chr_intermediates_dir}/${current_file_prefix}.1.renamed_body.vcf.gz"
    local header_step1_for_reheader="${chr_intermediates_dir}/${current_file_prefix}.2.header.txt"
    local vcf_step1_reheadered="${chr_intermediates_dir}/${current_file_prefix}.3.reheadered.vcf.gz"
    local vcf_step1_norm_final="${chr_intermediates_dir}/${current_file_prefix}.4.norm_final.vcf.gz"
    local temp_chrom_map_file_step1="${chr_intermediates_dir}/chrom_map_chr${current_chr_num}.txt"
    local aa_tsv_prepared_step2="${chr_intermediates_dir}/chr${current_chr_num}_AA_prepared.tsv.gz"
    local vcf_annotated_step3="${chr_intermediates_dir}/${current_file_prefix}.5.aa_annotated.vcf.gz"
    local final_vcf_for_smc_tool_input="${VCF_FINAL_FOR_SMC_DIR}/${current_file_prefix}.final_for_smc.vcf.gz"
    local mask_final_step5="${chr_intermediates_dir}/mask.${current_target_chr_id}.bed.gz"
    local smc_format_file_output_step6="${SMC_FORMAT_FILES_DIR}/${current_file_prefix}.smc.gz"
    local py_annotate_stdout_log="${chr_intermediates_dir}/py_annotate_stdout_chr${current_chr_num}.log"
    local py_annotate_stderr_log="${chr_intermediates_dir}/py_annotate_stderr_chr${current_chr_num}.log"
    local py_filter_stdout_log="${chr_intermediates_dir}/py_filter_stdout_chr${current_chr_num}.log"
    local py_filter_stderr_log="${chr_intermediates_dir}/py_filter_stderr_chr${current_chr_num}.log"

    # --- Step 0: Read Samples and Define Distinguished Individuals FROM THIS VCF ---
    log_chr "--- Step 0: Define Samples & DIs from VCF ---"
    log_chr "Input VCF for sample extraction: ${original_input_vcf}"
    if [ ! -f "${original_input_vcf}" ]; then log_chr "ERROR: Original VCF ${original_input_vcf} not found!" >&2; return 1; fi
    if [ ! -f "${original_input_vcf}.tbi" ] && [ ! -f "${original_input_vcf}.csi" ]; then
        log_chr "INFO: Index for ${original_input_vcf} not found, creating..."
        if ! tabix -f -p vcf "${original_input_vcf}"; then log_chr "ERROR: Indexing ${original_input_vcf} failed." >&2; return 1; fi
    fi
    local samples_this_vcf_list_file="${chr_intermediates_dir}/samples_from_${current_file_prefix}.list" 
    if ! bcftools query -l "${original_input_vcf}" > "${samples_this_vcf_list_file}"; then 
        log_chr "ERROR: bcftools query -l failed for ${original_input_vcf}" >&2; return 1; 
    fi
    mapfile -t SAMPLES_FROM_THIS_VCF < "${samples_this_vcf_list_file}"
    if [ ${#SAMPLES_FROM_THIS_VCF[@]} -lt 2 ]; then 
        log_chr "ERROR: Less than 2 samples found in VCF ${original_input_vcf}. Found: ${#SAMPLES_FROM_THIS_VCF[@]}" >&2; return 1; 
    fi
    local DI_1_this_vcf="${SAMPLES_FROM_THIS_VCF[0]}"
    local DI_2_this_vcf="${SAMPLES_FROM_THIS_VCF[1]}"
    log_chr "Samples in this VCF: ${#SAMPLES_FROM_THIS_VCF[@]}. DIs for this VCF: ${DI_1_this_vcf}, ${DI_2_this_vcf}"
    log_chr "--- Step 0 finished. ---"

    # Step 1: VCF Normalization
    log_chr "--- Step 1: VCF Normalization ---"
    local vcf_current_processing_step1="${original_input_vcf}" 
    # CORRECTED: Get chromosome name from first DATA line
    local original_chr_name_from_vcf_records_step1=$(zcat "${vcf_current_processing_step1}" | grep -v "^#" | head -n1 | cut -f1)
    log_chr "Original chromosome name in VCF records: '${original_chr_name_from_vcf_records_step1}'"

    if [ "${original_chr_name_from_vcf_records_step1}" != "${current_target_chr_id}" ]; then
        log_chr "Renaming VCF records: ${original_chr_name_from_vcf_records_step1} -> ${current_target_chr_id}. Output: ${vcf_step1_renamed_body}"
        echo -e "${original_chr_name_from_vcf_records_step1}\t${current_target_chr_id}" > "${temp_chrom_map_file_step1}"
        bcftools annotate --rename-chrs "${temp_chrom_map_file_step1}" "${vcf_current_processing_step1}" -O z -o "${vcf_step1_renamed_body}" --threads "${SMC_CORES}"
        if [ $? -ne 0 ] || [ ! -s "${vcf_step1_renamed_body}" ]; then log_chr "ERROR: bcftools annotate --rename-chrs failed." >&2; return 1; fi
        vcf_current_processing_step1="${vcf_step1_renamed_body}"
    fi

    bcftools view -h "${vcf_current_processing_step1}" > "${header_step1_for_reheader}.original"
    declare -A hg19_lengths_s1=( ["chr1"]="249250621" ["chr2"]="243199373" ["chr3"]="198022430" ["chr4"]="191154276" ["chr5"]="180915260" ["chr6"]="171115067" ["chr7"]="159138663" ["chrX"]="155270560" ["chr8"]="146364022" ["chr9"]="141213431" ["chr10"]="135534747" ["chr11"]="135006516" ["chr12"]="133851895" ["chr13"]="115169878" ["chr14"]="107349540" ["chr15"]="102531392" ["chr16"]="90354753"  ["chr17"]="81195210"  ["chr18"]="78077248"  ["chr19"]="59128983" ["chr20"]="63025520"  ["chr21"]="48129895"  ["chr22"]="51304566" )
    local target_len_s1="${hg19_lengths_s1[${current_target_chr_id}]}"
    local correct_contig_hdr_line_s1="##contig=<ID=${current_target_chr_id}"; if [ -n "$target_len_s1" ]; then correct_contig_hdr_line_s1+=",length=${target_len_s1}"; fi; correct_contig_hdr_line_s1+=">"
    
    local ffline_s1=$(grep "^##fileformat" "${header_step1_for_reheader}.original" || echo "##fileformat=VCFv4.2")
    echo "${ffline_s1}" > "${header_step1_for_reheader}" 
    echo "${correct_contig_hdr_line_s1}" >> "${header_step1_for_reheader}" 
    grep -E "^##(INFO|FORMAT|FILTER|ALT|SAMPLE|PEDIGREE|assembly|source|fileDate|reference)" "${header_step1_for_reheader}.original" >> "${header_step1_for_reheader}"
    grep "^#CHROM" "${header_step1_for_reheader}.original" >> "${header_step1_for_reheader}" 
    rm -f "${header_step1_for_reheader}.original"

    log_chr "Reheadering ${vcf_current_processing_step1} -> ${vcf_step1_reheadered}"
    bcftools reheader -h "${header_step1_for_reheader}" "${vcf_current_processing_step1}" -o "${vcf_step1_reheadered}"
    if [ $? -ne 0 ] || [ ! -s "${vcf_step1_reheadered}" ]; then log_chr "ERROR: bcftools reheader failed." >&2; return 1; fi
    vcf_current_processing_step1="${vcf_step1_reheadered}"

    log_chr "DEBUG: Validating reheadered file: ${vcf_current_processing_step1}"
    if ! gunzip -t "${vcf_current_processing_step1}"; then log_chr "ERROR_DEBUG: gzip integrity FAIL (post-reheader): ${vcf_current_processing_step1}"; return 1; fi
    if ! bcftools view -h "${vcf_current_processing_step1}" > /dev/null; then log_chr "ERROR_DEBUG: bcftools view -h FAIL (post-reheader): ${vcf_current_processing_step1}"; zcat "${vcf_current_processing_step1}" | head -n 10 >&2; return 1; fi
    log_chr "DEBUG: Post-reheader VCF appears valid."
    
    if [ ! -f "${vcf_current_processing_step1}.tbi" ] && [ ! -f "${vcf_current_processing_step1}.csi" ]; then
        if ! tabix -f -p vcf "${vcf_current_processing_step1}"; then log_chr "ERROR: tabix (pre-norm) failed." >&2; return 1; fi
    fi

    log_chr "Normalizing ${vcf_current_processing_step1} -> ${vcf_step1_norm_final}"
    bcftools norm -m -any --no-version "${vcf_current_processing_step1}" -O z -o "${vcf_step1_norm_final}" --threads "${SMC_CORES}"
    if [ $? -ne 0 ] || [ ! -s "${vcf_step1_norm_final}" ]; then log_chr "ERROR: bcftools norm failed." >&2; return 1; fi
    
    if ! tabix -f -p vcf "${vcf_step1_norm_final}"; then log_chr "ERROR: Indexing ${vcf_step1_norm_final} failed." >&2; return 1; fi
    local input_vcf_for_next_step="${vcf_step1_norm_final}"
    log_chr "--- Step 1 finished. Output: ${input_vcf_for_next_step} ---"

    # --- Step 2: Prepare Ancestral Allele TSV ---
    log_chr "--- Step 2: Prepare AA TSV ---"
    if [ ! -f "${aa_tsv_file_original}" ]; then log_chr "ERROR: AA TSV ${aa_tsv_file_original} not found" >&2; return 1; fi
    # CORRECTED AWK for AA TSV
    zcat "${aa_tsv_file_original}" | \
    awk -v t="${current_target_chr_id}" '
        BEGIN {OFS="\t"; print "#CHROM_AA\tPOS_AA\tAA_VALUE"} # Print our desired header first
        FNR == 1 && $1 ~ /^(#|CHR|CHROM)/ {next} # Skip first line if it looks like any kind of header
        NF >= 3 {print t, $2, $3} # Process data lines
    ' | sort -k1,1V -k2,2n | bgzip -c > "${aa_tsv_prepared_step2}"
    ps2=("${PIPESTATUS[@]}"); 
    if [[ ${ps2[0]} -ne 0 && ${ps2[0]} -ne 141 ]]; then log_chr "ERR:AA_TSV:zcat:${ps2[0]}";return 1;fi; 
    if [ ${ps2[1]} -ne 0 ];then log_chr "ERR:AA_TSV:awk:${ps2[1]}";return 1;fi; 
    if [ ${ps2[2]} -ne 0 ];then log_chr "ERR:AA_TSV:sort:${ps2[2]}";return 1;fi; 
    if [ ${ps2[3]} -ne 0 ];then log_chr "ERR:AA_TSV:bgzip:${ps2[3]}";return 1;fi
    if ! tabix -f -S 1 -s 1 -b 2 -e 2 "${aa_tsv_prepared_step2}"; then log_chr "ERROR: tabix for ${aa_tsv_prepared_step2} failed." >&2; return 1; fi
    log_chr "--- Step 2 finished. Output: ${aa_tsv_prepared_step2} ---"
    log_chr "DEBUG: First 5 lines of prepared AA TSV (${aa_tsv_prepared_step2}):"
    zcat "${aa_tsv_prepared_step2}" | head -n 5 | sed 's/\t/ [TAB] /g' 

    # --- Step 3: Annotate VCF with Ancestral Alleles (Python) ---
    # ... (Step 3 code unchanged, should now work with correct VCF chromosome names) ...
    log_chr "--- Step 3: Annotate VCF with AA (Python) ---"
    log_chr "DEBUG: First 5 data lines of VCF input to Python (${input_vcf_for_next_step}):"
    zcat "${input_vcf_for_next_step}" | grep -v "^#" | head -n 5 | cut -f 1-5 | sed 's/\t/ [TAB] /g'
    python3 "${PYTHON_ANNOTATE_SCRIPT_PATH}" "${input_vcf_for_next_step}" "${aa_tsv_prepared_step2}" "${vcf_annotated_step3}" 1>"${py_annotate_stdout_log}" 2>"${py_annotate_stderr_log}"
    py_exit_status=$?; if [ -s "${py_annotate_stderr_log}" ]; then log_chr "Python annotate stderr output:"; cat "${py_annotate_stderr_log}" >&2; fi 
    if [ ${py_exit_status} -ne 0 ]; then log_chr "ERROR: Python annotation script failed (Status: ${py_exit_status}). See above stderr." >&2; return 1; fi
    if [ ! -s "${vcf_annotated_step3}" ]; then log_chr "ERROR: Python annotation output ${vcf_annotated_step3} is empty." >&2; return 1; fi
    if ! tabix -f -p vcf "${vcf_annotated_step3}"; then log_chr "ERROR: Indexing ${vcf_annotated_step3} failed." >&2; return 1; fi
    input_vcf_for_next_step="${vcf_annotated_step3}"
    log_chr "--- Step 3 finished. Output: ${input_vcf_for_next_step} ---"

    # --- Step 4: Final VCF Filtering for SMC++ (Python script) ---
    # ... (Step 4 code unchanged) ...
    log_chr "--- Step 4: Final VCF Filtering for SMC++ (Python script) ---"
    local filter_args_s4=("${input_vcf_for_next_step}" "${final_vcf_for_smc_tool_input}")
    if [ -f "${EXCLUDE_RSID_FILE_ABS_PATH}" ]; then
        log_chr "Using RSID exclusion list: ${EXCLUDE_RSID_FILE_ABS_PATH}"
        filter_args_s4+=(--exclude_rsids "${EXCLUDE_RSID_FILE_ABS_PATH}")
    fi
    log_chr "Cmd: python3 ${PYTHON_FILTER_SCRIPT_PATH} ${filter_args_s4[*]}"
    python3 "${PYTHON_FILTER_SCRIPT_PATH}" "${filter_args_s4[@]}" 1>"${py_filter_stdout_log}" 2>"${py_filter_stderr_log}"
    py_filt_exit_status=$?; if [ -s "${py_filter_stderr_log}" ]; then log_chr "Python filter stderr:"; cat "${py_filter_stderr_log}" >&2; fi
    if [ ${py_filt_exit_status} -ne 0 ]; then log_chr "ERROR: Python filter script failed (Status: ${py_filt_exit_status})." >&2; return 1; fi
    if [ ! -s "${final_vcf_for_smc_tool_input}" ]; then 
        log_chr "ERROR: Python filter output ${final_vcf_for_smc_tool_input} is empty. This likely means no variants passed filters (e.g., all AA='.' or all excluded by rsID)." >&2; 
        return 1; 
    fi
    if ! tabix -f -p vcf "${final_vcf_for_smc_tool_input}"; then log_chr "ERROR: Indexing ${final_vcf_for_smc_tool_input} failed." >&2; return 1; fi
    input_vcf_for_next_step="${final_vcf_for_smc_tool_input}" 
    log_chr "--- Step 4 finished. Output for vcf2smc: ${input_vcf_for_next_step} ---"

    # --- Step 5: Prepare Genome Mask ---
    # ... (Step 5 code unchanged, but added more robust empty mask handling) ...
    log_chr "--- Step 5: Prepare Genome Mask ---"
    local MASK_PLAIN_TMP_S5="${chr_intermediates_dir}/mask.${current_target_chr_id}.plain.bed" 
    if [ ! -f "${FULL_MASK_FILE_ABS_PATH}" ]; then
        log_chr "WARN: Full mask file ${FULL_MASK_FILE_ABS_PATH} not found. Proceeding without mask for this chromosome."
        touch "${MASK_PLAIN_TMP_S5}" 
    else
        grep "^${current_target_chr_id}[[:space:]]" "${FULL_MASK_FILE_ABS_PATH}" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' | sort -k1,1V -k2,2n > "${MASK_PLAIN_TMP_S5}"
        ps5=("${PIPESTATUS[@]}"); 
        if [ ${ps5[0]} -gt 1 ];then log_chr "ERR:Mask:grep:${ps5[0]}";return 1;fi; 
        if [ ${ps5[1]} -ne 0 ];then log_chr "ERR:Mask:awk:${ps5[1]}";return 1;fi; 
        if [ ${ps5[2]} -ne 0 ];then log_chr "ERR:Mask:sort:${ps5[2]}";return 1;fi
        if [ ! -s "${MASK_PLAIN_TMP_S5}" ] && [ ${ps5[0]} -eq 1 ]; then 
            log_chr "WARN: No mask regions found for ${current_target_chr_id} in ${FULL_MASK_FILE_ABS_PATH}. Mask file will be empty."
        elif [ ! -s "${MASK_PLAIN_TMP_S5}" ]; then 
            log_chr "WARN: Mask file ${MASK_PLAIN_TMP_S5} is empty after processing."
        fi
    fi
    bgzip -f -c "${MASK_PLAIN_TMP_S5}" > "${mask_final_step5}"; rm -f "${MASK_PLAIN_TMP_S5}" 
    if [ -s "${mask_final_step5}" ] && [ "$(wc -c < "${mask_final_step5}")" -gt 28 ]; then # Only index if not effectively empty gzip
      if ! tabix -f -p bed "${mask_final_step5}"; then 
          log_chr "WARN: Indexing mask ${mask_final_step5} failed. This might be okay if mask is effectively empty."
      fi
    else
        log_chr "WARN: Mask file ${mask_final_step5} is empty or only gzip header. Skipping tabix."
        # Create a dummy index for empty file so smc++ doesn't complain about missing index if it checks
        touch "${mask_final_step5}.tbi" 2>/dev/null || true 
    fi
    log_chr "--- Step 5 finished. Output: ${mask_final_step5} ---"
    
    # --- Step 6: Run smc++ vcf2smc ---
    # ... (Step 6 code unchanged) ...
    log_chr "--- Step 6: smc++ vcf2smc ---"
    local SMC_SAMPLE_LIST_STR_S6="${POP_LABEL_TO_USE}:$(IFS=,; echo "${SAMPLES_FROM_THIS_VCF[*]}")" 
    log_chr "Input VCF for vcf2smc: ${input_vcf_for_next_step}"
    num_variants_for_smc=$(zcat "${input_vcf_for_next_step}" | grep -vc "^#")
    if [ "${num_variants_for_smc}" -eq 0 ]; then
        log_chr "ERROR: VCF file ${input_vcf_for_next_step} has no variants. Cannot run smc++ vcf2smc." >&2
        return 1
    fi

    smc++ vcf2smc --cores "${SMC_CORES}" -m "${mask_final_step5}" -d "${DI_1_this_vcf}" "${DI_2_this_vcf}" \
        "${input_vcf_for_next_step}" \
        "${smc_format_file_output_step6}" \
        "${current_target_chr_id}" \
        "${SMC_SAMPLE_LIST_STR_S6}"
    if [ $? -ne 0 ] || [ ! -s "${smc_format_file_output_step6}" ]; then log_chr "ERROR: smc++ vcf2smc failed or produced empty ${smc_format_file_output_step6}." >&2; return 1; fi
    log_chr "--- Step 6 finished. Output: ${smc_format_file_output_step6} ---"

    # --- Step 7: Run smc++ estimate ---
    # ... (Step 7 code unchanged) ...
    log_chr "--- Step 7: smc++ estimate ---"
    smc++ estimate --cores "${SMC_CORES}" -o "${chr_smc_model_dir}" \
        "${UNIFORM_RECOMBINATION_RATE}" \
        "${smc_format_file_output_step6}"
    if [ $? -ne 0 ] || [ ! -f "${chr_smc_model_dir}/model.final.json" ]; then log_chr "ERROR: smc++ estimate failed or model.final.json not found." >&2; return 1; fi
    log_chr "--- Step 7 finished. Model in ${chr_smc_model_dir} ---"

    # Cleanup
    log_chr "--- Cleaning up temporary files ---"
    # Add more files to cleanup if necessary from Step 1
    rm -f "${header_step1_for_reheader}" "${header_step1_for_reheader}.new" "${header_step1_for_reheader}.bak"* \
          "${temp_chrom_map_file_step1}" \
          "${vcf_step1_renamed_body}" "${vcf_step1_renamed_body}.tbi" "${vcf_step1_renamed_body}.csi" \
          "${vcf_step1_reheadered}" "${vcf_step1_reheadered}.tbi" "${vcf_step1_reheadered}.csi" \
          "${SAMPLES_TEMP_LIST}" \
          "${samples_this_vcf_list_file}" \
          2>/dev/null || log_chr "Note: Some temp files already cleaned or not created."
    log_chr "--- SUCCESS ---"
    return 0
}

# --- Main Script Logic ---
MAIN_LOG_PREFIX_SCRIPT="POP: ${POP_LABEL_TO_USE}, MainScript"
log_main_script() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] (${MAIN_LOG_PREFIX_SCRIPT}@$(hostname)) $1"; }

log_main_script "--- Initializing ---"
# ... (rest of main script logic is the same as previous version) ...
log_main_script "RUN_MODE: ${RUN_MODE}"
log_main_script "BASE_OPERATING_DIR (PWD of launch): ${BASE_OPERATING_DIR}"
log_main_script "SMC_CORES: ${SMC_CORES}"
log_main_script "Script actual location: ${SCRIPT_ACTUAL_LOCATION}"
log_main_script "Python annotate script: ${PYTHON_ANNOTATE_SCRIPT_PATH}"
log_main_script "Python filter script: ${PYTHON_FILTER_SCRIPT_PATH}"

CHROMOSOME_LIST_TO_PROCESS=()
if [ -n "${SPECIFIC_CHR_NUM_ARG_FROM_CMD}" ]; then 
    CHROMOSOME_LIST_TO_PROCESS+=("${SPECIFIC_CHR_NUM_ARG_FROM_CMD}")
    log_main_script "Processing specific chromosome(s): ${CHROMOSOME_LIST_TO_PROCESS[*]}"
else
    for i in $(seq 1 22); do CHROMOSOME_LIST_TO_PROCESS+=("$i"); done
    log_main_script "Processing all autosomes (1-22)."
fi

OVERALL_JOB_STATUS=0
SUCCESSFUL_CHROMS=0
FAILED_CHROMS=0

for chr_num_to_process in "${CHROMOSOME_LIST_TO_PROCESS[@]}"; do
    echo 
    if process_single_chromosome "${chr_num_to_process}"; then
        SUCCESSFUL_CHROMS=$((SUCCESSFUL_CHROMS + 1))
    else
        log_main_script "ERROR: Processing FAILED for chr${chr_num_to_process}." >&2
        FAILED_CHROMS=$((FAILED_CHROMS + 1))
        OVERALL_JOB_STATUS=1 
    fi
done

log_main_script "Finished processing designated chromosomes."
log_main_script "Summary: SUCCESSFUL = ${SUCCESSFUL_CHROMS}, FAILED = ${FAILED_CHROMS}"
if [ ${OVERALL_JOB_STATUS} -ne 0 ]; then
    log_main_script "--- SCRIPT FINISHED WITH ERRORS ---" >&2
fi
exit ${OVERALL_JOB_STATUS}
