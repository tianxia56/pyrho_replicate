#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status
set -o pipefail # Causes a pipeline to return the exit status of the last command in the pipe that failed

# --- Configuration ---
TARGET_POP_LABEL="IND"
ORIGINAL_INPUT_VCF_RELPATH="${TARGET_POP_LABEL}/${TARGET_POP_LABEL}.chr22.vcf.gz" 
MASTER_SAMPLES_FILE_ABS_PATH="/home/tx56/palmer_scratch/100kga/100kga/sample_lists/${TARGET_POP_LABEL}.IID.list"
EXCLUDE_RSID_FILE_ABS_PATH="/home/tx56/1KGP/exclude_rsid.txt" 

AA_TSV_DIR_ABS_PATH="/vast/palmer/pi/reilly/jfa38/datasets_for_annotation/ErinG_ARG_hg19_AncestralAlleles/tidied_AA_tsvs"
FULL_MASK_FILE_RELPATH="20141020.strict_mask.whole_genome.bed" 
TARGET_CHROMOSOME_ID="chr22" 
UNIFORM_RECOMBINATION_RATE="1.25e-8"
SMC_CORES=4
REFERENCE_FASTA="/path/to/your/hg19.fa" 
PYTHON_ANNOTATE_SCRIPT_NAME="annotate_aa.py" 
PYTHON_FILTER_SCRIPT_NAME="filter_vcf_for_smc.py" 
# --- End Configuration ---

TARGET_CHROMOSOME_NUMERIC_PART=$(echo "${TARGET_CHROMOSOME_ID}" | sed 's/[^0-9]*//g') 
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
PYTHON_ANNOTATE_SCRIPT_PATH="${SCRIPT_DIR}/${PYTHON_ANNOTATE_SCRIPT_NAME}"
PYTHON_FILTER_SCRIPT_PATH="${SCRIPT_DIR}/${PYTHON_FILTER_SCRIPT_NAME}" 

MAIN_OUTPUT_DIR="${SCRIPT_DIR}/smc_processing_${TARGET_POP_LABEL}_chr${TARGET_CHROMOSOME_NUMERIC_PART}"
VCF2SMC_INTERMEDIATE_DIR="${MAIN_OUTPUT_DIR}/vcf2smc_intermediate_files"
SMC_ESTIMATE_INITIAL_DIR="${MAIN_OUTPUT_DIR}/smc_estimate_initial_model"

mkdir -p "${VCF2SMC_INTERMEDIATE_DIR}"
mkdir -p "${SMC_ESTIMATE_INITIAL_DIR}"

ORIGINAL_INPUT_VCF_ABS_PATH="${SCRIPT_DIR}/${ORIGINAL_INPUT_VCF_RELPATH}" 
AA_TSV_ORIGINAL_BASENAME="chr${TARGET_CHROMOSOME_NUMERIC_PART}_ErinG_AA_hg19.tsv.gz"
AA_TSV_ORIGINAL="${AA_TSV_DIR_ABS_PATH}/${AA_TSV_ORIGINAL_BASENAME}"
FULL_MASK_FILE_ABS_PATH="${SCRIPT_DIR}/${FULL_MASK_FILE_RELPATH}"

FILE_PREFIX="${TARGET_POP_LABEL}.chr${TARGET_CHROMOSOME_NUMERIC_PART}"
VCF_BODY_RENAMED="${VCF2SMC_INTERMEDIATE_DIR}/${FILE_PREFIX}.phased.bodyRenamed.vcf.gz"
TEMP_RENAMED_HEADER_FILE="${VCF2SMC_INTERMEDIATE_DIR}/temp_renamed_header.txt"
TEMP_CHROM_MAP_FILE="${VCF2SMC_INTERMEDIATE_DIR}/temp_chromosome_map.txt"
VCF_CHR_NORMALIZED_RAW="${VCF2SMC_INTERMEDIATE_DIR}/${FILE_PREFIX}.phased.chrNormalized.raw.vcf.gz"
VCF_NORM_SPLIT_TEMP="${VCF2SMC_INTERMEDIATE_DIR}/${FILE_PREFIX}.phased.chrNormalized.norm_split_temp.vcf.gz" # New temp file
VCF_CHR_NORMALIZED_CLEANED="${VCF2SMC_INTERMEDIATE_DIR}/${FILE_PREFIX}.phased.chrNormalized.cleaned.vcf.gz"
AA_TSV_PREPARED="${VCF2SMC_INTERMEDIATE_DIR}/chr${TARGET_CHROMOSOME_NUMERIC_PART}_ErinG_AA_hg19_prepared.tsv.gz"
VCF_PYTHON_ANNOTATED="${VCF2SMC_INTERMEDIATE_DIR}/${FILE_PREFIX}.phased.annotated.py.vcf.gz" 
VCF_FINAL_FILTERED="${VCF2SMC_INTERMEDIATE_DIR}/${FILE_PREFIX}.phased.annotated.final_filtered.vcf.gz"
MASK_FILE_FINAL="${VCF2SMC_INTERMEDIATE_DIR}/mask.${TARGET_CHROMOSOME_ID}.bed.gz" 
SMC_FORMAT_FILE="${VCF2SMC_INTERMEDIATE_DIR}/${FILE_PREFIX}.smc.gz"

VCF_SAMPLES_TEMP="${VCF2SMC_INTERMEDIATE_DIR}/vcf_samples.tmp.txt"
MASTER_SAMPLES_SORTED_TEMP="${VCF2SMC_INTERMEDIATE_DIR}/master_samples_sorted.tmp.txt"
COMMON_SAMPLES_LIST_FILE="${VCF2SMC_INTERMEDIATE_DIR}/${TARGET_POP_LABEL}.common_samples.list"

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] (${TARGET_POP_LABEL}-${TARGET_CHROMOSOME_ID}) $1"
}

log "--- Sanity Checks ---"
for tool in bcftools tabix bgzip smc++ samtools python3 comm sort awk sed grep; do 
    if ! command -v $tool &> /dev/null; then
        log "ERROR: Required tool '$tool' not found. Please ensure it is installed and in your PATH."
        exit 1
    fi
done
log "All required command-line tools found."
if [ ! -f "${PYTHON_ANNOTATE_SCRIPT_PATH}" ]; then
    log "ERROR: Python annotation script ${PYTHON_ANNOTATE_SCRIPT_PATH} not found."
    exit 1
fi
log "Python annotation script found at ${PYTHON_ANNOTATE_SCRIPT_PATH}."
if [ ! -f "${PYTHON_FILTER_SCRIPT_PATH}" ]; then
    log "ERROR: Python filtering script ${PYTHON_FILTER_SCRIPT_PATH} not found."
    exit 1
fi
log "Python filtering script found at ${PYTHON_FILTER_SCRIPT_PATH}."
EXCLUDE_RSID_FILE_TO_USE="" 
if [ ! -f "${EXCLUDE_RSID_FILE_ABS_PATH}" ]; then
    log "WARNING: RSID exclusion file ${EXCLUDE_RSID_FILE_ABS_PATH} not found. Proceeding without excluding rsIDs from this list."
else
    log "RSID exclusion file found at ${EXCLUDE_RSID_FILE_ABS_PATH}."
    EXCLUDE_RSID_FILE_TO_USE="${EXCLUDE_RSID_FILE_ABS_PATH}"
fi
log "--- End Sanity Checks ---"

log "--- Starting Step 0: Determine Usable Samples and Define Distinguished Individuals ---"
# ... (Step 0 - exactly as in your last provided log - it was working) ...
if [ ! -f "${ORIGINAL_INPUT_VCF_ABS_PATH}" ]; then
    log "ERROR: Original VCF not found at Step 0: ${ORIGINAL_INPUT_VCF_ABS_PATH}"
    exit 1
fi
log "Extracting sample IDs from VCF header: ${ORIGINAL_INPUT_VCF_ABS_PATH}"
zcat "${ORIGINAL_INPUT_VCF_ABS_PATH}" | grep "^#CHROM" | head -n 1 | cut -f 10- | tr '\t' '\n' | sort > "${VCF_SAMPLES_TEMP}"
if [ ! -s "${VCF_SAMPLES_TEMP}" ]; then
    log "ERROR: No sample IDs extracted from VCF header, or VCF is malformed: ${ORIGINAL_INPUT_VCF_ABS_PATH}"
    exit 1
fi
log "VCF contains $(wc -l < "${VCF_SAMPLES_TEMP}") sample IDs."
if [ ! -f "${MASTER_SAMPLES_FILE_ABS_PATH}" ]; then
    log "ERROR: Master population sample list not found: ${MASTER_SAMPLES_FILE_ABS_PATH}"
    exit 1
fi
log "Sorting master sample list: ${MASTER_SAMPLES_FILE_ABS_PATH}"
sort "${MASTER_SAMPLES_FILE_ABS_PATH}" > "${MASTER_SAMPLES_SORTED_TEMP}"
log "Master list contains $(wc -l < "${MASTER_SAMPLES_SORTED_TEMP}") sample IDs."
log "Finding common samples between VCF and master list..."
comm -12 "${VCF_SAMPLES_TEMP}" "${MASTER_SAMPLES_SORTED_TEMP}" > "${COMMON_SAMPLES_LIST_FILE}"
if [ ! -s "${COMMON_SAMPLES_LIST_FILE}" ]; then
    log "ERROR: No common samples found between the VCF file and the master sample list (${MASTER_SAMPLES_FILE_ABS_PATH})."
    exit 1
fi
mapfile -t SAMPLES_FOR_SMC < <(cat "${COMMON_SAMPLES_LIST_FILE}")
if [ ${#SAMPLES_FOR_SMC[@]} -lt 2 ]; then
    log "ERROR: Not enough common samples found (need at least 2 for DIs) to proceed. Found: ${#SAMPLES_FOR_SMC[@]}"
    exit 1
fi
DISTINGUISHED_INDIVIDUALS_TARGET_POP=("${SAMPLES_FOR_SMC[0]}" "${SAMPLES_FOR_SMC[1]}")
log "Total common samples to be used for ${TARGET_POP_LABEL} in this VCF: ${#SAMPLES_FOR_SMC[@]}"
log "Distinguished Individuals for ${TARGET_POP_LABEL} (from common set): ${DISTINGUISHED_INDIVIDUALS_TARGET_POP[0]}, ${DISTINGUISHED_INDIVIDUALS_TARGET_POP[1]}"
log "List of common samples being used is in: ${COMMON_SAMPLES_LIST_FILE}"
log "--- Step 0 finished. ---"

log "--- Starting Step 1: Normalize Input VCF Chromosome Names & Clean ---"
# ... (Step 1a and 1b - exactly as in your last provided log - they were working) ...
log "Path to original VCF: '${ORIGINAL_INPUT_VCF_ABS_PATH}'"
ORIGINAL_INPUT_VCF_INDEX_TBI="${ORIGINAL_INPUT_VCF_ABS_PATH}.tbi"
ORIGINAL_INPUT_VCF_INDEX_CSI="${ORIGINAL_INPUT_VCF_ABS_PATH}.csi"
if [ ! -f "${ORIGINAL_INPUT_VCF_INDEX_TBI}" ] && [ ! -f "${ORIGINAL_INPUT_VCF_INDEX_CSI}" ]; then
    log "INFO: Index for original VCF (${ORIGINAL_INPUT_VCF_ABS_PATH}) not found (.tbi or .csi)."
    log "Attempting to index ${ORIGINAL_INPUT_VCF_ABS_PATH} with tabix..."
    if tabix -f -p vcf "${ORIGINAL_INPUT_VCF_ABS_PATH}"; then
        log "Successfully indexed ${ORIGINAL_INPUT_VCF_ABS_PATH}."
         if [ ! -f "${ORIGINAL_INPUT_VCF_INDEX_TBI}" ] && [ ! -f "${ORIGINAL_INPUT_VCF_INDEX_CSI}" ]; then
            log "ERROR: Indexing attempt seemed to succeed but index file still not found."
            exit 1
         fi
    else
        log "ERROR: Failed to index ${ORIGINAL_INPUT_VCF_ABS_PATH}."
        exit 1
    fi
fi
log "Original VCF and its index seem OK."
original_vcf_chrom_name_sample=""
_command_output=""
_exit_status=0
set +e 
_command_output=$(zcat "${ORIGINAL_INPUT_VCF_ABS_PATH}" 2>/dev/null | grep -v "^#" | head -n 1 | awk '{print $1}')
_exit_status=$?
set -e 
if [ ${_exit_status} -eq 0 ] && [ -n "$_command_output" ]; then
    original_vcf_chrom_name_sample="$_command_output"
    log "Cmd to get chrom name sample succeeded. Name: '${original_vcf_chrom_name_sample}'"
elif [ ${_exit_status} -ne 0 ] && [ -n "$_command_output" ]; then 
    original_vcf_chrom_name_sample="$_command_output"
    log "Cmd to get chrom name sample had non-zero exit (${_exit_status}) but produced output. Name: '${original_vcf_chrom_name_sample}'. Proceeding cautiously."
else
    log "ERROR: Cmd to get chrom name sample failed. Exit: ${_exit_status}. Output: '$_command_output'."
    zcat "${ORIGINAL_INPUT_VCF_ABS_PATH}" | head -n 20
    exit 1
fi
if [ -z "${original_vcf_chrom_name_sample}" ]; then
    log "ERROR: Could not determine chromosome name sample from ${ORIGINAL_INPUT_VCF_ABS_PATH}."
    exit 1
fi
log "Successfully determined original VCF chromosome name sample: '${original_vcf_chrom_name_sample}'"
known_length_for_target_chr=""
if [ "${TARGET_CHROMOSOME_ID}" == "chr22" ]; then known_length_for_target_chr="51304566"; fi
if [ "${original_vcf_chrom_name_sample}" != "${TARGET_CHROMOSOME_ID}" ]; then
    log "Normalization required: VCF has '${original_vcf_chrom_name_sample}', target is '${TARGET_CHROMOSOME_ID}'."
    echo -e "${original_vcf_chrom_name_sample}\t${TARGET_CHROMOSOME_ID}" > "${TEMP_CHROM_MAP_FILE}"
    log "Step 1a: Renaming chromosome names in VCF records..."
    bcftools annotate --rename-chrs "${TEMP_CHROM_MAP_FILE}" "${ORIGINAL_INPUT_VCF_ABS_PATH}" -O z -o "${VCF_BODY_RENAMED}"
    log "Step 1b: Fixing header of body-renamed VCF to include correct ID and length..."
    bcftools view -h "${VCF_BODY_RENAMED}" > "${TEMP_RENAMED_HEADER_FILE}"
    correct_contig_line="##contig=<ID=${TARGET_CHROMOSOME_ID}"
    if [ -n "${known_length_for_target_chr}" ]; then
        correct_contig_line="${correct_contig_line},length=${known_length_for_target_chr}"
    fi
    correct_contig_line="${correct_contig_line}>"
    awk -v orig_id="${original_vcf_chrom_name_sample}" \
        -v target_id="${TARGET_CHROMOSOME_ID}" \
        -v replacement_line="${correct_contig_line}" \
        'BEGIN{FS=OFS="\t"; replaced=0}
         $0 ~ ("^##contig=<ID=" orig_id "[,>]") || (orig_id != target_id && $0 ~ ("^##contig=<ID=" target_id "[,>]")) {
             if (!replaced) {
                 print replacement_line;
                 replaced=1;
             }
             next;
         }
         { print }' \
        "${TEMP_RENAMED_HEADER_FILE}" > "${TEMP_RENAMED_HEADER_FILE}.final"
    mv "${TEMP_RENAMED_HEADER_FILE}.final" "${TEMP_RENAMED_HEADER_FILE}"
    if ! grep -Fxq "${correct_contig_line}" "${TEMP_RENAMED_HEADER_FILE}"; then
        log "Correct contig line ('${correct_contig_line}') was not found after awk processing. Adding/Fixing."
        sed -i".bak" "/^##contig=<ID=${TARGET_CHROMOSOME_ID}/d" "${TEMP_RENAMED_HEADER_FILE}" 
        echo "${correct_contig_line}" >> "${TEMP_RENAMED_HEADER_FILE}"
        log "Added/Ensured contig definition: ${correct_contig_line}"
    else
        log "Correct contig line found in header: ${correct_contig_line}"
    fi
    log "Reheadering body-renamed VCF with corrected header to ${VCF_CHR_NORMALIZED_RAW}..."
    bcftools reheader -h "${TEMP_RENAMED_HEADER_FILE}" "${VCF_BODY_RENAMED}" -o "${VCF_CHR_NORMALIZED_RAW}"
else
    log "Chromosome names in original VCF ('${original_vcf_chrom_name_sample}') already match target ('${TARGET_CHROMOSOME_ID}'). Checking/Ensuring contig length."
    bcftools view -h "${ORIGINAL_INPUT_VCF_ABS_PATH}" > "${TEMP_RENAMED_HEADER_FILE}"
    needs_reheader=false
    correct_contig_line_check="##contig=<ID=${TARGET_CHROMOSOME_ID}"
    if [ -n "${known_length_for_target_chr}" ]; then
        correct_contig_line_check="${correct_contig_line_check},length=${known_length_for_target_chr}"
    fi
    correct_contig_line_check="${correct_contig_line_check}>"
    if ! grep -Fxq "${correct_contig_line_check}" "${TEMP_RENAMED_HEADER_FILE}"; then
        log "Correct contig line ('${correct_contig_line_check}') not found. Attempting to fix/add."
        sed -i".bak" "/^##contig=<ID=${TARGET_CHROMOSOME_ID}/d" "${TEMP_RENAMED_HEADER_FILE}"
        echo "${correct_contig_line_check}" >> "${TEMP_RENAMED_HEADER_FILE}"
        needs_reheader=true
        log "Ensured contig definition: ${correct_contig_line_check}"
    else
        log "Correct contig line with length already present."
    fi
    if [ "$needs_reheader" = true ]; then
        log "Reheadering original VCF to ensure correct contig length for ${TARGET_CHROMOSOME_ID}..."
        bcftools reheader -h "${TEMP_RENAMED_HEADER_FILE}" "${ORIGINAL_INPUT_VCF_ABS_PATH}" -o "${VCF_CHR_NORMALIZED_RAW}"
    else
        log "Contig definition appears correct. Copying to ${VCF_CHR_NORMALIZED_RAW}."
        cp "${ORIGINAL_INPUT_VCF_ABS_PATH}" "${VCF_CHR_NORMALIZED_RAW}"
    fi
fi
rm -f "${TEMP_RENAMED_HEADER_FILE}.bak"

# --- Step 1c: bcftools norm ---
log "Step 1c: Cleaning the VCF (${VCF_CHR_NORMALIZED_RAW}) using bcftools norm..."
if [ ! -f "${REFERENCE_FASTA}" ] || [ "${REFERENCE_FASTA}" == "/path/to/your/hg19.fa" ]; then 
    log "WARNING: Reference FASTA ${REFERENCE_FASTA} not found or default path used. Running bcftools norm without -f."
    log "Normalizing VCF: Splitting multiallelics..."
    bcftools norm -m -any "${VCF_CHR_NORMALIZED_RAW}" -O z -o "${VCF_NORM_SPLIT_TEMP}" --threads "${SMC_CORES}"
    if [ $? -ne 0 ]; then log "ERROR: bcftools norm (splitting multiallelics) failed."; exit 1; fi
    
    log "Normalizing VCF: Removing exact duplicates from split VCF..."
    bcftools norm -d exact "${VCF_NORM_SPLIT_TEMP}" -O z -o "${VCF_CHR_NORMALIZED_CLEANED}" --threads "${SMC_CORES}"
    if [ $? -ne 0 ]; then log "ERROR: bcftools norm (removing exact duplicates) failed."; exit 1; fi
else
    log "Using reference FASTA for normalization: ${REFERENCE_FASTA}"
    if [ ! -f "${REFERENCE_FASTA}.fai" ]; then
        log "INFO: FASTA index ${REFERENCE_FASTA}.fai not found. Creating index..."
        samtools faidx "${REFERENCE_FASTA}"
        if [ $? -ne 0 ] || [ ! -f "${REFERENCE_FASTA}.fai" ]; then
            log "ERROR: Failed to create FASTA index for ${REFERENCE_FASTA}."
            exit 1
        fi
        log "FASTA index created successfully."
    fi
    log "Normalizing VCF: Splitting multiallelics and left-aligning indels..."
    bcftools norm -m -any -f "${REFERENCE_FASTA}" "${VCF_CHR_NORMALIZED_RAW}" -O z -o "${VCF_NORM_SPLIT_TEMP}" --threads "${SMC_CORES}"
    if [ $? -ne 0 ]; then log "ERROR: bcftools norm (splitting multiallelics with -f) failed."; exit 1; fi

    log "Normalizing VCF: Removing exact duplicates from split and aligned VCF..."
    bcftools norm -d exact -f "${REFERENCE_FASTA}" "${VCF_NORM_SPLIT_TEMP}" -O z -o "${VCF_CHR_NORMALIZED_CLEANED}" --threads "${SMC_CORES}"
    if [ $? -ne 0 ]; then log "ERROR: bcftools norm (removing exact duplicates with -f) failed."; exit 1; fi
fi
rm -f "${VCF_NORM_SPLIT_TEMP}" "${VCF_NORM_SPLIT_TEMP}.tbi" 2>/dev/null || true 
log "bcftools norm operations completed."
log "Indexing cleaned VCF: ${VCF_CHR_NORMALIZED_CLEANED}"
tabix -f -p vcf "${VCF_CHR_NORMALIZED_CLEANED}"
INPUT_VCF_FOR_AA_ANNOTATION="${VCF_CHR_NORMALIZED_CLEANED}" 
log "--- Step 1 (Chrom Name Normalization & Cleaning) finished. VCF for AA annotation: ${INPUT_VCF_FOR_AA_ANNOTATION} ---"

# --- Step 2: Prepare AA TSV ---
# ... (Step 2 logic - exactly as in your last provided log - it was working) ...
log "--- Starting Step 2: Prepare Ancestral Allele TSV ---"
if [ ! -f "${AA_TSV_ORIGINAL}" ]; then
    log "ERROR: Original AA TSV file not found: ${AA_TSV_ORIGINAL}"
    exit 1
fi
log "DEBUG_STEP2: First 5 lines of ORIGINAL AA TSV (${AA_TSV_ORIGINAL}):"
zcat "${AA_TSV_ORIGINAL}" | head -n 5 || log "DEBUG_STEP2_WARNING: Failed to read head of ${AA_TSV_ORIGINAL}"
log "DEBUG_STEP2: Line count of ORIGINAL AA TSV (${AA_TSV_ORIGINAL}):"
original_aa_line_count=$(zcat "${AA_TSV_ORIGINAL}" | wc -l) || original_aa_line_count="Error"
log "DEBUG_STEP2: Original AA TSV line count: ${original_aa_line_count}"
log "Preparing AA TSV: ${AA_TSV_PREPARED} from ${AA_TSV_ORIGINAL}"
zcat "${AA_TSV_ORIGINAL}" | \
    awk -v target_chrom_awk="${TARGET_CHROMOSOME_ID}" '
        BEGIN {
            OFS="\t";
            header_printed=0; 
        }
        ($0 ~ /^#/) { 
            if (!header_printed && NR==1) { 
                print $0; 
                header_printed=1;
            }
            next 
        }
        (!header_printed) {
            print "#CHROM_AA", "POS_AA", "AA_VALUE";
            header_printed=1;
        }
        { if (NF >=3) print target_chrom_awk, $2, $3 }
    ' | \
    sort -k1,1V -k2,2n | \
    bgzip -c > "${AA_TSV_PREPARED}"
pipe_status_step2=("${PIPESTATUS[@]}")
if [ ${pipe_status_step2[0]} -ne 0 ] || [ ${pipe_status_step2[1]} -ne 0 ] || [ ${pipe_status_step2[2]} -ne 0 ] || [ ${pipe_status_step2[3]} -ne 0 ]; then
    log "ERROR: A command in the AA TSV preparation pipeline failed."
    log "PIPESTATUS for (zcat | awk | sort | bgzip): ${pipe_status_step2[*]}"
    exit 1
fi
log "AA TSV preparation pipeline completed."
log "DEBUG_STEP2: First 5 lines of PREPARED AA TSV (${AA_TSV_PREPARED}) after creation:"
zcat "${AA_TSV_PREPARED}" | head -n 5 || log "DEBUG_STEP2_WARNING: Failed to read head of ${AA_TSV_PREPARED} after creation"
log "DEBUG_STEP2: Line count of PREPARED AA TSV (${AA_TSV_PREPARED}) after creation:"
prepared_aa_line_count=$(zcat "${AA_TSV_PREPARED}" | wc -l) || prepared_aa_line_count="Error"
log "DEBUG_STEP2: Prepared AA TSV line count: ${prepared_aa_line_count}"
log "Indexing prepared AA TSV: ${AA_TSV_PREPARED}"
tabix -f -S 1 -s 1 -b 2 -e 2 "${AA_TSV_PREPARED}" 
log "--- Step 2 (Prepare AA TSV) finished. Prepared AA TSV is: ${AA_TSV_PREPARED} ---"

# --- Step 3: AA Annotation (Python) ---
# ... (Step 3 logic - exactly as in your last provided log - it was working) ...
log "--- Starting Step 3: Annotate VCF with Ancestral Alleles using Python ---"
log "Input VCF for Python annotation: ${INPUT_VCF_FOR_AA_ANNOTATION}"
log "AA TSV for Python annotation: ${AA_TSV_PREPARED}"
log "Output VCF from Python annotation: ${VCF_PYTHON_ANNOTATED}"
PYTHON_STDOUT_LOG="${VCF2SMC_INTERMEDIATE_DIR}/python_annotate_${FILE_PREFIX}_stdout.log"
PYTHON_STDERR_LOG="${VCF2SMC_INTERMEDIATE_DIR}/python_annotate_${FILE_PREFIX}_stderr.log"
log "Python script stdout will be logged to: ${PYTHON_STDOUT_LOG}"
log "Python script stderr (diagnostics) will be logged to: ${PYTHON_STDERR_LOG}"
python3 "${PYTHON_ANNOTATE_SCRIPT_PATH}" \
    "${INPUT_VCF_FOR_AA_ANNOTATION}" \
    "${AA_TSV_PREPARED}" \
    "${VCF_PYTHON_ANNOTATED}" \
    1> "${PYTHON_STDOUT_LOG}" \
    2> "${PYTHON_STDERR_LOG}"
PYTHON_EXIT_STATUS=$?
log "--- Content of Python stderr log (${PYTHON_STDERR_LOG}) ---"
if [ -f "${PYTHON_STDERR_LOG}" ]; then
    cat "${PYTHON_STDERR_LOG}" 
else
    log "WARNING: Python stderr log file not found: ${PYTHON_STDERR_LOG}"
fi
log "--- End of Python stderr log content ---"
if [ -s "${PYTHON_STDOUT_LOG}" ]; then 
    log "--- Content of Python stdout log (${PYTHON_STDOUT_LOG}) ---"
    cat "${PYTHON_STDOUT_LOG}"
    log "--- End of Python stdout log content ---"
fi
if [ ${PYTHON_EXIT_STATUS} -ne 0 ]; then
    log "ERROR: Python annotation script failed with exit status ${PYTHON_EXIT_STATUS}."
    exit 1
fi
log "Python annotation script appears to have finished successfully (exit status 0)."
VCF_TO_BE_FILTERED="${VCF_PYTHON_ANNOTATED}" 

# --- Step 4: Final VCF Filtering using Python Script ---
log "--- Starting Step 4: Final VCF Filtering ---"
log "Input VCF for final filtering: ${VCF_TO_BE_FILTERED}"
log "Output VCF after final filtering: ${VCF_FINAL_FILTERED}"
log "RSID exclusion list for filtering: ${EXCLUDE_RSID_FILE_TO_USE}" 
PYTHON_FINAL_FILTER_STDOUT_LOG="${VCF2SMC_INTERMEDIATE_DIR}/python_final_filter_${FILE_PREFIX}_stdout.log"
PYTHON_FINAL_FILTER_STDERR_LOG="${VCF2SMC_INTERMEDIATE_DIR}/python_final_filter_${FILE_PREFIX}_stderr.log"
log "Python final filter script stdout will be logged to: ${PYTHON_FINAL_FILTER_STDOUT_LOG}"
log "Python final filter script stderr (diagnostics) will be logged to: ${PYTHON_FINAL_FILTER_STDERR_LOG}"
PYTHON_FILTER_ARGS=("${VCF_TO_BE_FILTERED}" "${VCF_FINAL_FILTERED}")
if [ -n "${EXCLUDE_RSID_FILE_TO_USE}" ] && [ -f "${EXCLUDE_RSID_FILE_TO_USE}" ]; then
    PYTHON_FILTER_ARGS+=(--exclude_rsids "${EXCLUDE_RSID_FILE_TO_USE}")
    log "Passing --exclude_rsids ${EXCLUDE_RSID_FILE_TO_USE} to filter script."
else
    log "Not passing --exclude_rsids to filter script (no valid file or file path empty)."
fi
python3 "${PYTHON_FILTER_SCRIPT_PATH}" \
    "${PYTHON_FILTER_ARGS[@]}" \
    1> "${PYTHON_FINAL_FILTER_STDOUT_LOG}" \
    2> "${PYTHON_FINAL_FILTER_STDERR_LOG}"
PYTHON_FILTER_EXIT_STATUS=$?
log "--- Content of Python final filter stderr log (${PYTHON_FINAL_FILTER_STDERR_LOG}) ---"
if [ -f "${PYTHON_FINAL_FILTER_STDERR_LOG}" ]; then
    cat "${PYTHON_FINAL_FILTER_STDERR_LOG}" 
else
    log "WARNING: Python final filter stderr log file not found."
fi
log "--- End of Python final filter stderr log content ---"
if [ -s "${PYTHON_FINAL_FILTER_STDOUT_LOG}" ]; then 
    log "--- Content of Python final filter stdout log (${PYTHON_FINAL_FILTER_STDOUT_LOG}) ---"
    cat "${PYTHON_FINAL_FILTER_STDOUT_LOG}"
    log "--- End of Python final filter stdout log content ---"
fi
if [ ${PYTHON_FILTER_EXIT_STATUS} -ne 0 ]; then
    log "ERROR: Python final filtering script failed with exit status ${PYTHON_FILTER_EXIT_STATUS}."
    exit 1
fi
if [ ! -s "${VCF_FINAL_FILTERED}" ]; then
    log "ERROR: Python final filtering script completed but output file ${VCF_FINAL_FILTERED} is empty or does not exist."
    log "This might happen if ALL variants were filtered out. Check the filter script's log."
    exit 1
fi
log "Python final filtering script completed successfully."
log "Indexing final filtered VCF: ${VCF_FINAL_FILTERED}"
tabix -f -p vcf "${VCF_FINAL_FILTERED}"
INPUT_VCF_CURRENT_STEP="${VCF_FINAL_FILTERED}"
log "--- Step 4 (Final VCF Filtering) finished. VCF for smc++ is now: ${INPUT_VCF_CURRENT_STEP} ---"

# --- Step 5: Prepare Genome Mask ---
# ... (Step 5 logic - exactly as in your last provided log - it was working) ...
log "--- Starting Step 5: Prepare Genome Mask ---"
MASK_FILE_PLAIN_TEMP="${VCF2SMC_INTERMEDIATE_DIR}/mask.${TARGET_CHROMOSOME_ID}.temp.plain.bed"
if [ ! -f "${FULL_MASK_FILE_ABS_PATH}" ]; then
    log "ERROR: Full genome mask file not found: ${FULL_MASK_FILE_ABS_PATH}"
    exit 1
fi
log "Extracting and sorting mask regions for ${TARGET_CHROMOSOME_ID} from ${FULL_MASK_FILE_ABS_PATH}..."
grep "^${TARGET_CHROMOSOME_ID}[[:space:]]" "${FULL_MASK_FILE_ABS_PATH}" | awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' | sort -k1,1V -k2,2n > "${MASK_FILE_PLAIN_TEMP}"
pipe_status_step5_prep=("${PIPESTATUS[@]}")
if [ ${pipe_status_step5_prep[0]} -gt 1 ] || [ ${pipe_status_step5_prep[1]} -ne 0 ] || [ ${pipe_status_step5_prep[2]} -ne 0 ]; then 
    log "ERROR: Mask preparation pipeline failed. Grep: ${pipe_status_step5_prep[0]}, Awk: ${pipe_status_step5_prep[1]}, Sort: ${pipe_status_step5_prep[2]}."
    exit 1
fi
if [ ! -s "${MASK_FILE_PLAIN_TEMP}" ]; then
    log "WARNING: Temporary plain mask file for ${TARGET_CHROMOSOME_ID} (${MASK_FILE_PLAIN_TEMP}) is empty."
fi
log "Compressing mask file ${MASK_FILE_PLAIN_TEMP} to ${MASK_FILE_FINAL}..."
bgzip -f -c "${MASK_FILE_PLAIN_TEMP}" > "${MASK_FILE_FINAL}"
log "Indexing gzipped mask file: ${MASK_FILE_FINAL}"
tabix -f -p bed "${MASK_FILE_FINAL}"
rm -f "${MASK_FILE_PLAIN_TEMP}"
log "--- Step 5 (Prepare Genome Mask) finished. Mask file is: ${MASK_FILE_FINAL} ---"

# --- Step 6: Run smc++ vcf2smc ---
# ... (Step 6 logic - exactly as in your last provided log - it was working, including DEBUGs) ...
log "--- Starting Step 6: Run smc++ vcf2smc ---"
log "DEBUG STEP6: Value of SAMPLES_FOR_SMC array count: ${#SAMPLES_FOR_SMC[@]}"
if [ ${#SAMPLES_FOR_SMC[@]} -eq 0 ]; then
    log "ERROR DEBUG STEP6: SAMPLES_FOR_SMC array is empty! Cannot proceed."
    exit 1
fi
log "DEBUG STEP6: First sample in SAMPLES_FOR_SMC: ${SAMPLES_FOR_SMC[0]}"
log "DEBUG STEP6: Value of DISTINGUISHED_INDIVIDUALS_TARGET_POP array: (${DISTINGUISHED_INDIVIDUALS_TARGET_POP[*]})"
log "DEBUG STEP6: DI1: '${DISTINGUISHED_INDIVIDUALS_TARGET_POP[0]}', DI2: '${DISTINGUISHED_INDIVIDUALS_TARGET_POP[1]}'"
SMC_SAMPLE_LIST_STR="${TARGET_POP_LABEL}:$(IFS=,; echo "${SAMPLES_FOR_SMC[*]}")" 
DI1="${DISTINGUISHED_INDIVIDUALS_TARGET_POP[0]}"
DI2="${DISTINGUISHED_INDIVIDUALS_TARGET_POP[1]}"
if [ -z "${DI1}" ] || [ -z "${DI2}" ]; then
    log "ERROR DEBUG STEP6: Distinguished individuals DI1 ('${DI1}') or DI2 ('${DI2}') are empty!"
    exit 1
fi
log "Using VCF for vcf2smc: ${INPUT_VCF_CURRENT_STEP}" 
log "This VCF contains ALL samples from the original input VCF (but variants are heavily filtered)."
log "smc++ will only consider samples specified in SMC_SAMPLE_LIST_STR."
log "Output SMC file: ${SMC_FORMAT_FILE}"
log "Population and samples for SMC (these are the common samples): ${SMC_SAMPLE_LIST_STR}"
log "Distinguished individuals for -d option (from common set): ${DI1} ${DI2}"
log "Mask file for -m option: ${MASK_FILE_FINAL}"
smc++ vcf2smc \
    --cores "${SMC_CORES}" \
    -m "${MASK_FILE_FINAL}" \
    -d "${DI1}" "${DI2}" \
    "${INPUT_VCF_CURRENT_STEP}" \
    "${SMC_FORMAT_FILE}" \
    "${TARGET_CHROMOSOME_ID}" \
    "${SMC_SAMPLE_LIST_STR}"
log "smc++ vcf2smc completed. Output: ${SMC_FORMAT_FILE}"
log "--- Step 6 (smc++ vcf2smc) finished. ---"

# --- Step 7: Run smc++ estimate ---
# ... (Step 7 logic - exactly as in your last provided log - it was working) ...
log "--- Starting Step 7: Run smc++ estimate for Initial Demographic Model ---"
if [ ! -f "${SMC_FORMAT_FILE}" ] || [ ! -s "${SMC_FORMAT_FILE}" ]; then
    log "ERROR: SMC format file from vcf2smc not found or is empty: ${SMC_FORMAT_FILE}"
    exit 1
fi
log "Running smc++ estimate with uniform recombination rate: ${UNIFORM_RECOMBINATION_RATE}"
smc++ estimate \
    --cores "${SMC_CORES}" \
    -o "${SMC_ESTIMATE_INITIAL_DIR}" \
    "${UNIFORM_RECOMBINATION_RATE}" \
    "${SMC_FORMAT_FILE}"
log "smc++ estimate (initial model) completed. Model files should be in ${SMC_ESTIMATE_INITIAL_DIR}"
log "Key output file: ${SMC_ESTIMATE_INITIAL_DIR}/model.final.json"
if [ -f "${SMC_ESTIMATE_INITIAL_DIR}/model.final.json" ]; then
    log "Initial model JSON found. Details:"
    ls -lh "${SMC_ESTIMATE_INITIAL_DIR}/model.final.json"
else
    log "ERROR: ${SMC_ESTIMATE_INITIAL_DIR}/model.final.json not created by smc++ estimate."
    log "Listing contents of ${SMC_ESTIMATE_INITIAL_DIR}:"
    ls -lR "${SMC_ESTIMATE_INITIAL_DIR}"
    exit 1
fi
log "--- Step 7 (smc++ estimate) finished. ---"

# --- Cleanup ---
log "--- Cleaning up temporary files ---"
rm -f "${TEMP_RENAMED_HEADER_FILE}" \
      "${TEMP_RENAMED_HEADER_FILE}.final" \
      "${TEMP_RENAMED_HEADER_FILE}.bak" \ 
      "${TEMP_CHROM_MAP_FILE}" \
      "${VCF_BODY_RENAMED}" "${VCF_BODY_RENAMED}.tbi" "${VCF_BODY_RENAMED}.csi" \
      "${VCF_CHR_NORMALIZED_RAW}" "${VCF_CHR_NORMALIZED_RAW}.tbi" "${VCF_CHR_NORMALIZED_RAW}.csi" \
      "${VCF_NORM_SPLIT_TEMP}" "${VCF_NORM_SPLIT_TEMP}.tbi" \ # Added VCF_NORM_SPLIT_TEMP
      "${VCF_SAMPLES_TEMP}" "${MASTER_SAMPLES_SORTED_TEMP}" \
      # Keep COMMON_SAMPLES_LIST_FILE for records, but could be added to rm
      # Keep Python log files
      2>/dev/null || log "Note: Some temporary files might have already been cleaned or not created."
log "--- Temporary file cleanup finished. ---"
log "--- SUCCESS: All processing steps completed for ${TARGET_POP_LABEL} on ${TARGET_CHROMOSOME_ID}. ---"
log "SMC++ formatted data is: ${SMC_FORMAT_FILE}"
log "Initial demographic model (JSON) is in: ${SMC_ESTIMATE_INITIAL_DIR}/model.final.json"
log "Intermediate VCF processing files are in ${VCF2SMC_INTERMEDIATE_DIR}"
log "Python script diagnostic logs are in relevant log files in ${VCF2SMC_INTERMEDIATE_DIR}"
