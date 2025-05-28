#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status
set -o pipefail # Causes a pipeline to return the exit status of the last command in the pipe that failed

# --- Configuration ---
ORIGINAL_INPUT_VCF_RELPATH="PJL/PJL.chr22.vcf.gz" 
AA_TSV_DIR_ABS_PATH="/vast/palmer/pi/reilly/jfa38/datasets_for_annotation/ErinG_ARG_hg19_AncestralAlleles/tidied_AA_tsvs"
FULL_MASK_FILE_RELPATH="20141020.strict_mask.whole_genome.bed"
PJL_SAMPLES_FILE_RELPATH="PJL.samples.list"
TARGET_CHROMOSOME_ID="chr22" 
UNIFORM_RECOMBINATION_RATE="1.25e-8"
SMC_CORES=4
REFERENCE_FASTA="/path/to/your/hg19.fa" # <--- !!! SET THIS PATH or leave as is for norm without -f !!!
PYTHON_ANNOTATE_SCRIPT_NAME="annotate_aa.py" # Ensure this script is present and executable
# --- End Configuration ---

TARGET_CHROMOSOME_NUMERIC_PART=$(echo "${TARGET_CHROMOSOME_ID}" | sed 's/[^0-9]*//g') 

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
PYTHON_ANNOTATE_SCRIPT_PATH="${SCRIPT_DIR}/${PYTHON_ANNOTATE_SCRIPT_NAME}"

MAIN_OUTPUT_DIR="${SCRIPT_DIR}/smc_processing_chr${TARGET_CHROMOSOME_NUMERIC_PART}"
VCF2SMC_INTERMEDIATE_DIR="${MAIN_OUTPUT_DIR}/vcf2smc_intermediate_files"
SMC_ESTIMATE_INITIAL_DIR="${MAIN_OUTPUT_DIR}/smc_estimate_initial_model"

mkdir -p "${VCF2SMC_INTERMEDIATE_DIR}"
mkdir -p "${SMC_ESTIMATE_INITIAL_DIR}"

ORIGINAL_INPUT_VCF_ABS_PATH="${SCRIPT_DIR}/${ORIGINAL_INPUT_VCF_RELPATH}" 
AA_TSV_ORIGINAL_BASENAME="chr${TARGET_CHROMOSOME_NUMERIC_PART}_ErinG_AA_hg19.tsv.gz"
AA_TSV_ORIGINAL="${AA_TSV_DIR_ABS_PATH}/${AA_TSV_ORIGINAL_BASENAME}"
FULL_MASK_FILE_ABS_PATH="${SCRIPT_DIR}/${FULL_MASK_FILE_RELPATH}" # Use absolute path for consistency
PJL_SAMPLES_FILE_ABS_PATH="${SCRIPT_DIR}/${PJL_SAMPLES_FILE_RELPATH}" # Use absolute path

VCF_BODY_RENAMED="${VCF2SMC_INTERMEDIATE_DIR}/PJL.chr${TARGET_CHROMOSOME_NUMERIC_PART}.phased.bodyRenamed.vcf.gz"
TEMP_RENAMED_HEADER_FILE="${VCF2SMC_INTERMEDIATE_DIR}/temp_renamed_header.txt"
TEMP_CHROM_MAP_FILE="${VCF2SMC_INTERMEDIATE_DIR}/temp_chromosome_map.txt"
VCF_CHR_NORMALIZED_RAW="${VCF2SMC_INTERMEDIATE_DIR}/PJL.chr${TARGET_CHROMOSOME_NUMERIC_PART}.phased.chrNormalized.raw.vcf.gz"
VCF_CHR_NORMALIZED_CLEANED="${VCF2SMC_INTERMEDIATE_DIR}/PJL.chr${TARGET_CHROMOSOME_NUMERIC_PART}.phased.chrNormalized.cleaned.vcf.gz"
AA_TSV_PREPARED="${VCF2SMC_INTERMEDIATE_DIR}/chr${TARGET_CHROMOSOME_NUMERIC_PART}_ErinG_AA_hg19_prepared.tsv.gz"
VCF_ANNOTATED_PY="${VCF2SMC_INTERMEDIATE_DIR}/PJL.chr${TARGET_CHROMOSOME_NUMERIC_PART}.phased.annotated.py.vcf.gz"
MASK_FILE_FINAL="${VCF2SMC_INTERMEDIATE_DIR}/mask.${TARGET_CHROMOSOME_ID}.bed.gz"
SMC_FORMAT_FILE="${VCF2SMC_INTERMEDIATE_DIR}/PJL.chr${TARGET_CHROMOSOME_NUMERIC_PART}.smc.gz"
POP_LABEL="PJL"

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] (${TARGET_CHROMOSOME_ID}) $1"
}

# --- Tool Check and Script Existence ---
log "--- Sanity Checks ---"
for tool in bcftools tabix bgzip smc++ samtools python3; do
    if ! command -v $tool &> /dev/null; then
        log "ERROR: Required tool '$tool' not found. Please ensure it is installed and in your PATH."
        exit 1
    fi
done
log "All required command-line tools found."

if [ ! -f "${PYTHON_ANNOTATE_SCRIPT_PATH}" ]; then
    log "ERROR: Python annotation script ${PYTHON_ANNOTATE_SCRIPT_PATH} not found."
    log "Please ensure it exists in the same directory as this bash script, or update PYTHON_ANNOTATE_SCRIPT_NAME."
    exit 1
fi
log "Python annotation script found at ${PYTHON_ANNOTATE_SCRIPT_PATH}."
log "--- End Sanity Checks ---"

# --- Step 0: Read Samples and Define Distinguished Individuals ---
log "--- Starting Step 0: Read Samples and Define Distinguished Individuals ---"
if [ ! -f "${PJL_SAMPLES_FILE_ABS_PATH}" ]; then
    log "ERROR: PJL samples list not found: ${PJL_SAMPLES_FILE_ABS_PATH}"
    exit 1
fi
mapfile -t SAMPLES_PJL < <(cat "${PJL_SAMPLES_FILE_ABS_PATH}")
if [ ${#SAMPLES_PJL[@]} -lt 2 ]; then
    log "ERROR: Not enough samples found in ${PJL_SAMPLES_FILE_ABS_PATH} for Distinguished Individuals. Need at least 2."
    exit 1
fi
DISTINGUISHED_INDIVIDUALS_PJL=("${SAMPLES_PJL[0]}" "${SAMPLES_PJL[1]}") # Takes the first two samples
log "Total samples for ${POP_LABEL}: ${#SAMPLES_PJL[@]}"
log "Distinguished Individuals for ${POP_LABEL}: ${DISTINGUISHED_INDIVIDUALS_PJL[0]}, ${DISTINGUISHED_INDIVIDUALS_PJL[1]}"
log "--- Step 0 finished. ---"

# --- Step 1: Normalize Input VCF Chromosome Names ---
log "--- Starting Step 1: Normalize Input VCF Chromosome Names & Clean ---"
log "Path to original VCF: '${ORIGINAL_INPUT_VCF_ABS_PATH}'"
if [ ! -f "${ORIGINAL_INPUT_VCF_ABS_PATH}" ]; then
    log "ERROR: Original VCF not found: ${ORIGINAL_INPUT_VCF_ABS_PATH}"
    exit 1
fi
if [ ! -f "${ORIGINAL_INPUT_VCF_ABS_PATH}.tbi" ] && [ ! -f "${ORIGINAL_INPUT_VCF_ABS_PATH}.csi" ]; then
    log "ERROR: Index for original VCF (${ORIGINAL_INPUT_VCF_ABS_PATH}) not found (.tbi or .csi)."
    exit 1
fi
log "Original VCF and its index seem OK."

original_vcf_chrom_name_sample=""
_command_output=""
_exit_status=0
set +e # Temporarily disable exit on error for this command
_command_output=$(zcat "${ORIGINAL_INPUT_VCF_ABS_PATH}" 2>/dev/null | grep -v "^#" | head -n 1 | awk '{print $1}')
_exit_status=$?
set -e # Re-enable exit on error

if [ ${_exit_status} -eq 0 ] && [ -n "$_command_output" ]; then
    original_vcf_chrom_name_sample="$_command_output"
    log "Cmd to get chrom name sample succeeded. Name: '${original_vcf_chrom_name_sample}'"
elif [ ${_exit_status} -ne 0 ] && [ -n "$_command_output" ]; then # e.g. SIGPIPE (141) but still got output
    original_vcf_chrom_name_sample="$_command_output"
    log "Cmd to get chrom name sample had non-zero exit (${_exit_status}) but produced output. Name: '${original_vcf_chrom_name_sample}'. Proceeding cautiously."
else
    log "ERROR: Cmd to get chrom name sample failed. Exit: ${_exit_status}. Output: '$_command_output'."
    log "First 20 lines of ${ORIGINAL_INPUT_VCF_ABS_PATH}:"
    zcat "${ORIGINAL_INPUT_VCF_ABS_PATH}" | head -n 20
    exit 1
fi
if [ -z "${original_vcf_chrom_name_sample}" ]; then
    log "ERROR: Could not determine chromosome name sample from ${ORIGINAL_INPUT_VCF_ABS_PATH}."
    exit 1
fi
log "Successfully determined original VCF chromosome name sample: '${original_vcf_chrom_name_sample}'"

if [ "${original_vcf_chrom_name_sample}" != "${TARGET_CHROMOSOME_ID}" ]; then
    log "Normalization required: VCF has '${original_vcf_chrom_name_sample}', target is '${TARGET_CHROMOSOME_ID}'."
    echo -e "${original_vcf_chrom_name_sample}\t${TARGET_CHROMOSOME_ID}" > "${TEMP_CHROM_MAP_FILE}"
    log "Step 1a: Renaming chromosome names in VCF records..."
    bcftools annotate --rename-chrs "${TEMP_CHROM_MAP_FILE}" "${ORIGINAL_INPUT_VCF_ABS_PATH}" -O z -o "${VCF_BODY_RENAMED}"
    log "Step 1b: Fixing header of body-renamed VCF..."
    bcftools view -h "${VCF_BODY_RENAMED}" > "${TEMP_RENAMED_HEADER_FILE}"
    # More robust sed: handle cases where length might be present or not
    sed -E -e "s/##contig=<ID=${original_vcf_chrom_name_sample}([,>])/##contig=<ID=${TARGET_CHROMOSOME_ID}\\1/g" "${TEMP_RENAMED_HEADER_FILE}" > "${TEMP_RENAMED_HEADER_FILE}.final"
    mv "${TEMP_RENAMED_HEADER_FILE}.final" "${TEMP_RENAMED_HEADER_FILE}"
    
    # Check if target contig exists, if not, add a basic one (common for chr22 length)
    if ! grep -q "##contig=<ID=${TARGET_CHROMOSOME_ID}" "${TEMP_RENAMED_HEADER_FILE}"; then
        log "WARNING: Contig '${TARGET_CHROMOSOME_ID}' not found in header after sed rename. Adding a basic definition."
        known_length_target=""
        if [ "${TARGET_CHROMOSOME_ID}" == "chr22" ]; then known_length_target="51304566"; fi # Approx hg19 chr22 length
        
        if [ -n "$known_length_target" ]; then
            echo "##contig=<ID=${TARGET_CHROMOSOME_ID},length=${known_length_target}>" >> "${TEMP_RENAMED_HEADER_FILE}"
        else
            echo "##contig=<ID=${TARGET_CHROMOSOME_ID}>" >> "${TEMP_RENAMED_HEADER_FILE}" # Fallback if length unknown
        fi
        log "Added contig definition for ${TARGET_CHROMOSOME_ID} to header."
    fi
    log "Reheadering body-renamed VCF with corrected header to ${VCF_CHR_NORMALIZED_RAW}..."
    bcftools reheader -h "${TEMP_RENAMED_HEADER_FILE}" "${VCF_BODY_RENAMED}" -o "${VCF_CHR_NORMALIZED_RAW}"
else
    log "Chromosome names in original VCF ('${original_vcf_chrom_name_sample}') already match target ('${TARGET_CHROMOSOME_ID}'). Copying to ${VCF_CHR_NORMALIZED_RAW}."
    cp "${ORIGINAL_INPUT_VCF_ABS_PATH}" "${VCF_CHR_NORMALIZED_RAW}"
fi

log "Step 1c: Cleaning the VCF (${VCF_CHR_NORMALIZED_RAW}) using bcftools norm..."
if [ ! -f "${REFERENCE_FASTA}" ] || [ "${REFERENCE_FASTA}" == "/path/to/your/hg19.fa" ]; then 
    log "WARNING: Reference FASTA ${REFERENCE_FASTA} not found or default path used. Running bcftools norm without -f (reference genome). This is less effective for normalization (e.g. for multi-allelic splitting, left-alignment)."
    bcftools norm -m -any "${VCF_CHR_NORMALIZED_RAW}" -O z -o "${VCF_CHR_NORMALIZED_CLEANED}" --threads "${SMC_CORES}"
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
    bcftools norm -m -any -f "${REFERENCE_FASTA}" "${VCF_CHR_NORMALIZED_RAW}" -O z -o "${VCF_CHR_NORMALIZED_CLEANED}" --threads "${SMC_CORES}"
fi
log "bcftools norm completed."
log "Indexing cleaned VCF: ${VCF_CHR_NORMALIZED_CLEANED}"
tabix -f -p vcf "${VCF_CHR_NORMALIZED_CLEANED}"
INPUT_VCF_CURRENT_STEP="${VCF_CHR_NORMALIZED_CLEANED}"
log "--- Step 1 (Chrom Name Normalization & Cleaning) finished. Final VCF for this step: ${INPUT_VCF_CURRENT_STEP} ---"

# --- Step 2: Prepare Ancestral Allele TSV ---
log "--- Starting Step 2: Prepare Ancestral Allele TSV ---"
if [ ! -f "${AA_TSV_ORIGINAL}" ]; then
    log "ERROR: Original AA TSV file not found: ${AA_TSV_ORIGINAL}"
    log "Please check AA_TSV_DIR_ABS_PATH and AA_TSV_ORIGINAL_BASENAME settings."
    exit 1
fi

log "DEBUG_STEP2: First 5 lines of ORIGINAL AA TSV (${AA_TSV_ORIGINAL}):"
zcat "${AA_TSV_ORIGINAL}" | head -n 5 || log "DEBUG_STEP2_WARNING: Failed to read head of ${AA_TSV_ORIGINAL}"
log "DEBUG_STEP2: Line count of ORIGINAL AA TSV (${AA_TSV_ORIGINAL}):"
original_aa_line_count=$(zcat "${AA_TSV_ORIGINAL}" | wc -l) || original_aa_line_count="Error"
log "DEBUG_STEP2: Original AA TSV line count: ${original_aa_line_count}"

log "Preparing AA TSV: ${AA_TSV_PREPARED} from ${AA_TSV_ORIGINAL}"
# Robust awk: Handles if original TSV has a header or not. Assumes $1=chrom (may be ignored), $2=pos, $3=AA in original data lines
# The Python script will skip any line starting with #, so the exact header content here is mostly for human readability of the prepared TSV.
zcat "${AA_TSV_ORIGINAL}" | \
    awk -v target_chrom_awk="${TARGET_CHROMOSOME_ID}" '
        BEGIN {OFS="\t"}
        # If the first line of the original file is NOT a comment, print our desired header.
        (NR==1 && $0 !~ /^#/) { print "#CHROM_AA", "POS_AA", "AA_VALUE" }
        # If the line is NOT a comment, print the data using target_chrom_awk and original $2, $3
        # This assumes original data has at least 3 fields, and $2 is pos, $3 is AA.
        # It also assumes that if $1 was a chromosome, it might not match target_chrom_awk, so we force target_chrom_awk.
        ($0 !~ /^#/) { if (NF >=3) print target_chrom_awk, $2, $3 }
        # If the line IS a comment (e.g. existing header in original), print it as is.
        ($0 ~ /^#/) { print $0 }
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
tabix -f -S 1 -s 1 -b 2 -e 2 "${AA_TSV_PREPARED}" # Preset for TSV: skip 1 line (header), col1=seq, col2=start, col2=end
log "--- Step 2 (Prepare AA TSV) finished. Prepared AA TSV is: ${AA_TSV_PREPARED} ---"

# --- Step 3: Annotate VCF with Ancestral Alleles using Python script ---
log "--- Starting Step 3: Annotate VCF with Ancestral Alleles using Python ---"
log "Input VCF for Python annotation: ${INPUT_VCF_CURRENT_STEP}"
log "AA TSV for Python annotation: ${AA_TSV_PREPARED}"
log "Output VCF from Python annotation: ${VCF_ANNOTATED_PY}"

PYTHON_STDOUT_LOG="${VCF2SMC_INTERMEDIATE_DIR}/python_annotate_stdout.log"
PYTHON_STDERR_LOG="${VCF2SMC_INTERMEDIATE_DIR}/python_annotate_stderr.log"

log "Python script stdout will be logged to: ${PYTHON_STDOUT_LOG}"
log "Python script stderr (diagnostics) will be logged to: ${PYTHON_STDERR_LOG}"

# Execute Python script and capture its stdout and stderr
python3 "${PYTHON_ANNOTATE_SCRIPT_PATH}" \
    "${INPUT_VCF_CURRENT_STEP}" \
    "${AA_TSV_PREPARED}" \
    "${VCF_ANNOTATED_PY}" \
    1> "${PYTHON_STDOUT_LOG}" \
    2> "${PYTHON_STDERR_LOG}"
PYTHON_EXIT_STATUS=$?

# Display the content of stderr log to the main script log for easier viewing
log "--- Content of Python stderr log (${PYTHON_STDERR_LOG}) ---"
if [ -f "${PYTHON_STDERR_LOG}" ]; then
    cat "${PYTHON_STDERR_LOG}" # This will make it part of the main bash log
else
    log "WARNING: Python stderr log file not found: ${PYTHON_STDERR_LOG}"
fi
log "--- End of Python stderr log content ---"

# Display the content of stdout log if it's not empty
if [ -s "${PYTHON_STDOUT_LOG}" ]; then # -s checks if file exists and is not empty
    log "--- Content of Python stdout log (${PYTHON_STDOUT_LOG}) ---"
    cat "${PYTHON_STDOUT_LOG}"
    log "--- End of Python stdout log content ---"
fi

if [ ${PYTHON_EXIT_STATUS} -ne 0 ]; then
    log "ERROR: Python annotation script failed with exit status ${PYTHON_EXIT_STATUS}."
    log "Check ${PYTHON_STDERR_LOG} and ${PYTHON_STDOUT_LOG} (if any content) for details."
    exit 1
fi
log "Python annotation script appears to have finished successfully (exit status 0)."

log "Indexing Python-annotated VCF: ${VCF_ANNOTATED_PY}"
tabix -f -p vcf "${VCF_ANNOTATED_PY}"

INPUT_VCF_CURRENT_STEP="${VCF_ANNOTATED_PY}" 
log "--- Step 3 (Python AA Annotation) finished. Annotated VCF is: ${INPUT_VCF_CURRENT_STEP} ---"

# --- Step 5: Prepare Genome Mask for the target chromosome ---
# (Step 4 was effectively part of Step 3 with the Python script logic)
log "--- Starting Step 5: Prepare Genome Mask ---"
MASK_FILE_PLAIN_TEMP="${VCF2SMC_INTERMEDIATE_DIR}/mask.${TARGET_CHROMOSOME_ID}.temp.plain.bed"

if [ ! -f "${FULL_MASK_FILE_ABS_PATH}" ]; then
    log "ERROR: Full genome mask file not found: ${FULL_MASK_FILE_ABS_PATH}"
    exit 1
fi

log "Extracting and sorting mask regions for ${TARGET_CHROMOSOME_ID} from ${FULL_MASK_FILE_ABS_PATH}..."
# grep for lines starting with the target chromosome, followed by a tab or space
# awk to ensure correct 3-column BED format (though original is likely already fine)
# sort as required by tabix for BED files
grep "^${TARGET_CHROMOSOME_ID}[[:space:]]" "${FULL_MASK_FILE_ABS_PATH}" | awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' | sort -k1,1V -k2,2n > "${MASK_FILE_PLAIN_TEMP}"
pipe_status_step5_prep=("${PIPESTATUS[@]}")

if [ ${pipe_status_step5_prep[0]} -gt 1 ]; then 
    log "ERROR: grep for mask regions failed with status ${pipe_status_step5_prep[0]}."
    exit 1
elif [ ${pipe_status_step5_prep[1]} -ne 0 ]; then
    log "ERROR: awk for mask regions failed with status ${pipe_status_step5_prep[1]}."
    exit 1
elif [ ${pipe_status_step5_prep[2]} -ne 0 ]; then
    log "ERROR: sort for mask regions failed with status ${pipe_status_step5_prep[2]}."
    exit 1
fi

if [ ! -s "${MASK_FILE_PLAIN_TEMP}" ]; then
    log "WARNING: Temporary plain mask file for ${TARGET_CHROMOSOME_ID} (${MASK_FILE_PLAIN_TEMP}) is empty after extraction and sorting."
    log "This means no mask regions were found for this chromosome in ${FULL_MASK_FILE_ABS_PATH}."
    log "smc++ will proceed without masking for this chromosome if this file remains effectively empty."
fi

log "Compressing mask file ${MASK_FILE_PLAIN_TEMP} to ${MASK_FILE_FINAL}..."
bgzip -f -c "${MASK_FILE_PLAIN_TEMP}" > "${MASK_FILE_FINAL}"
log "Indexing gzipped mask file: ${MASK_FILE_FINAL}"
tabix -f -p bed "${MASK_FILE_FINAL}"
rm -f "${MASK_FILE_PLAIN_TEMP}"
log "--- Step 5 (Prepare Genome Mask) finished. Mask file is: ${MASK_FILE_FINAL} ---"

# --- Step 6: Run smc++ vcf2smc ---
log "--- Starting Step 6: Run smc++ vcf2smc ---"
SMC_SAMPLE_LIST_STR="${POP_LABEL}:$(IFS=,; echo "${SAMPLES_PJL[*]}")" 
DI1="${DISTINGUISHED_INDIVIDUALS_PJL[0]}"
DI2="${DISTINGUISHED_INDIVIDUALS_PJL[1]}"

log "Using VCF for vcf2smc: ${INPUT_VCF_CURRENT_STEP}"
log "Output SMC file: ${SMC_FORMAT_FILE}"
log "Population and samples for SMC: ${SMC_SAMPLE_LIST_STR}"
log "Distinguished individuals for -d option: ${DI1} ${DI2}"
log "Mask file for -m option: ${MASK_FILE_FINAL}"

smc++ vcf2smc \
    --cores "${SMC_CORES}" \
    -m "${MASK_FILE_FINAL}" \
    -d "${DI1}" "${DI2}" \
    "${INPUT_VCF_CURRENT_STEP}" \
    "${SMC_FORMAT_FILE}" \
    "${TARGET_CHROMOSOME_ID}" \
    "${SMC_SAMPLE_LIST_STR}"
# smc++ should auto-detect INFO/AA if the tag is present and correctly formatted (base or '.')

log "smc++ vcf2smc completed. Output: ${SMC_FORMAT_FILE}"
log "--- Step 6 (smc++ vcf2smc) finished. ---"

# --- Step 7: Run smc++ estimate for Initial Demographic Model ---
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

# --- Cleanup Temporary Files ---
log "--- Cleaning up temporary files ---"
rm -f "${TEMP_RENAMED_HEADER_FILE}" \
      "${TEMP_RENAMED_HEADER_FILE}.final" \
      "${TEMP_CHROM_MAP_FILE}" \
      "${VCF_BODY_RENAMED}" "${VCF_BODY_RENAMED}.tbi" "${VCF_BODY_RENAMED}.csi" \
      "${VCF_CHR_NORMALIZED_RAW}" "${VCF_CHR_NORMALIZED_RAW}.tbi" "${VCF_CHR_NORMALIZED_RAW}.csi" \
      2>/dev/null || log "Note: Some temporary files might have already been cleaned or not created."
# MASK_FILE_PLAIN_TEMP is removed within Step 5.
# Python log files PYTHON_STDOUT_LOG and PYTHON_STDERR_LOG are kept for debugging.
log "--- Temporary file cleanup finished. ---"

log "--- SUCCESS: All processing steps completed. ---"
log "SMC++ formatted data for ${TARGET_CHROMOSOME_ID} is: ${SMC_FORMAT_FILE}"
log "Initial demographic model (JSON) is in: ${SMC_ESTIMATE_INITIAL_DIR}/model.final.json"
log "Intermediate files for VCF to SMC format conversion are in ${VCF2SMC_INTERMEDIATE_DIR}"
log "Python script diagnostic logs are in ${PYTHON_STDOUT_LOG} and ${PYTHON_STDERR_LOG}"
