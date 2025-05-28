import pysam
import sys
import gzip

def load_ancestral_alleles(aa_tsv_gz_path):
    aa_map = {}
    log_prefix = "PY_TSV_LOADER" 
    sys.stderr.write(f"{log_prefix}: Attempting to load ancestral alleles from: {aa_tsv_gz_path}\n")
    count = 0
    valid_bases = {'A', 'C', 'G', 'T'}
    header_lines_skipped = 0
    malformed_data_lines = 0

    try:
        with gzip.open(aa_tsv_gz_path, 'rt') as f_aa:
            for line_num, line in enumerate(f_aa):
                line_stripped = line.strip()
                if not line_stripped: continue
                if line_stripped.startswith('#'):
                    if header_lines_skipped == 0:
                        sys.stderr.write(f"  {log_prefix}: Skipping header line in TSV: {line_stripped}\n")
                    header_lines_skipped += 1
                    continue
                
                parts = line_stripped.split('\t')
                if len(parts) >= 3:
                    chrom, pos_str, ancestral_allele = parts[0], parts[1], parts[2].upper()
                    try:
                        pos = int(pos_str)
                    except ValueError:
                        sys.stderr.write(f"  {log_prefix}_WARNING: Malformed position '{pos_str}' at chr='{chrom}' in TSV (line {line_num + 1}): \"{line_stripped}\"\n")
                        malformed_data_lines += 1
                        continue

                    if ancestral_allele in valid_bases:
                        aa_map[(chrom, pos)] = ancestral_allele
                        count += 1
                    elif ancestral_allele == 'N' or ancestral_allele == '.':
                        aa_map[(chrom, pos)] = '.'
                        count += 1
                    else:
                        sys.stderr.write(f"  {log_prefix}_WARNING: Invalid AA '{ancestral_allele}' at {chrom}:{pos} in TSV (line {line_num + 1}). Expected A,C,G,T,N, or '. Value was: '{parts[2]}'\n")
                        malformed_data_lines += 1
                else:
                    sys.stderr.write(f"  {log_prefix}_WARNING: Malformed data line {line_num + 1} in AA TSV (expected >=3 columns, got {len(parts)}): \"{line_stripped}\"\n")
                    malformed_data_lines += 1
    except FileNotFoundError:
        sys.stderr.write(f"  {log_prefix}_ERROR: Ancestral allele TSV file not found: {aa_tsv_gz_path}\n")
        return None 
    except Exception as e:
        sys.stderr.write(f"  {log_prefix}_ERROR: Failed to read/process ancestral allele TSV {aa_tsv_gz_path}: {e}\n")
        return None

    sys.stderr.write(f"  {log_prefix}: Total header lines skipped: {header_lines_skipped}\n")
    if malformed_data_lines > 0:
        sys.stderr.write(f"  {log_prefix}_WARNING: Total malformed/skipped data lines: {malformed_data_lines}\n")
    sys.stderr.write(f"  {log_prefix}: Loaded {count} valid ancestral allele entries from TSV.\n")
    if count == 0: # Further refined this warning
        if header_lines_skipped > 0 and malformed_data_lines == 0 and count == 0 :
             sys.stderr.write(f"  {log_prefix}_WARNING: No valid data entries loaded from {aa_tsv_gz_path}, though headers were present. Check data format.\n")
        else:
            sys.stderr.write(f"  {log_prefix}_WARNING: No ancestral allele entries were loaded at all from {aa_tsv_gz_path}. Check file format and content.\n")
    return aa_map

def main(input_vcf_path, aa_tsv_path, output_vcf_path):
    log_prefix = "PY_ANNOTATE" 
    sys.stderr.write(f"{log_prefix}: Input VCF: {input_vcf_path}\n")
    sys.stderr.write(f"{log_prefix}: AA TSV: {aa_tsv_path}\n")
    sys.stderr.write(f"{log_prefix}: Output VCF: {output_vcf_path}\n")

    aa_data = load_ancestral_alleles(aa_tsv_path)
    if aa_data is None: 
        sys.stderr.write(f"{log_prefix}_ERROR: Aborting annotation due to failure in loading ancestral alleles.\n")
        sys.exit(1)
    
    try:
        vcf_in = pysam.VariantFile(input_vcf_path, "rb") 
    except Exception as e:
        sys.stderr.write(f"{log_prefix}_ERROR: Could not open input VCF file {input_vcf_path}: {e}\n")
        sys.exit(1)

    header_out = vcf_in.header.copy()
    expected_aa_description = "Ancestral Allele (from TSV or '.' if unknown/not in TSV)"
    
    # Ensure AA INFO field is defined in the output header
    if 'AA' not in header_out.info:
        try:
            header_out.info.add('AA', 1, 'String', expected_aa_description)
            sys.stderr.write(f"{log_prefix}: Added new 'AA' INFO definition to output VCF header.\n")
        except Exception as e:
            sys.stderr.write(f"{log_prefix}_ERROR: Could not add 'AA' INFO field to header: {e}\n")
            vcf_in.close()
            sys.exit(1)
    else: # AA is already defined in the input header (e.g., for PJL VCFs)
        current_aa_def = header_out.info['AA']
        if not (current_aa_def.number == 1 and current_aa_def.type == 'String'):
            sys.stderr.write(f"{log_prefix}_ERROR: Pre-existing AA INFO definition is incompatible! Number={current_aa_def.number}, Type='{current_aa_def.type}'. Aborting.\n")
            vcf_in.close()
            sys.exit(1) 
        else:
            sys.stderr.write(f"{log_prefix}: Compatible 'AA' INFO definition already exists in input VCF header. Original description: '{current_aa_def.description}'.\n")
            # Optional: update description if desired, though usually not critical
            # header_out.info['AA'].description = expected_aa_description 

    try:
        # Open output VCF with the (potentially) modified header_out
        vcf_out = pysam.VariantFile(output_vcf_path, "w", header=header_out)
    except Exception as e:
        sys.stderr.write(f"{log_prefix}_ERROR: Could not open output VCF file {output_vcf_path} for writing: {e}\n")
        vcf_in.close()
        sys.exit(1)

    # Counters
    total_variants = 0
    aa_from_tsv_written_stat = 0
    aa_from_tsv_base_stat = 0
    aa_from_tsv_dot_stat = 0
    aa_set_to_dot_no_tsv_entry_stat = 0
    # Detailed stats about overwriting/adding (relevant if input VCF *had* AA tags)
    overwritten_existing_aa_with_tsv_val_stat = 0 
    added_new_aa_with_tsv_val_stat = 0
    overwritten_existing_aa_with_dot_no_tsv_stat = 0
    added_new_aa_as_dot_no_tsv_stat = 0


    for record_in in vcf_in:
        total_variants += 1
        if total_variants % 200000 == 0: 
            sys.stderr.write(f"  {log_prefix}_PROGRESS: Processed {total_variants} variants...\n")

        # Create a new record using the output file's header.
        # This ensures the record is aware of all INFO fields defined in header_out.
        record_out = vcf_out.new_record()
        record_out.chrom = record_in.chrom
        record_out.pos = record_in.pos
        if record_in.id is not None: record_out.id = record_in.id
        record_out.ref = record_in.ref
        record_out.alts = record_in.alts
        if record_in.qual is not None: record_out.qual = record_in.qual
        
        for f_id in record_in.filter:
            record_out.filter.add(f_id)
        
        # Copy FORMAT fields and sample genotypes
        # This assumes the sample list in header_out matches record_in.
        # Since header_out is a copy of vcf_in.header, this should be fine.
        for sample_name_in, sample_data_in in record_in.samples.items():
            for fmt_key, fmt_val in sample_data_in.items():
                 record_out.samples[sample_name_in][fmt_key] = fmt_val
        
        # Check if the original input record had an AA tag (for stats)
        # Note: record_in.info['AA'] would fail if INFO was '.' or AA not present
        original_aa_present_in_input_record = False
        if record_in.info is not None : # Check if INFO field itself exists (not just '.')
            try:
                if 'AA' in record_in.info:
                    original_aa_present_in_input_record = True
            except TypeError: # Handle if record_in.info is not dict-like (e.g. if it was just '.')
                 pass


        key = (record_out.chrom, record_out.pos) # Use record_out.chrom/pos
        ancestral_allele_from_tsv = aa_data.get(key)
        
        # Set the AA tag on the new record_out
        if ancestral_allele_from_tsv is not None:
            record_out.info['AA'] = ancestral_allele_from_tsv
            aa_from_tsv_written_stat += 1
            if ancestral_allele_from_tsv == '.': aa_from_tsv_dot_stat +=1
            else: aa_from_tsv_base_stat +=1
            
            if original_aa_present_in_input_record: overwritten_existing_aa_with_tsv_val_stat +=1
            else: added_new_aa_with_tsv_val_stat +=1
        else: 
            record_out.info['AA'] = '.'
            aa_set_to_dot_no_tsv_entry_stat += 1

            if original_aa_present_in_input_record: overwritten_existing_aa_with_dot_no_tsv_stat +=1
            else: added_new_aa_as_dot_no_tsv_stat +=1
            
        vcf_out.write(record_out)

    vcf_in.close()
    vcf_out.close()

    sys.stderr.write(f"\n--- Annotation Summary (Python Script v3.3) ---\n") 
    stat_prefix = f"{log_prefix}_STATS"
    sys.stderr.write(f"{stat_prefix}: Processed {total_variants} variants.\n")
    
    sys.stderr.write(f"\n{stat_prefix}: [TSV MATCH FOUND - AA tag set from TSV data]\n")
    sys.stderr.write(f"{stat_prefix}: Total sites with AA tag derived from TSV: {aa_from_tsv_written_stat}\n")
    sys.stderr.write(f"{stat_prefix}:   - TSV provided a base (A/C/G/T): {aa_from_tsv_base_stat}\n")
    sys.stderr.write(f"{stat_prefix}:   - TSV provided unknown (N or .) and AA set to '.': {aa_from_tsv_dot_stat}\n")
    sys.stderr.write(f"{stat_prefix}:     - Details: Overwrote existing AA tag in input with TSV value: {overwritten_existing_aa_with_tsv_val_stat}\n")
    sys.stderr.write(f"{stat_prefix}:     - Details: Added new AA tag (input had no AA) with TSV value: {added_new_aa_with_tsv_val_stat}\n")

    sys.stderr.write(f"\n{stat_prefix}: [NO TSV MATCH - AA tag set to '.']\n")
    sys.stderr.write(f"{stat_prefix}: Total sites where AA set to '.' because variant was NOT in TSV: {aa_set_to_dot_no_tsv_entry_stat}\n")
    sys.stderr.write(f"{stat_prefix}:   - Details: Overwrote existing AA tag in input with '.': {overwritten_existing_aa_with_dot_no_tsv_stat}\n")
    sys.stderr.write(f"{stat_prefix}:   - Details: Added new AA='.' tag (input had no AA): {added_new_aa_as_dot_no_tsv_stat}\n")
    
    sys.stderr.write(f"\n{log_prefix}: Output VCF written to: {output_vcf_path}\n")


if __name__ == "__main__":
    log_prefix_main = "PY_MAIN"
    if len(sys.argv) != 4:
        sys.stderr.write(f"{log_prefix_main}_ERROR: Usage: python3 annotate_aa.py <input.vcf.gz> <ancestral_alleles.tsv.gz> <output.vcf.gz>\n")
        sys.exit(1)
    
    input_vcf, aa_tsv, output_vcf = sys.argv[1], sys.argv[2], sys.argv[3]
    try:
        main(input_vcf, aa_tsv, output_vcf)
        sys.stderr.write(f"{log_prefix_main}: Python annotation script finished successfully.\n")
    except Exception as e:
        sys.stderr.write(f"{log_prefix_main}_ERROR: An unhandled error occurred in main(): {e}\n") # Python 3.6+ f-string
        # For older Python, use: sys.stderr.write("{}_ERROR: An unhandled error occurred in main(): {}\n".format(log_prefix_main, e))
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
