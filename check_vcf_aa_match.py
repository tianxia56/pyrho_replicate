import pysam
import sys
import gzip
import os # For temp file name

def load_ancestral_alleles_diagnostic(aa_tsv_gz_path):
    """Loads ancestral alleles from a gzipped TSV and prints diagnostics."""
    aa_map = {}
    print(f"\n--- DIAGNOSTIC: Loading Ancestral Alleles from: {aa_tsv_gz_path} ---")
    loaded_count = 0
    header_lines_found = 0
    malformed_lines = 0
    valid_bases = {'A', 'C', 'G', 'T'}

    try:
        with gzip.open(aa_tsv_gz_path, 'rt') as f_aa:
            for line_num, line in enumerate(f_aa):
                line_stripped = line.strip()
                if not line_stripped:
                    continue
                if line_stripped.startswith('#'):
                    if header_lines_found == 0: 
                        print(f"  AA_TSV_HEADER: Found header: {line_stripped}")
                    header_lines_found += 1
                    continue
                
                parts = line_stripped.split('\t')
                if len(parts) >= 3:
                    chrom, pos_str, ancestral_allele_val = parts[0], parts[1], parts[2].upper()
                    try:
                        pos = int(pos_str)
                    except ValueError:
                        print(f"  AA_TSV_WARNING: Malformed position '{pos_str}' at {chrom} (line {line_num + 1}): {line_stripped}")
                        malformed_lines += 1
                        continue

                    if ancestral_allele_val in valid_bases:
                        aa_map[(chrom, pos)] = ancestral_allele_val
                        loaded_count += 1
                    elif ancestral_allele_val == 'N' or ancestral_allele_val == '.':
                        aa_map[(chrom, pos)] = '.'
                        loaded_count += 1
                    else:
                        print(f"  AA_TSV_WARNING: Invalid AA '{ancestral_allele_val}' at {chrom}:{pos} (line {line_num + 1}).")
                        malformed_lines +=1
                else:
                    print(f"  AA_TSV_WARNING: Malformed line {line_num + 1} (expected >=3 columns): {line_stripped}")
                    malformed_lines += 1
        
        print(f"  AA_TSV_SUMMARY: Total header lines skipped: {header_lines_found}")
        print(f"  AA_TSV_SUMMARY: Total malformed/skipped data lines: {malformed_lines}")
        print(f"  AA_TSV_SUMMARY: Loaded {loaded_count} valid ancestral allele entries.")
        if loaded_count == 0 and header_lines_found > 0 and malformed_lines == 0:
             print(f"  AA_TSV_WARNING: No data lines found after header(s) in {aa_tsv_gz_path}")

    except FileNotFoundError:
        print(f"  AA_TSV_ERROR: File not found: {aa_tsv_gz_path}")
        return None
    except Exception as e:
        print(f"  AA_TSV_ERROR: Could not read {aa_tsv_gz_path}: {e}")
        return None
    
    print(f"--- DIAGNOSTIC: Finished loading Ancestral Alleles ---\n")
    return aa_map

def analyze_annotated_vcf(annotated_vcf_path, num_examples=5):
    print(f"\n--- DIAGNOSTIC: Analyzing Temporary Annotated VCF: {annotated_vcf_path} ---")
    
    vcf_records_analyzed = 0
    records_with_aa_base = 0
    records_with_aa_dot = 0
    records_with_no_aa_tag = 0 # Should be 0 if annotation worked

    examples_aa_base = []
    examples_aa_dot = []

    try:
        vcf_in_annotated = pysam.VariantFile(annotated_vcf_path, "rb")
    except Exception as e:
        print(f"DIAGNOSTIC_ERROR: Could not open temporary annotated VCF {annotated_vcf_path}: {e}")
        return

    for record in vcf_in_annotated:
        vcf_records_analyzed += 1
        
        aa_tag_value = record.info.get('AA')

        if aa_tag_value is not None:
            if aa_tag_value == '.':
                records_with_aa_dot += 1
                if len(examples_aa_dot) < num_examples:
                    examples_aa_dot.append(
                        f"  AA='.' : VCF {record.chrom}:{record.pos} REF={record.ref} ALT={','.join(str(x) for x in record.alts)} INFO_AA='{aa_tag_value}'"
                    )
            elif aa_tag_value in ['A', 'C', 'G', 'T']:
                records_with_aa_base += 1
                if len(examples_aa_base) < num_examples:
                    examples_aa_base.append(
                        f"  AA=<base>: VCF {record.chrom}:{record.pos} REF={record.ref} ALT={','.join(str(x) for x in record.alts)} INFO_AA='{aa_tag_value}'"
                    )
            else: # Should not happen if annotation script logic is correct
                print(f"  WARNING: Unexpected AA tag value '{aa_tag_value}' in {record.chrom}:{record.pos}")
                records_with_no_aa_tag +=1 # Count as unclassified for now
        else:
            records_with_no_aa_tag += 1
            print(f"  WARNING: No AA tag found in {record.chrom}:{record.pos} in the temp annotated VCF.")

    vcf_in_annotated.close()
    
    print(f"\n--- DIAGNOSTIC: Summary from Reading Temporary Annotated VCF ---")
    print(f"  Total VCF records analyzed from temp file: {vcf_records_analyzed}")
    print(f"  Records with AA=<base> (A/C/G/T): {records_with_aa_base}")
    print(f"  Records with AA='.' (unknown): {records_with_aa_dot}")
    if records_with_no_aa_tag > 0:
        print(f"  WARNING: Records found with NO AA tag or unexpected AA tag: {records_with_no_aa_tag}")

    print(f"\n--- DIAGNOSTIC: Examples from Temporary Annotated VCF (max {num_examples}) ---")
    print(f"  Examples with AA=<base>:")
    if examples_aa_base:
        for ex in examples_aa_base: print(ex)
    else: print("    None found or sampled.")
    
    print(f"\n  Examples with AA='.' (unknown):")
    if examples_aa_dot:
        for ex in examples_aa_dot: print(ex)
    else: print("    None found or sampled.")

def main_diagnostic(vcf_file_path, aa_tsv_file_path, num_examples=5):
    print(f"--- DIAGNOSTIC: Starting VCF vs AA TSV Annotation and Analysis ---")
    print(f"  Input VCF: {vcf_file_path}")
    print(f"  AA TSV: {aa_tsv_file_path}")

    aa_data = load_ancestral_alleles_diagnostic(aa_tsv_file_path)

    if aa_data is None or not aa_data:
        print(f"DIAGNOSTIC_ERROR: No ancestral allele data loaded. Cannot proceed.")
        return

    temp_annotated_vcf_path = "diagnostic_temp_annotated.vcf.gz"
    print(f"  Temporary annotated VCF will be written to: {temp_annotated_vcf_path}")

    vcf_records_processed_for_writing = 0
    vcf_records_matched_in_tsv_for_writing = 0 # For initial pre-annotation count

    try:
        vcf_in = pysam.VariantFile(vcf_file_path, "rb")
        header_out = vcf_in.header.copy()
        if 'AA' not in header_out.info:
            header_out.info.add('AA', 1, 'String', "Ancestral Allele (Diagnostic)")
            print("  Added 'AA' INFO field to header for temp VCF.")
        
        with pysam.VariantFile(temp_annotated_vcf_path, "w", header=header_out) as vcf_out:
            print(f"\n--- DIAGNOSTIC: Annotating and Writing to Temporary VCF ---")
            for record_in in vcf_in:
                vcf_records_processed_for_writing += 1
                record_out = record_in.copy()
                
                key = (record_out.chrom, record_out.pos)
                ancestral_allele_from_tsv = aa_data.get(key)

                # Remove pre-existing AA tag if any, to ensure clean write
                if 'AA' in record_out.info:
                    record_out.info.pop('AA')

                if ancestral_allele_from_tsv is not None:
                    record_out.info['AA'] = ancestral_allele_from_tsv
                    vcf_records_matched_in_tsv_for_writing +=1
                else:
                    record_out.info['AA'] = '.'
                
                vcf_out.write(record_out)
        print(f"  Finished writing temporary VCF. Processed {vcf_records_processed_for_writing} records.")
        print(f"  During writing, {vcf_records_matched_in_tsv_for_writing} records found a match in the TSV.")

    except Exception as e:
        print(f"DIAGNOSTIC_ERROR: Error during VCF processing or writing temp VCF: {e}")
        return
    finally:
        if 'vcf_in' in locals() and vcf_in: # Ensure vcf_in is defined and not None
             vcf_in.close()


    # Now analyze the temporary VCF file we just created
    analyze_annotated_vcf(temp_annotated_vcf_path, num_examples)
    
    # Clean up the temporary file
    try:
        os.remove(temp_annotated_vcf_path)
        # Also remove index if it was created by pysam implicitly (usually for .vcf.gz)
        if os.path.exists(temp_annotated_vcf_path + ".tbi"):
            os.remove(temp_annotated_vcf_path + ".tbi")
        print(f"\n  Cleaned up temporary file: {temp_annotated_vcf_path}")
    except OSError as e:
        print(f"  Warning: Could not remove temporary file {temp_annotated_vcf_path}: {e}")
            
    print(f"\n--- DIAGNOSTIC: Full Analysis Complete ---")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 check_vcf_aa_match.py <input.vcf.gz> <ancestral_alleles_prepared.tsv.gz>")
        sys.exit(1)
    
    vcf_file_arg = sys.argv[1]
    aa_tsv_file_arg = sys.argv[2]
    
    main_diagnostic(vcf_file_arg, aa_tsv_file_arg)
