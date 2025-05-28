import pysam
import sys
import gzip
import argparse # Keep argparse for consistency if you expand it later

def load_rsids_to_exclude(filepath):
    """Loads a set of rsIDs from a file (one rsID per line)."""
    exclude_set = set()
    try:
        with open(filepath, 'r') as f:
            for line in f:
                rsid = line.strip()
                if rsid and not rsid.startswith('#'): # Ignore empty lines and comments
                    exclude_set.add(rsid)
        sys.stderr.write(f"FILTER_VCF_LOADER: Loaded {len(exclude_set)} rsIDs to exclude from {filepath}\n")
    except FileNotFoundError:
        sys.stderr.write(f"FILTER_VCF_LOADER_WARNING: Exclude RSID file not found: {filepath}. No rsIDs will be excluded by this criterion.\n")
    except Exception as e:
        sys.stderr.write(f"FILTER_VCF_LOADER_ERROR: Could not read exclude RSID file {filepath}: {e}\n")
        # Optionally, exit or return None to indicate a critical failure if the file is mandatory but unreadable
    return exclude_set

def main(input_vcf_path, output_vcf_path, exclude_rsid_filepath):
    log_prefix = "FILTER_VCF"
    sys.stderr.write(f"{log_prefix}: Input VCF: {input_vcf_path}\n")
    sys.stderr.write(f"{log_prefix}: Output filtered VCF: {output_vcf_path}\n")
    
    rsids_to_exclude = set()
    if exclude_rsid_filepath: # Check if a path was actually provided
        sys.stderr.write(f"{log_prefix}: Using Exclude RSID file: {exclude_rsid_filepath}\n")
        rsids_to_exclude = load_rsids_to_exclude(exclude_rsid_filepath)
    else:
        sys.stderr.write(f"{log_prefix}: No exclude RSID file provided or path was empty.\n")

    seen_positions = set() 
    seen_rsids = set()     

    try:
        vcf_in = pysam.VariantFile(input_vcf_path, "rb")
    except Exception as e:
        sys.stderr.write(f"{log_prefix}_ERROR: Could not open input VCF {input_vcf_path}: {e}\n")
        sys.exit(1)

    try:
        vcf_out = pysam.VariantFile(output_vcf_path, "w", header=vcf_in.header)
    except Exception as e:
        sys.stderr.write(f"{log_prefix}_ERROR: Could not open output VCF {output_vcf_path}: {e}\n")
        vcf_in.close()
        sys.exit(1)

    variants_processed = 0
    variants_written = 0
    excluded_by_duplicate_pos = 0
    excluded_by_duplicate_rsid = 0
    excluded_by_aa_dot = 0
    excluded_by_rsid_file = 0

    sys.stderr.write(f"{log_prefix}: Starting filtering...\n")
    for record in vcf_in:
        variants_processed += 1
        if variants_processed % 200000 == 0:
            sys.stderr.write(f"  {log_prefix}_PROGRESS: Processed {variants_processed} variants...\n")

        # Filter 1: AA tag being '.'
        aa_tag = record.info.get('AA')
        if aa_tag is None:
            sys.stderr.write(f"  {log_prefix}_WARNING: Variant at {record.chrom}:{record.pos} has no AA tag. Keeping it (not treating as AA='.')\n")
        elif aa_tag == '.':
            excluded_by_aa_dot += 1
            continue

        # Filter 2: rsID from exclude file
        current_rsid = record.id # This can be None
        if current_rsid and current_rsid in rsids_to_exclude:
            excluded_by_rsid_file += 1
            continue

        # Filter 3: Duplicate position (CHROM, POS) - keep first
        position_key = (record.chrom, record.pos)
        if position_key in seen_positions:
            excluded_by_duplicate_pos += 1
            continue 
        seen_positions.add(position_key)

        # Filter 4: Duplicate rsID (if rsID is present and not '.') - keep first
        if current_rsid and current_rsid != "." and current_rsid is not None: 
            if current_rsid in seen_rsids:
                excluded_by_duplicate_rsid += 1
                continue 
            seen_rsids.add(current_rsid)
        
        vcf_out.write(record)
        variants_written += 1

    vcf_in.close()
    vcf_out.close()

    sys.stderr.write(f"\n--- {log_prefix} Filtering Summary ---\n")
    stats_prefix = f"{log_prefix}_STATS"
    sys.stderr.write(f"{stats_prefix}: Total variants processed from input: {variants_processed}\n")
    sys.stderr.write(f"{stats_prefix}: Variants written to output: {variants_written}\n")
    sys.stderr.write(f"{stats_prefix}: Variants excluded due to AA='.' : {excluded_by_aa_dot}\n")
    sys.stderr.write(f"{stats_prefix}: Variants excluded by rsID file list: {excluded_by_rsid_file}\n")
    sys.stderr.write(f"{stats_prefix}: Variants excluded as duplicate position (CHROM+POS): {excluded_by_duplicate_pos}\n")
    sys.stderr.write(f"{stats_prefix}: Variants excluded as duplicate rsID: {excluded_by_duplicate_rsid}\n")
    
    total_excluded_by_filters = excluded_by_aa_dot + excluded_by_rsid_file + excluded_by_duplicate_pos + excluded_by_duplicate_rsid
    sys.stderr.write(f"{stats_prefix}: Total variants explicitly excluded by these filters: {total_excluded_by_filters}\n")
    
    # Sanity check: processed = written + explicitly excluded
    if variants_processed != (variants_written + total_excluded_by_filters):
        # This might happen if a record has no AA tag and is kept, for example.
        # Or if other implicit skips occurred.
        other_skipped = variants_processed - (variants_written + total_excluded_by_filters)
        sys.stderr.write(f"{stats_prefix}_WARNING: Discrepancy in counts. Processed: {variants_processed}, Written: {variants_written}, Explicitly Excluded: {total_excluded_by_filters}. Other unaccounted/kept: {other_skipped}\n")

    sys.stderr.write(f"{log_prefix}: Filtered VCF written to {output_vcf_path}\n")

if __name__ == "__main__":
    log_prefix_main = "FILTER_VCF_MAIN" # Consistent prefix
    parser = argparse.ArgumentParser(
        description="Filter a VCF file: remove duplicate positions/rsIDs, sites with AA='.', and specified rsIDs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter) # Added formatter
    parser.add_argument("input_vcf", help="Input VCF file (must be .vcf.gz and preferably indexed).")
    parser.add_argument("output_vcf", help="Output filtered VCF file (will be .vcf.gz).")
    parser.add_argument("--exclude_rsids", metavar="FILE", help="File containing list of rsIDs to exclude (one per line). If not provided, no rsIDs are excluded by list.", default=None) # Clarified help
    
    args = parser.parse_args()

    try:
        main(args.input_vcf, args.output_vcf, args.exclude_rsids)
        sys.stderr.write(f"{log_prefix_main}: VCF filtering script finished successfully.\n")
    except Exception as e:
        sys.stderr.write(f"{log_prefix_main}_ERROR: An unhandled error occurred: {e}\n")
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
