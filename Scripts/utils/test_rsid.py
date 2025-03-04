#!/usr/bin/env python3
import os
import sys
import subprocess
import pandas as pd

def process_chromosome_sample(input_file, output_file, dbsnp_vcf, sample_size=5):
    """Process a sample of lines from a chromosome file to retrieve RSIDs."""
    print(f"Processing first {sample_size} variants from {input_file}...")
    
    # Read just the header and the first few lines of the annotation file
    df = pd.read_csv(input_file, sep='\t', nrows=sample_size)
    
    # Initialize new rsid column if it doesn't exist
    if 'rsid' not in df.columns:
        df['rsid'] = '.'
    
    # Get chromosome number from the first row
    chr_num = str(df.iloc[0]['chromosome'])
    # Remove "chr" prefix if present for consistency
    chr_num = chr_num.replace('chr', '')
    
    print(f"Testing retrieval for chromosome {chr_num}, first {len(df)} variants")
    print("-" * 80)
    
    # Track statistics
    total_variants = len(df)
    rsids_found = 0
    
    # First, examine the dbSNP VCF header to understand its format
    print("Checking dbSNP VCF format:")
    cmd = f"tabix -H {dbsnp_vcf} | head -5"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    
    # Process each variant
    for idx, row in df.iterrows():
        print(f"\nVariant {idx+1}/{total_variants}:")
        print(f"  Position: {row['pos']}")
        print(f"  Ref/Alt: {row['ref_vcf']}/{row['alt_vcf']}")
        
        if 'rsid_dbSNP150' in df.columns:
            print(f"  Existing rsid_dbSNP150: {row['rsid_dbSNP150']}")
        
        pos = int(row['pos']) if pd.notna(row['pos']) else None
        if pos is None:
            print("  Skipping: Invalid position")
            continue
            
        ref = row['ref_vcf']
        alt = row['alt_vcf']
        
        # Test query without chr prefix
        cmd1 = f"tabix {dbsnp_vcf} {chr_num}:{pos}-{pos}"
        print(f"  Testing command: {cmd1}")
        result1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True, check=False)
        
        if result1.returncode == 0 and result1.stdout.strip():
            lines = result1.stdout.strip().split('\n')
            print(f"  SUCCESS! Found {len(lines)} entries:")
            for line in lines[:3]:  # Show first 3 entries
                print(f"    {line}")
            if len(lines) > 3:
                print(f"    ... and {len(lines) - 3} more")
        else:
            print("  No results without 'chr' prefix")
            print(f"  Error: {result1.stderr}" if result1.stderr else "  No error message")
            
            # Try with chr prefix
            cmd2 = f"tabix {dbsnp_vcf} chr{chr_num}:{pos}-{pos}"
            print(f"  Testing alternative command: {cmd2}")
            result2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True, check=False)
            
            if result2.returncode == 0 and result2.stdout.strip():
                lines = result2.stdout.strip().split('\n')
                print(f"  SUCCESS with 'chr' prefix! Found {len(lines)} entries:")
                for line in lines[:3]:  # Show first 3 entries
                    print(f"    {line}")
                if len(lines) > 3:
                    print(f"    ... and {len(lines) - 3} more")
            else:
                print("  No results with 'chr' prefix either")
                print(f"  Error: {result2.stderr}" if result2.stderr else "  No error message")
        
        print("-" * 80)
    
    # Now try processing with the actual function
    print("\nNow testing the actual processing function on these variants...")
    rsids_found = process_sample(df, dbsnp_vcf)
    
    # Save results
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\nTest completed. Found {rsids_found}/{total_variants} RSIDs. Saved to {output_file}")
    
    return df


def process_sample(df, dbsnp_vcf):
    """Actual processing function to retrieve RSIDs."""
    rsids_found = 0
    
    for idx, row in df.iterrows():
        if row['rsid'] != '.':
            continue
            
        # Try to use existing rsid_dbSNP150 if available
        if 'rsid_dbSNP150' in df.columns and pd.notna(row['rsid_dbSNP150']) and row['rsid_dbSNP150'] != '.':
            df.at[idx, 'rsid'] = row['rsid_dbSNP150']
            rsids_found += 1
            continue
            
        chr_num = str(row['chromosome']).replace('chr', '')
        pos = int(row['pos']) if pd.notna(row['pos']) else None
        if pos is None:
            continue
            
        ref = row['ref_vcf']
        alt = row['alt_vcf']
        
        # Try without chr prefix first
        cmd = f"tabix {dbsnp_vcf} {chr_num}:{pos}-{pos}"
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
            
            # If first query fails, try with chr prefix
            if result.returncode != 0 or not result.stdout.strip():
                cmd = f"tabix {dbsnp_vcf} chr{chr_num}:{pos}-{pos}"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
            
            if result.returncode == 0 and result.stdout.strip():
                # Parse dbSNP output
                for line in result.stdout.strip().split('\n'):
                    fields = line.split('\t')
                    if len(fields) < 5:
                        continue
                        
                    dbsnp_chr = fields[0].replace('chr', '')
                    try:
                        dbsnp_pos = int(fields[1])
                        dbsnp_rsid = fields[2]
                        dbsnp_ref = fields[3]
                        dbsnp_alt = fields[4]
                    except (ValueError, IndexError):
                        continue
                    
                    # Match by position and alleles
                    if dbsnp_pos == pos and dbsnp_ref == ref and (dbsnp_alt == alt or alt in dbsnp_alt.split(',')):
                        df.at[idx, 'rsid'] = dbsnp_rsid
                        rsids_found += 1
                        break
        except Exception as e:
            print(f"Error processing position {pos} on chr{chr_num}: {e}")
    
    return rsids_found


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input_file> <output_file> <dbsnp_vcf>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    dbsnp_vcf = sys.argv[3]
    
    process_chromosome_sample(input_file, output_file, dbsnp_vcf)


if __name__ == "__main__":
    main()