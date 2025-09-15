import subprocess
import pandas as pd
import logging
import os

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_average_depth(bam_file):
    """Calculate average read depth from a BAM file using samtools."""
    result = subprocess.run(
        ['samtools', 'depth', bam_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )
    
    if result.returncode != 0:
        raise Exception("Error calculating depth from {}: {}".format(bam_file, result.stderr))
    
    depth_values = [int(line.split()[2]) for line in result.stdout.splitlines() if line]
    return sum(depth_values) / len(depth_values) if depth_values else 0

def parse_snp_positions(snp_file):
    """Parse an SNP file to extract the list of positions."""
    positions = set()
    with open(snp_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) > 0:
                positions.add(parts[0])  # First element is the position
    return positions

def parse_coverage(coverage_str):
    """Parse coverage string and return the coverage value."""
    try:
        return int(coverage_str.split(':')[1])
    except (IndexError, ValueError) as e:
        logging.error("Error parsing coverage from {}: {}".format(coverage_str, e))
        return 0
    
def parse_ratio(cov_reads_str):
    """Extract the Cov and (Reads1 + Reads2) from the given column and return their ratio as total_reads / cov."""
    try:
        parts = cov_reads_str.split(':')
        cov = int(parts[1])
        reads1 = parts[2]
        reads2 = parts[3]

        if reads1 == '-' or reads2 == '-' or reads1 == 'N' or reads2 == 'N':
            return 0

        reads1 = int(reads1)
        reads2 = int(reads2)
        total_reads = reads1 + reads2

        if cov == 0:
            return 0
        
        return total_reads / cov
    except (IndexError, ValueError) as e:
        logging.error("Error parsing ratio from {}: {}".format(cov_reads_str, e))
        return 0
    
def calculate_strand_bias(strand_filter_str):
    """Calculate strand bias from the given strand filter string."""
    try:
        parts = strand_filter_str.split(':')
        R1_plus = int(parts[1])
        R1_minus = int(parts[2])
        R2_plus = int(parts[3])
        R2_minus = int(parts[4])

        total_plus = R1_plus + R2_plus
        total_minus = R1_minus + R2_minus
        total_reads = total_plus + total_minus
        
        if total_reads == 0:
            return True

        plus_proportion = total_plus / total_reads
        minus_proportion = total_minus / total_reads

        return plus_proportion > 0.7 or minus_proportion > 0.7

    except (IndexError, ValueError) as e:
        logging.error("Error calculating strand bias from {}: {}".format(strand_filter_str, e))
        return True
    
def check_freq(freq_str):
    """See if the allele frequency passes our threshold"""
    try:
        parts = freq_str.split(':')
        reads1 = int(parts[2])
        reads2 = int(parts[3])

        # Ensure no division by zero
        total_reads = reads1 + reads2
        if total_reads == 0:
            return False  # Handle cases where no reads are available
        
        # Threshold for allele frequency
        freq_thresh = 0.9
        return (reads1/total_reads > freq_thresh or reads2/total_reads > freq_thresh)  

    except (IndexError, ValueError) as e:
        logging.error(f"Error checking allele frequencies from {freq_str}: {e}")
        return False  # If there was an error, consider it invalid.
    
def diff_alleles(s1_str, s2_str):
    """Check if the alleles from both samples are different"""
    try:
        parts1 = s1_str.split(':')
        parts2 = s2_str.split(':')

        # Get the allele information
        a1 = parts1[0]
        a2 = parts2[0]

        return a1 != a2  # Return True if alleles are different, False if they are the same

    except (IndexError, ValueError) as e:
        logging.error(f"Error checking if alleles from both samples are different: {e}")
        return False  # If there was an error, consider it invalid
    
def is_real_snp(row, avg_depth_sample1, avg_depth_sample2, ratio_threshold):
    """Determine if a row represents a real SNP based on coverage, allele differences, ratio, and strand bias."""
    ref1, var1 = row.get('Ref_sample1'), row.get('Var_sample1')
    ref2, var2 = row.get('Ref_sample2'), row.get('Var_sample2')

    # splitting by the colons and getting the coverage for the row
    cov_sample1 = parse_coverage(row.get('Cons:Cov:Reads1:Reads2:Freq:P-value_sample1', '0:0'))
    cov_sample2 = parse_coverage(row.get('Cons:Cov:Reads1:Reads2:Freq:P-value_sample2', '0:0'))

    # splitting by the colons and then getting the ratio of coverage to number of reads for both alleles
    ratio_sample1 = parse_ratio(row.get('Cons:Cov:Reads1:Reads2:Freq:P-value_sample1', '0:0:0:0:0'))
    ratio_sample2 = parse_ratio(row.get('Cons:Cov:Reads1:Reads2:Freq:P-value_sample2', '0:0:0:0:0'))

    # assuming we dont have strand bias, then if there is another allele called other than the reference, check to see if there is strand bias 
    strand_bias_sample1 = False
    strand_bias_sample2 = False

    # if there are any reads that do not match the reference, check the strand bias. if there is strand bias, true is returned
    if var1 != '.':
        strand_bias_sample1 = calculate_strand_bias(row.get('StrandFilter:R1+:R1-:R2+:R2-:pval_sample1', '0:0:0:0:0'))

    if var2 != '.':
        strand_bias_sample2 = calculate_strand_bias(row.get('StrandFilter:R1+:R1-:R2+:R2-:pval_sample2', '0:0:0:0:0'))

    # return false if there is strand bias or if the ratio of reads to coverage is less than the threshol 
    if (ratio_sample1 < ratio_threshold or ratio_sample2 < ratio_threshold or
        (var1 != '.' and strand_bias_sample1) or (var2 != '.' and strand_bias_sample2)):
        return False
    
    # checking to see if the allele frequency passes our threshold
    #freq1 = check_freq(row.get('Cons:Cov:Reads1:Reads2:Freq:P-value_sample1', '0:0:0:0:0'))
    #freq2 = check_freq(row.get('Cons:Cov:Reads1:Reads2:Freq:P-value_sample2', '0:0:0:0:0'))

    # what if var1 is A and var2 is A but there is only one read of A for var2
    alleles_are_diff = diff_alleles(row.get('Cons:Cov:Reads1:Reads2:Freq:P-value_sample1', '0:0'), row.get('Cons:Cov:Reads1:Reads2:Freq:P-value_sample2', '0:0'))
    
    return (
        var1 != var2 and
        cov_sample1 > 0.2 * avg_depth_sample1 and
        cov_sample2 > 0.2 * avg_depth_sample2 and
    #    freq1 and freq2 and
        alleles_are_diff
    )

bam_file1 = '/path/to/file/' 
bam_file2 = '/path/to/file/' 
cns_file1 = '/path/to/file/'
cns_file2 = '/path/to/file/'
snp_file1 = '/path/to/file/'
snp_file2 = '/path/to/file/'

logging.info("Calculating average depth for sample 1: {}".format(bam_file1))
avg_depth_sample1 = calculate_average_depth(bam_file1)
logging.info("Calculating average depth for sample 2: {}".format(bam_file2))
avg_depth_sample2 = calculate_average_depth(bam_file2)
avg_depth_sample1

snp_positions1 = parse_snp_positions(snp_file1)
snp_positions2 = parse_snp_positions(snp_file2)
all_snp_positions = snp_positions1.union(snp_positions2)
logging.info("Loaded {} unique positions from both SNP files.".format(len(all_snp_positions)))

try:
    df1 = pd.read_csv(cns_file1, sep='\t')
    df2 = pd.read_csv(cns_file2, sep='\t')
except FileNotFoundError as e:
    logging.error("Error reading CNS files: {}".format(e))
    raise
df1.head()

df1_filtered = df1[df1['Position'].astype(str).isin(all_snp_positions)]
df2_filtered = df2[df2['Position'].astype(str).isin(all_snp_positions)]
df1_filtered.head()

if 'Chrom' not in df1.columns or 'Position' not in df1.columns or 'Chrom' not in df2.columns or 'Position' not in df2.columns:
    raise ValueError("One or both CNS files are missing 'Chrom' or 'Position' columns.")
    
merged = pd.merge(df1_filtered, df2_filtered, on=['Chrom', 'Position'], suffixes=('_sample1', '_sample2'))
merged.head()

merged['IsRealSNP'] = merged.apply(
    lambda row: is_real_snp(row, avg_depth_sample1, avg_depth_sample2, ratio_threshold=0.7), axis=1
)
len(merged)
merged.head()

len(merged)
real_snps = merged[merged['IsRealSNP']] 
len(real_snps)
if not real_snps.empty:
    logging.info("Real SNPs found:\n{}".format(real_snps[['Chrom', 'Position', 'Ref_sample1', 'Var_sample1', 'Cons:Cov:Reads1:Reads2:Freq:P-value_sample1', 'Ref_sample2', 'Var_sample2', 'Cons:Cov:Reads1:Reads2:Freq:P-value_sample2']]))
else:
    logging.info("No real SNPs found.")

sample1_id = os.path.basename(bam_file1).split('.')[0]
sample2_id = os.path.basename(bam_file2).split('.')[0]

snp_filename = "/path/to/outputs/{}_vs_{}.snp".format(sample1_id, sample2_id)
snp_data = real_snps[['Position', 'Cons:Cov:Reads1:Reads2:Freq:P-value_sample1', 'Cons:Cov:Reads1:Reads2:Freq:P-value_sample2']]
split_cols = snp_data['Cons:Cov:Reads1:Reads2:Freq:P-value_sample1'].str.split(':', expand=True)
split_cols.columns = ['Cons_sample{}'.format(sample1_id), 'Cov_sample{}'.format(sample1_id), 'Reads1_sample{}'.format(sample1_id), 'Reads2_sample{}'.format(sample1_id), 'Freq_sample{}'.format(sample1_id), 'P-value_sample{}'.format(sample1_id)]
snp_data = snp_data.join(split_cols)
snp_data.head()
split_cols = snp_data['Cons:Cov:Reads1:Reads2:Freq:P-value_sample2'].str.split(':', expand=True)
split_cols.columns = ['Cons_sample{}'.format(sample2_id), 'Cov_sample{}'.format(sample2_id), 'Reads1_sample{}'.format(sample2_id), 'Reads2_sample{}'.format(sample2_id), 'Freq_sample{}'.format(sample2_id), 'P-value_sample{}'.format(sample2_id)]
snp_data = snp_data.join(split_cols)
snp_data.head()
snp_data.to_csv(snp_filename, sep='\t', header=True, index=False)

logging.info("SNP positions and alleles saved to: {}".format(snp_filename))
len(snp_data)
