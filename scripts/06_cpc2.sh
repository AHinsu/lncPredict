#!/bin/bash
#SBATCH --job-name=06-cpc2
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

###############################################################################
# 06_cpc2.sh - CPC2 coding potential assessment
# 
# This script runs CPC2 (Coding Potential Calculator 2) to identify
# non-coding transcripts using an SVM-based classifier.
#
# Usage: sbatch 06_cpc2.sh
#
# Outputs:
#   ${OUTDIR}/coding_potential/cpc2/cpc2_output.txt
#   ${OUTDIR}/coding_potential/cpc2/cpc2_stats.txt
#   ${OUTDIR}/coding_potential/cpc2/cpc2_noncoding.txt
###############################################################################

printf "\n\nStarted Job: 06-cpc2\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results" # Make sure to provide the same path for all scripts

# CPC2 threshold (default: 0.5, non-coding < 0.5 < coding)
CPC2_THRESHOLD=0.5

# Conda environment name (created from environment.yml)
CONDA_ENV="lncrna"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

# Create output directory
mkdir -p "${OUTDIR}/coding_potential/cpc2"
mkdir -p logs

# Check input file
INPUT_FASTA="${OUTDIR}/filtered_transcripts/filtered.fasta"
if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "ERROR: Filtered FASTA not found: ${INPUT_FASTA}"
    echo "Please run 04_filter_transcripts.sh first"
    exit 1
fi

echo "CPC2 threshold: ${CPC2_THRESHOLD}"
echo ""

# Run CPC2
echo "Running CPC2..."
CPC2.py \
    -i "${INPUT_FASTA}" \
    -o "${OUTDIR}/coding_potential/cpc2/cpc2_output"


# Generate statistics
echo ""
echo "Generating CPC2 statistics..."
python3 << EOF
import pandas as pd
import re

threshold = ${CPC2_THRESHOLD}

# Read CPC2 results
df = pd.read_csv('${OUTDIR}/coding_potential/cpc2/cpc2_output.txt', sep='\t')

# Ensure proper column names
if 'score' not in df.columns and 'coding_probability' in df.columns:
    df['score'] = df['coding_probability']

# Handle case where score column might be named differently
score_col = 'score' if 'score' in df.columns else 'coding_probability' if 'coding_probability' in df.columns else df.columns[6]
df['score'] = df[score_col]

# Classify transcripts based on score
df['classification'] = df['score'].apply(
    lambda x: 'coding' if x >= threshold else 'noncoding'
)

# Statistics
total = len(df)
noncoding = (df['classification'] == 'noncoding').sum()
coding = (df['classification'] == 'coding').sum()

# Write statistics
with open('${OUTDIR}/coding_potential/cpc2/cpc2_stats.txt', 'w') as out:
    out.write("CPC2 Coding Potential Assessment\n")
    out.write("=" * 60 + "\n\n")
    out.write(f"Threshold: {threshold}\n")
    out.write(f"Total transcripts: {total}\n")
    out.write(f"Non-coding (< {threshold}): {noncoding} ({noncoding/total*100:.2f}%)\n")
    out.write(f"Coding (>= {threshold}): {coding} ({coding/total*100:.2f}%)\n\n")
    
    out.write("Score Distribution:\n")
    out.write(f"  Min:  {df['score'].min():.4f}\n")
    out.write(f"  Max:  {df['score'].max():.4f}\n")
    out.write(f"  Mean: {df['score'].mean():.4f}\n")
    out.write(f"  Median: {df['score'].median():.4f}\n")
    
    if 'label' in df.columns:
        label_counts = df['label'].value_counts()
        out.write("\nPredicted Labels:\n")
        for label, count in label_counts.items():
            out.write(f"  {label}: {count}\n")

# Write non-coding transcript IDs - ensure clean format
id_col = 'ID' if 'ID' in df.columns else df.columns[0]
noncoding_df = df[df['classification'] == 'noncoding']

# Clean transcript IDs: remove any commas, spaces, or special characters
# Extract only the transcript ID (first field if comma/tab/space-separated)
clean_ids = []
for seq_id in noncoding_df[id_col]:
    seq_id_str = str(seq_id).strip()
    # Remove '>' prefix if present
    if seq_id_str.startswith('>'):
        seq_id_str = seq_id_str[1:].strip()
    # Split by any delimiter (comma, tab, or space) and take first field
    seq_id_str = re.split(r'[,\t\s]+', seq_id_str)[0]
    clean_ids.append(seq_id_str)

# Write clean IDs one per line
with open('${OUTDIR}/coding_potential/cpc2/cpc2_noncoding.txt', 'w') as out:
    for transcript_id in clean_ids:
        out.write(f"{transcript_id}\n")

print(f"CPC2 completed: {noncoding} non-coding, {coding} coding (threshold: {threshold})")
EOF

printf "\n\nCompleted Job: 06-cpc2\n\n"
pwd
date
