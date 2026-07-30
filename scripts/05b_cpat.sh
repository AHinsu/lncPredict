#!/bin/bash
#SBATCH --job-name=05b-cpat
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

###############################################################################
# 05b_cpat.sh - CPAT coding potential assessment with pre-built models
# 
# This script runs CPAT (Coding-Potential Assessment Tool) using pre-built
# models to identify non-coding transcripts based on sequence features.
#
# Use this script when:
#   - Pre-built CPAT models are available for your organism
#   - You want to use standard models (Human, Mouse, Fly, Zebrafish)
#   - You have downloaded models from CPAT website or other sources
#
# CPAT provides pre-built models for:
#   - Human (Homo sapiens)
#   - Mouse (Mus musculus)
#   - Fly (Drosophila melanogaster)
#   - Zebrafish (Danio rerio)
#
# Download pre-built models from:
#   https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/
#
# Usage: sbatch 05b_cpat.sh
#
# Outputs:
#   ${OUTDIR}/coding_potential/cpat/cpat_output.txt
#   ${OUTDIR}/coding_potential/cpat/cpat_stats.txt
#   ${OUTDIR}/coding_potential/cpat/cpat_noncoding.txt
###############################################################################

printf "\n\nStarted Job: 05b-cpat\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results" # Make sure to provide the same path for all scripts

# CPAT threshold (species-specific)
# Human: 0.364, Mouse: 0.44, Fly: 0.39, Zebrafish: 0.38
CPAT_THRESHOLD=0.364

# Path to CPAT hexamer table (optional, leave empty to use default)
# For pre-built models, download from CPAT website:
# Example: /path/to/Human_Hexamer.tsv
CPAT_HEXAMER=""

# Path to CPAT logit model (optional, leave empty to use default)
# For pre-built models, download from CPAT website:
# Example: /path/to/Human_logitModel.RData
CPAT_LOGIT=""

# Conda environment name
CONDA_ENV="lncrna"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

# Create output directory
mkdir -p "${OUTDIR}/coding_potential/cpat"
mkdir -p logs

# Check input file
INPUT_FASTA="${OUTDIR}/filtered_transcripts/filtered.fasta"
if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "ERROR: Filtered FASTA not found: ${INPUT_FASTA}"
    echo "Please run 04_filter_transcripts.sh first"
    exit 1
fi

echo "CPAT threshold: ${CPAT_THRESHOLD}"
echo ""

# Build CPAT command
CPAT_CMD="cpat.py -g ${INPUT_FASTA} --min-orf=50 --top-orf=50 -o ${OUTDIR}/coding_potential/cpat/cpat_output.txt"

# Add hexamer table if provided
if [[ -n "${CPAT_HEXAMER:-}" && -f "${CPAT_HEXAMER}" ]]; then
    echo "Using hexamer table: ${CPAT_HEXAMER}"
    CPAT_CMD="${CPAT_CMD} -x ${CPAT_HEXAMER}"
else
    echo "No hexamer table provided - CPAT will use default"
fi

# Add logit model if provided
if [[ -n "${CPAT_LOGIT:-}" && -f "${CPAT_LOGIT}" ]]; then
    echo "Using logit model: ${CPAT_LOGIT}"
    CPAT_CMD="${CPAT_CMD} -d ${CPAT_LOGIT}"
else
    echo "No logit model provided - CPAT will use default"
fi

echo ""

# Run CPAT
echo "Running CPAT with pre-built models..."
echo "Command: ${CPAT_CMD}"
echo ""

eval "${CPAT_CMD}"

# Generate statistics
echo ""
echo "Generating CPAT statistics..."
python3 << EOF
import pandas as pd
import re

threshold = ${CPAT_THRESHOLD}

# Read CPAT results (using the best ORF probability file)
df = pd.read_csv('${OUTDIR}/coding_potential/cpat/cpat_output.txt.ORF_prob.best.tsv', sep='\t')

# Handle column name variations between CPAT versions
# The file has both 'seq_ID' (transcript ID) and 'ID' (ORF ID) columns
# We want to use 'seq_ID' which contains the transcript IDs (e.g., MSTRG.1.1)
# NOT 'ID' which contains ORF IDs (e.g., MSTRG.1.1_ORF_2)
if 'Coding_prob' in df.columns:
    df.rename(columns={'Coding_prob': 'coding_prob'}, inplace=True)
# Only rename 'ID' to 'seq_ID' if 'seq_ID' doesn't already exist
if 'ID' in df.columns and 'seq_ID' not in df.columns:
    df.rename(columns={'ID': 'seq_ID'}, inplace=True)

# Classify transcripts
df['classification'] = df['coding_prob'].apply(
    lambda x: 'coding' if x >= threshold else 'noncoding'
)

# Statistics
total = len(df)
noncoding = (df['classification'] == 'noncoding').sum()
coding = (df['classification'] == 'coding').sum()

# Write statistics
with open('${OUTDIR}/coding_potential/cpat/cpat_stats.txt', 'w') as out:
    out.write("CPAT Coding Potential Assessment (Pre-built Models)\n")
    out.write("=" * 60 + "\n\n")
    out.write(f"Threshold: {threshold}\n")
    out.write(f"Total transcripts: {total}\n")
    out.write(f"Non-coding (< {threshold}): {noncoding} ({noncoding/total*100:.2f}%)\n")
    out.write(f"Coding (>= {threshold}): {coding} ({coding/total*100:.2f}%)\n\n")
    
    out.write("Score Distribution:\n")
    out.write(f"  Min:  {df['coding_prob'].min():.4f}\n")
    out.write(f"  Max:  {df['coding_prob'].max():.4f}\n")
    out.write(f"  Mean: {df['coding_prob'].mean():.4f}\n")
    out.write(f"  Median: {df['coding_prob'].median():.4f}\n")

# Write non-coding transcript IDs - ensure clean format
id_col = 'seq_ID' if 'seq_ID' in df.columns else df.columns[0]
noncoding_df = df[df['classification'] == 'noncoding']

# Clean transcript IDs: remove any commas, spaces, or special characters
# Extract only the transcript ID (first field if comma/tab/space-separated)
clean_ids = []
for seq_id in noncoding_df[id_col].values:
    seq_id_str = str(seq_id).strip()
    # Remove '>' prefix if present
    if seq_id_str.startswith('>'):
        seq_id_str = seq_id_str[1:].strip()
    # Split by any delimiter (comma, tab, or space) and take first field
    seq_id_str = re.split(r'[,\t\s]+', seq_id_str)[0]
    clean_ids.append(seq_id_str)

# Write clean IDs one per line
with open('${OUTDIR}/coding_potential/cpat/cpat_noncoding.txt', 'w') as out:
    for transcript_id in clean_ids:
        out.write(f"{transcript_id}\n")

print(f"CPAT completed: {noncoding} non-coding, {coding} coding (threshold: {threshold})")
EOF

printf "\n\nCompleted Job: 05b-cpat\n\n"
pwd
date
