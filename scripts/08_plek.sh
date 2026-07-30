#!/bin/bash
#SBATCH --job-name=08-plek2
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

###############################################################################
# 08_plek.sh - PLEKv2 coding potential assessment
# 
# This script uses PLEKv2, an upgraded alignment-free method for predicting
# lncRNA and mRNA based on sequence intrinsic features and Coding-Net model.
#
# PLEKv2 uses pre-trained deep learning models and does not require custom
# model building. It provides two models:
#   - Vertebrate (ve): For vertebrate species
#   - Plant (pl): For plant species
#
# PLEKv2 GitHub: https://github.com/AHinsu/plek2
# Citation: Zhou H, Li A, et al. https://pmc.ncbi.nlm.nih.gov/articles/PMC11295476/
#
# Requirements:
#   - Python >= 3.8.5
#   - Conda environment "PLEK2" with dependencies (see installation instructions)
#   - PLEKv2 installed following https://github.com/AHinsu/plek2
#   - PLEKv2 models downloaded and decompressed to conda environment
#
# Usage: sbatch 08_plek.sh
#
# Outputs:
#   ${OUTDIR}/coding_potential/plek2/plek2_scores.txt       (sequence classifications and scores)
#   ${OUTDIR}/coding_potential/plek2/plek2_noncoding.txt    (list of non-coding sequences)
#   ${OUTDIR}/coding_potential/plek2/plek2_stats.txt        (summary statistics)
###############################################################################

printf "\n\nStarted Job: 08-plek2\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results" # Make sure to provide the same path for all scripts

# PLEKv2 model type (user-changeable)
# Options:
#   - "ve" for vertebrate species (default)
#   - "pl" for plant species
PLEK2_MODEL="ve"

# Conda environment name (should contain PLEKv2 and dependencies)
# Installation: Follow instructions in README.md or
# https://github.com/AHinsu/plek2/blob/main/INSTALLATION.md
CONDA_ENV="PLEK2"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

# Create output directories
mkdir -p "${OUTDIR}/coding_potential/plek2"
mkdir -p logs

# Check input file
INPUT_FASTA="${OUTDIR}/filtered_transcripts/filtered.fasta"
if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "ERROR: Filtered FASTA not found: ${INPUT_FASTA}"
    echo "Please run 04_filter_transcripts.sh first"
    exit 1
fi

# Check if PLEK2.py is available in PATH
if ! command -v PLEK2.py &> /dev/null; then
    echo "ERROR: PLEK2.py not found in PATH"
    echo "Please install PLEKv2 following the instructions:"
    echo "  https://github.com/AHinsu/plek2/blob/main/INSTALLATION.md"
    echo "Make sure to copy PLEK2.py to the conda environment bin directory"
    exit 1
fi

# Validate model selection
if [[ "${PLEK2_MODEL}" != "ve" && "${PLEK2_MODEL}" != "pl" ]]; then
    echo "ERROR: Invalid PLEK2_MODEL value: ${PLEK2_MODEL}"
    echo "Must be 've' (vertebrate) or 'pl' (plant)"
    exit 1
fi

echo "============================================================================"
echo "PLEKv2 Coding Potential Assessment"
echo "============================================================================"
echo ""
echo "Configuration:"
echo "  Model: ${PLEK2_MODEL} ($([ "${PLEK2_MODEL}" == "ve" ] && echo "vertebrate" || echo "plant"))"
echo "  Input: ${INPUT_FASTA}"
echo "  Output directory: ${OUTDIR}/coding_potential/plek2"
echo ""

# Run PLEKv2
echo "Running PLEKv2..."
echo "Command: PLEK2.py -i ${INPUT_FASTA} -m ${PLEK2_MODEL} -o ${OUTDIR}/coding_potential/plek2/plek2"
echo ""

PLEK2.py \
    -i "${INPUT_FASTA}" \
    -m "${PLEK2_MODEL}" \
    -o "${OUTDIR}/coding_potential/plek2/plek2"

# Check if PLEKv2 completed successfully
if [[ $? -ne 0 ]]; then
    echo "ERROR: PLEKv2 failed. Check output files for errors:"
    echo "  ${OUTDIR}/coding_potential/plek2/"
    exit 1
fi

echo ""
echo "PLEKv2 completed successfully"
echo ""

# Verify output files exist
SCORES_FILE="${OUTDIR}/coding_potential/plek2/plek2_scores.txt"
NONCODING_FILE="${OUTDIR}/coding_potential/plek2/plek2_noncoding.txt"
STATS_FILE="${OUTDIR}/coding_potential/plek2/plek2_stats.txt"

if [[ ! -f "${SCORES_FILE}" ]]; then
    echo "ERROR: Expected output file not found: ${SCORES_FILE}"
    exit 1
fi

if [[ ! -f "${NONCODING_FILE}" ]]; then
    echo "ERROR: Expected output file not found: ${NONCODING_FILE}"
    exit 1
fi

if [[ ! -f "${STATS_FILE}" ]]; then
    echo "ERROR: Expected output file not found: ${STATS_FILE}"
    exit 1
fi

echo "Output files:"
echo "  Scores: ${SCORES_FILE}"
echo "  Non-coding list: ${NONCODING_FILE}"
echo "  Statistics: ${STATS_FILE}"
echo ""

# Display statistics
echo "PLEKv2 Results Summary:"
cat "${STATS_FILE}"

# ============================================================================
# Ensure uniform output format for script 09
# ============================================================================
echo ""
echo "Post-processing PLEK2 output to ensure uniform format..."
python3 << EOF
import os
import re

output_dir = "${OUTDIR}/coding_potential/plek2"
noncoding_file = os.path.join(output_dir, "plek2_noncoding.txt")

# Read the noncoding file and clean it
if os.path.exists(noncoding_file):
    clean_ids = []
    with open(noncoding_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Remove '>' prefix if present (FASTA format)
            if line.startswith('>'):
                line = line[1:].strip()
            # Split by any delimiter (comma, tab, or space) and take first field
            line = re.split(r'[,\t\s]+', line)[0]
            clean_ids.append(line)
    
    # Write back clean IDs
    with open(noncoding_file, 'w') as f:
        for transcript_id in clean_ids:
            f.write(f"{transcript_id}\n")
    
    print(f"Cleaned {len(clean_ids)} transcript IDs in {noncoding_file}")
else:
    print(f"ERROR: {noncoding_file} not found")
EOF

echo "PLEK2 output post-processing complete"

printf "\n\nCompleted Job: 08-plek2\n\n"
pwd
date
