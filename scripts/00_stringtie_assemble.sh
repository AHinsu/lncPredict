#!/bin/bash
#SBATCH --job-name=00-stringtie_assemble
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%a.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%a.%x.err
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --array=0-134%10


###############################################################################
# 00_stringtie_assemble.sh - Assemble GTF transcripts from BAM files
# 
# This script runs StringTie on aligned BAM files to assemble transcripts
# for each sample using SLURM array jobs for parallel processing.
#
# Usage: 
#   1. Adjust --array parameter above to match number of BAM files
#      Format: --array=0-N%M where N=total files -1, M=max concurrent jobs
#   2. sbatch 00_stringtie_assemble.sh
#
# Outputs:
#   ${OUTDIR}/stringtie_assembly/*.lncRNA.gtf (one per BAM file)
#   ${OUTDIR}/stringtie_assembly/*.log (assembly logs)
###############################################################################

printf "\n\nStarted Job: 00-stringtie_assemble\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Directory containing BAM files (aligned reads)
INPUT_BAM_DIR="/storage/users/ahinsu/lncPredict/bam/files"

# Reference genome annotation GTF file
REFERENCE_GTF="/storage/users/ahinsu/lncPredict/reference/genes.gtf"

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results"  # Make sure to provide the same path for all scripts.

# Number of threads for StringTie
THREADS=20

# Minimum transcript length (nucleotides)
MIN_LENGTH=200

# StringTie strand specificity: --rf (reverse-forward), --fr (forward-reverse), or leave empty for unstranded
STRAND_SPECIFIC=""

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
mkdir -p "${OUTDIR}/stringtie_assembly"

# ============================================================================
# SLURM Array Job Processing
# ============================================================================

# Get list of all BAM files
mapfile -t BAM_FILES < <(find "${INPUT_BAM_DIR}" -maxdepth 1 -name "*.bam" 2>/dev/null | sort)

# Check if BAM files exist
if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found in ${INPUT_BAM_DIR}"
    exit 1
fi

# Get the BAM file for this array task (SLURM_ARRAY_TASK_ID is 1-indexed)
ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

if [[ ${ARRAY_INDEX} -ge ${#BAM_FILES[@]} ]]; then
    echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} exceeds number of BAM files (${#BAM_FILES[@]})"
    exit 1
fi

BAM_FILE="${BAM_FILES[$ARRAY_INDEX]}"
SAMPLE_NAME=$(basename "${BAM_FILE}" .bam)

echo ""
echo "SLURM Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing BAM file $((ARRAY_INDEX + 1)) of ${#BAM_FILES[@]}: ${SAMPLE_NAME}"
echo "----------------------------------------"

# Run StringTie assembly
echo "Running StringTie assembly for ${SAMPLE_NAME}..."
# Note: STRAND_SPECIFIC is only included if set (avoids empty argument)
stringtie "${BAM_FILE}" \
    -G "${REFERENCE_GTF}" \
    ${STRAND_SPECIFIC:+${STRAND_SPECIFIC}} \
    -o "${OUTDIR}/stringtie_assembly/${SAMPLE_NAME}.lncRNA.gtf" \
    -l "${SAMPLE_NAME}" \
    -p "${THREADS}" \
    -m "${MIN_LENGTH}" \
    -v \
    2>&1 | tee "${OUTDIR}/stringtie_assembly/${SAMPLE_NAME}.log"

if [[ $? -eq 0 ]]; then
    echo "Successfully assembled transcripts for ${SAMPLE_NAME}"
    
    # Count transcripts in assembled GTF
    transcript_count=$(grep -c $'\ttranscript\t' "${OUTDIR}/stringtie_assembly/${SAMPLE_NAME}.lncRNA.gtf" || echo "0")
    echo "  Assembled transcripts: ${transcript_count}"
else
    echo "ERROR: StringTie assembly failed for ${SAMPLE_NAME}"
    exit 1
fi

printf "\n\nCompleted Job: 00-stringtie_assemble\n\n"
pwd
date
