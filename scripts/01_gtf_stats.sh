#!/bin/bash
#SBATCH --job-name=01-gtf_stats
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-72

###############################################################################
# 01_gtf_stats.sh - Generate GTF statistics for individual samples
# 
# This script generates statistics for individual StringTie GTF files
# using SLURM array jobs for parallel processing.
#
# Usage:
#   1. Adjust --array parameter above to match number of GTF files
#      Format: --array=1-N%M where N=total files, M=max concurrent jobs
#   2. sbatch 01_gtf_stats.sh
#
# Outputs:
#   ${OUTDIR}/gtf_stats/individual/*.stats.txt
#   ${OUTDIR}/gtf_stats/individual/*.summary.json
###############################################################################

printf "\n\nStarted Job: 01-gtf_stats\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Directory containing StringTie GTF files (assembled transcripts) 
# Will be $OUTDIR/stringtie_assembly for output from script00 or provide the folder with all gtf files
INPUT_GTF_DIR="/storage/users/ahinsu/lncPredict/results/stringtie_assembly"

# GTF file pattern to match (*.gtf for standard files, *.lncRNA.gtf for script 00 output)
GTF_PATTERN="*.gtf"

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results/" # Make sure to provide the same path for all scripts

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
mkdir -p "${OUTDIR}/gtf_stats/individual"
mkdir -p logs

# ============================================================================
# SLURM Array Job Processing
# ============================================================================

# Get list of all GTF files
mapfile -t GTF_FILES < <(find "${INPUT_GTF_DIR}" -maxdepth 1 -name "${GTF_PATTERN}" 2>/dev/null | sort)

# Check if GTF files exist
if [[ ${#GTF_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No GTF files found in ${INPUT_GTF_DIR}"
    exit 1
fi

# Get the GTF file for this array task (SLURM_ARRAY_TASK_ID is 1-indexed)
ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

if [[ ${ARRAY_INDEX} -ge ${#GTF_FILES[@]} ]]; then
    echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} exceeds number of GTF files (${#GTF_FILES[@]})"
    exit 1
fi

GTF_FILE="${GTF_FILES[$ARRAY_INDEX]}"
# Extract sample name - handles both .gtf and .lncRNA.gtf extensions
SAMPLE_NAME=$(basename "${GTF_FILE}")
SAMPLE_NAME="${SAMPLE_NAME%.gtf}"        # Remove .gtf extension
SAMPLE_NAME="${SAMPLE_NAME%.lncRNA}"     # Remove .lncRNA if present (for .lncRNA.gtf files)

echo ""
echo "SLURM Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing GTF file $((ARRAY_INDEX + 1)) of ${#GTF_FILES[@]}: ${SAMPLE_NAME}"
echo "----------------------------------------"

python3 << EOF
import json
from collections import defaultdict

# Parse GTF file
feature_counts = defaultdict(int)
transcript_ids = set()
gene_ids = set()
exon_counts = defaultdict(int)
transcript_lengths = {}
chromosomes = set()

with open('${GTF_FILE}', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        chrom = fields[0]
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        attributes = fields[8]
        
        chromosomes.add(chrom)
        feature_counts[feature] += 1
        
        # Extract transcript and gene IDs
        transcript_id = None
        gene_id = None
        
        for attr in attributes.split(';'):
            attr = attr.strip()
            if attr.startswith('transcript_id'):
                transcript_id = attr.split('"')[1]
            elif attr.startswith('gene_id'):
                gene_id = attr.split('"')[1]
        
        if transcript_id:
            transcript_ids.add(transcript_id)
            
            if feature == 'transcript':
                transcript_lengths[transcript_id] = end - start + 1
            elif feature == 'exon':
                exon_counts[transcript_id] += 1
        
        if gene_id:
            gene_ids.add(gene_id)

# Calculate statistics
stats = {
    'total_features': sum(feature_counts.values()),
    'feature_counts': dict(feature_counts),
    'total_transcripts': len(transcript_ids),
    'total_genes': len(gene_ids),
    'chromosomes': len(chromosomes),
    'avg_exons_per_transcript': sum(exon_counts.values()) / len(transcript_ids) if transcript_ids else 0,
    'avg_transcript_length': sum(transcript_lengths.values()) / len(transcript_lengths) if transcript_lengths else 0,
    'transcripts_with_single_exon': sum(1 for count in exon_counts.values() if count == 1),
    'transcripts_with_multi_exon': sum(1 for count in exon_counts.values() if count > 1)
}

# Write detailed stats
with open('${OUTDIR}/gtf_stats/individual/${SAMPLE_NAME}.stats.txt', 'w') as out:
    out.write(f"GTF Statistics for ${SAMPLE_NAME} (Stage: individual)\n")
    out.write("=" * 60 + "\n\n")
    out.write(f"Total Features: {stats['total_features']}\n")
    out.write(f"Total Transcripts: {stats['total_transcripts']}\n")
    out.write(f"Total Genes: {stats['total_genes']}\n")
    out.write(f"Chromosomes: {stats['chromosomes']}\n")
    out.write(f"Average Exons per Transcript: {stats['avg_exons_per_transcript']:.2f}\n")
    out.write(f"Average Transcript Length: {stats['avg_transcript_length']:.2f} bp\n")
    out.write(f"Single-exon Transcripts: {stats['transcripts_with_single_exon']}\n")
    out.write(f"Multi-exon Transcripts: {stats['transcripts_with_multi_exon']}\n")
    out.write("\nFeature Counts:\n")
    for feature, count in sorted(stats['feature_counts'].items()):
        out.write(f"  {feature}: {count}\n")

# Write JSON summary
with open('${OUTDIR}/gtf_stats/individual/${SAMPLE_NAME}.summary.json', 'w') as out:
    json.dump(stats, out, indent=2)

print(f"GTF statistics generated for ${SAMPLE_NAME}")
EOF

printf "\n\nCompleted Job: 01-gtf_stats\n\n"
pwd
date
