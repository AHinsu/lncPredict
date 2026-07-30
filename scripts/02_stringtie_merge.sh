#!/bin/bash
#SBATCH --job-name=02-stringtie_merge
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12


###############################################################################
# 02_stringtie_merge.sh - Merge StringTie assemblies
# 
# This script merges multiple StringTie GTF assemblies into a unified
# transcriptome using stringtie --merge.
#
# Usage: sbatch 02_stringtie_merge.sh
#
# Outputs:
#   ${OUTDIR}/stringtie_merge/merged.gtf
#   ${OUTDIR}/stringtie_merge/merged_stats.txt
#   ${OUTDIR}/stringtie_merge/merged_stats.json
#   ${OUTDIR}/stringtie_merge/merge_list.txt
#   ${OUTDIR}/stringtie_merge/merge.log
###############################################################################

printf "\n\nStarted Job: 02-stringtie_merge\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Directory containing StringTie GTF files (assembled transcripts)
# Will be $OUTDIR/stringtie_assembly for output from script00 or provide the folder with all gtf files
INPUT_GTF_DIR="/storage/users/ahinsu/lncPredict/results/stringtie_assembly"

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results/" # Make sure to provide the same path for all scripts

# Conda environment name (created from environment.yml)
CONDA_ENV="lncrna"

# Reference genome annotation GTF file
REFERENCE_GTF="/storage/users/ahinsu/lncPredict/reference/Bos_taurus.cleaned.gtf"

# Minimum transcript length (nucleotides)
MIN_LENGTH=200

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

# Create output directory
mkdir -p "${OUTDIR}/stringtie_merge"
mkdir -p logs

# Create merge list
echo "Creating GTF merge list..."
ls "${INPUT_GTF_DIR}"/*.gtf > "${OUTDIR}/stringtie_merge/merge_list.txt"

echo "GTF files to merge:"
cat "${OUTDIR}/stringtie_merge/merge_list.txt"
echo ""

# Run StringTie merge
echo "Running StringTie merge..."
stringtie --merge \
    -G "${REFERENCE_GTF}" \
    -o "${OUTDIR}/stringtie_merge/merged.gtf" \
    -m "${MIN_LENGTH}" -p 12 \
    "${OUTDIR}/stringtie_merge/merge_list.txt" \
    2>&1 | tee "${OUTDIR}/stringtie_merge/merge.log"

# Generate merged GTF statistics
echo ""
echo "Generating merged GTF statistics..."
python3 << EOF
import json
from collections import defaultdict

# Parse merged GTF
feature_counts = defaultdict(int)
transcript_ids = set()
gene_ids = set()
exon_counts = defaultdict(int)
transcript_lengths = {}

with open('${OUTDIR}/stringtie_merge/merged.gtf', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        attributes = fields[8]
        
        feature_counts[feature] += 1
        
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

stats = {
    'total_transcripts': len(transcript_ids),
    'total_genes': len(gene_ids),
    'avg_exons_per_transcript': sum(exon_counts.values()) / len(transcript_ids) if transcript_ids else 0,
    'avg_transcript_length': sum(transcript_lengths.values()) / len(transcript_lengths) if transcript_lengths else 0
}

# Write statistics
with open('${OUTDIR}/stringtie_merge/merged_stats.txt', 'w') as out:
    out.write("Merged GTF Statistics\n")
    out.write("=" * 60 + "\n\n")
    out.write(f"Total Transcripts: {stats['total_transcripts']}\n")
    out.write(f"Total Genes: {stats['total_genes']}\n")
    out.write(f"Average Exons per Transcript: {stats['avg_exons_per_transcript']:.2f}\n")
    out.write(f"Average Transcript Length: {stats['avg_transcript_length']:.2f} bp\n")

with open('${OUTDIR}/stringtie_merge/merged_stats.json', 'w') as out:
    json.dump(stats, out, indent=2)

print(f"Merged GTF contains {stats['total_transcripts']} transcripts from {stats['total_genes']} genes")
EOF

printf "\n\nCompleted Job: 02-stringtie_merge\n\n"
pwd
date
