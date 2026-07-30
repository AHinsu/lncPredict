#!/bin/bash
#SBATCH --job-name=14-create_STAR_reference
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=240G

###############################################################################
# 14_create_STAR_db.sh - Build STAR genome index with unified annotation
#
# This script builds a STAR genome index using the merged reference GTF
# produced by script 12 (12_create_unified_gtf.sh). The resulting index
# can be used for alignment of RNA-seq reads against the augmented reference
# that includes both known annotations and novel lncRNAs.
#
# Prerequisites:
#   - Run 12_create_unified_gtf.sh first to generate merged_filter.gtf
#   - Sufficient RAM (at least 200 GB recommended for large mammalian genomes)
#
# Usage: sbatch 14_create_STAR_db.sh
#
# Outputs:
#   ${OUTDIR}/STAR/                  (STAR genome index directory)
#   ${OUTDIR}/STAR.tar.gz            (compressed archive of the STAR index)
###############################################################################

printf "\n\nStarted Job: 14-create_STAR_reference\n\n"
pwd
date

printf "\n\n"

# ================================
# User configuration
# ================================

# Reference files
GENOME_FASTA="/storage/users/ahinsu/lncPredict/reference/Bos_taurus.chromsomal.fa"
ANNOTATION_GTF="/storage/users/ahinsu/lncPredict/results/final_reference/merged_filter.gtf"

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results" # Make sure to provide the same path for all scripts


THREADS=${SLURM_CPUS_PER_TASK:-20}

# ================================
# Create output directory
# ================================
mkdir -p "$OUTDIR/STAR"

echo "Building STAR genome index..."
echo "Genome: $GENOME_FASTA"
echo "Annotation: $ANNOTATION_GTF"
echo "Threads: $THREADS"
echo "Memory limit: $MEMORY_LIMIT bytes"
echo ""

# ================================
# Run STAR genomeGenerate
# ================================
STAR \
  --runThreadN $THREADS \
  --runMode genomeGenerate \
  --genomeDir "$OUTDIR/STAR" \
  --genomeFastaFiles "$GENOME_FASTA" \
  --sjdbGTFfile "$ANNOTATION_GTF" \
  --sjdbOverhang 100 \
  --genomeSAindexNbases 14

echo ""
echo "STAR index built successfully."

# ================================
# Compress the index
# ================================
echo "Compressing STAR index into tar.gz..."
cd "${OUTDIR}"
tar -czf STAR.tar.gz STAR/

echo "Done."
echo "Output: ${OUTDIR}.tar.gz"


printf "\n\nCompleted Job 14\n\n"
pwd
date
