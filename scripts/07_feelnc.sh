#!/bin/bash
#SBATCH --job-name=07-feelnc
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

###############################################################################
# 07_feelnc.sh - FEELnc coding potential assessment with model training
# 
# This script runs FEELnc (FlExible Extraction of LncRNAs) to identify
# non-coding transcripts using k-mer composition and ORF characteristics.
# It trains a model using known mRNA and known lncRNA from the reference,
# then uses it for prediction with automatic cutoff optimization.
#
# PERFORMANCE OPTIMIZATION:
#   This script creates FASTA files for known mRNA and lncRNA from the
#   reference GTF, then uses these FASTA files as input to FEELnc.
#   This avoids the slow GTF-to-FASTA conversion that FEELnc does internally.
#
# Requirements:
#   - FEELnc (installed via conda or from source)
#   - Reference genome FASTA file
#   - Reference annotation GTF file with gene_biotype annotations
#
# Usage: sbatch 07_feelnc.sh
#
# Outputs:
#   ${OUTDIR}/coding_potential/feelnc/known_mRNA.gtf
#   ${OUTDIR}/coding_potential/feelnc/known_mRNA.fasta
#   ${OUTDIR}/coding_potential/feelnc/known_lncRNA.gtf
#   ${OUTDIR}/coding_potential/feelnc/known_lncRNA.fasta
#   ${OUTDIR}/coding_potential/feelnc/feelnc.lncRNA.fa
#   ${OUTDIR}/coding_potential/feelnc/feelnc.mRNA.fa
#   ${OUTDIR}/coding_potential/feelnc/feelnc.noORF.fa
#   ${OUTDIR}/coding_potential/feelnc/feelnc_RF_summary.txt
#   ${OUTDIR}/coding_potential/feelnc/feelnc_stats.txt
#   ${OUTDIR}/coding_potential/feelnc/feelnc_noncoding.txt
###############################################################################

printf "\n\nStarted Job: 07-feelnc\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results" # Make sure to provide the same path for all scripts

# Reference files for FEELnc model building
# Used to extract known mRNA and known lncRNA for model training
REFERENCE_FASTA="/storage/users/ahinsu/lncPredict/reference/Bos_taurus.chromsomal.fa"
REFERENCE_GTF="/storage/users/ahinsu/lncPredict/reference/Bos_taurus.cleaned.gtf"

# Number of CPU threads to use
THREADS=${SLURM_CPUS_PER_TASK:-12}

# Conda environment name (created from environment.yml)
CONDA_ENV="FEELnc"


# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
module load libs/bioperl-1.6.924.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

# Create output directories
mkdir -p "${OUTDIR}/coding_potential/feelnc"
mkdir -p logs

# Check reference files
if [[ ! -f "${REFERENCE_FASTA}" ]]; then
    echo "ERROR: Reference FASTA not found: ${REFERENCE_FASTA}"
    echo "Please provide a valid reference genome FASTA file"
    exit 1
fi

if [[ ! -f "${REFERENCE_GTF}" ]]; then
    echo "ERROR: Reference GTF not found: ${REFERENCE_GTF}"
    echo "Please provide a valid reference annotation GTF file"
    exit 1
fi

# Check input file - prefer FASTA from step 04
INPUT_FASTA="${OUTDIR}/filtered_transcripts/filtered.fasta"

if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "ERROR: Filtered FASTA not found: ${INPUT_FASTA}"
    echo "Please run 04_filter_transcripts.sh first"
    exit 1
fi

echo "============================================================================"
echo "Step 1: Extracting known mRNA and lncRNA from reference"
echo "============================================================================"
echo ""
echo "Extracting known protein-coding mRNA..."
python3 - "${REFERENCE_GTF}" "${OUTDIR}/coding_potential/feelnc/known_mRNA.gtf" << 'EOF_EXTRACT_MRNA'
import sys
import re

gtf_file = sys.argv[1]
output_gtf = sys.argv[2]

# Filter GTF to extract only protein-coding genes
count = 0
with open(gtf_file, 'r') as f_in, open(output_gtf, 'w') as f_out:
    for line in f_in:
        if line.startswith('#'):
            f_out.write(line)
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        attributes = fields[8]
        
        # Extract only protein_coding transcripts
        if 'gene_biotype "protein_coding"' in attributes or \
           'transcript_biotype "protein_coding"' in attributes:
            f_out.write(line)
            count += 1

print(f"Wrote {count} entries for known mRNA (protein_coding)")

EOF_EXTRACT_MRNA

echo ""
echo "Extracting known lncRNA (lncRNA biotype)..."
python3 - "${REFERENCE_GTF}" "${OUTDIR}/coding_potential/feelnc/known_lncRNA.gtf" << 'EOF_EXTRACT_LNCRNA'
import sys
import re

gtf_file = sys.argv[1]
output_gtf = sys.argv[2]

# Filter GTF to extract lncRNA genes
count = 0
with open(gtf_file, 'r') as f_in, open(output_gtf, 'w') as f_out:
    for line in f_in:
        if line.startswith('#'):
            f_out.write(line)
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        attributes = fields[8]
        
        # Extract only lncRNA biotype
        if 'gene_biotype "lncRNA"' in attributes or \
           'transcript_biotype "lncRNA"' in attributes:
            f_out.write(line)
            count += 1

print(f"Wrote {count} entries for known lncRNA")

EOF_EXTRACT_LNCRNA

echo ""
echo "============================================================================"
echo "Step 2: Converting GTF to FASTA using gffread"
echo "============================================================================"
echo ""
echo "Converting known mRNA GTF to FASTA..."
gffread "${OUTDIR}/coding_potential/feelnc/known_mRNA.gtf" \
    -g "${REFERENCE_FASTA}" \
    -w "${OUTDIR}/coding_potential/feelnc/known_mRNA.fasta"

MRNA_COUNT=$(grep -c "^>" "${OUTDIR}/coding_potential/feelnc/known_mRNA.fasta" || echo "0")
echo "Extracted ${MRNA_COUNT} mRNA sequences"

echo ""
echo "Converting known lncRNA GTF to FASTA..."
gffread "${OUTDIR}/coding_potential/feelnc/known_lncRNA.gtf" \
    -g "${REFERENCE_FASTA}" \
    -w "${OUTDIR}/coding_potential/feelnc/known_lncRNA.fasta"

LNCRNA_COUNT=$(grep -c "^>" "${OUTDIR}/coding_potential/feelnc/known_lncRNA.fasta" || echo "0")
echo "Extracted ${LNCRNA_COUNT} lncRNA sequences"

echo ""
echo "============================================================================"
echo "Step 3: Running FEELnc_codpot.pl with FASTA inputs"
echo "============================================================================"
echo ""
echo "FEELnc Coding Potential Assessment"
echo "This will train a model using known mRNA and lncRNA, then predict on candidates"
echo "FEELnc will automatically find the best cutoff for the newly trained model"
echo ""
echo "Configuration:"
echo "  Candidate transcripts: ${INPUT_FASTA}"
echo "  Known mRNA: ${OUTDIR}/coding_potential/feelnc/known_mRNA.fasta"
echo "  Known lncRNA: ${OUTDIR}/coding_potential/feelnc/known_lncRNA.fasta"
echo "  Threads: ${THREADS}"
echo ""

# Run FEELnc_codpot.pl with FASTA inputs
# -i: input candidate transcripts FASTA
# -a: known mRNA FASTA (protein-coding)
# -l: known lncRNA FASTA (non-coding)
# By using FASTA inputs, we skip the slow GTF-to-FASTA conversion in FEELnc

echo "Running FEELnc_codpot.pl..."
echo ""

# --numtx specifies number of transcripts used for training. Currently it is 10000(mRNA),10000(lncRNA). This is taken from the reference dataset. 
FEELnc_codpot.pl \
    -i "${INPUT_FASTA}" \
    -a "${OUTDIR}/coding_potential/feelnc/known_mRNA.fasta" \
    -l "${OUTDIR}/coding_potential/feelnc/known_lncRNA.fasta" \
    --numtx=10000,10000 \
    --outdir="${OUTDIR}/coding_potential/feelnc" --outname "feelnc" \
    --verbosity=0 --testorftype=4 \
    2>&1 | tee "${OUTDIR}/coding_potential/feelnc/feelnc_codpot.log"

# Check if FEELnc completed successfully by checking for FASTA output files
FEELNC_LNCRNA="${OUTDIR}/coding_potential/feelnc/feelnc.lncRNA.fa"
FEELNC_MRNA="${OUTDIR}/coding_potential/feelnc/feelnc.mRNA.fa"
FEELNC_NOORF="${OUTDIR}/coding_potential/feelnc/feelnc.noORF.fa"
FEELNC_SUMMARY="${OUTDIR}/coding_potential/feelnc/feelnc_RF_summary.txt"

if [[ ! -f "${FEELNC_LNCRNA}" ]] && [[ ! -f "${FEELNC_MRNA}" ]]; then
    echo "ERROR: FEELnc output files not found"
    echo "Expected files: ${FEELNC_LNCRNA}, ${FEELNC_MRNA}"
    echo "Check the log file for errors: ${OUTDIR}/coding_potential/feelnc/feelnc_codpot.log"
    exit 1
fi

echo ""
echo "FEELnc coding potential assessment completed"
echo "Output files:"
echo "  lncRNA: ${FEELNC_LNCRNA}"
echo "  mRNA: ${FEELNC_MRNA}"
echo "  noORF: ${FEELNC_NOORF}"
echo "  Summary: ${FEELNC_SUMMARY}"
echo ""

# Generate statistics and extract non-coding transcripts
echo "Generating FEELnc statistics..."
python3 << EOF
import os
import re

output_dir = "${OUTDIR}/coding_potential/feelnc"

# File paths
lncrna_file = os.path.join(output_dir, "feelnc.lncRNA.fa")
mrna_file = os.path.join(output_dir, "feelnc.mRNA.fa")
noorf_file = os.path.join(output_dir, "feelnc.noORF.fa")
summary_file = os.path.join(output_dir, "feelnc_RF_summary.txt")

# Function to extract transcript IDs from FASTA file
def extract_fasta_ids(fasta_file):
    """Extract transcript IDs from FASTA headers (without '>')"""
    ids = []
    if os.path.exists(fasta_file):
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Extract ID from header, remove '>' and take first field
                    header = line[1:].strip()
                    # Get the ID (first word before space)
                    transcript_id = header.split()[0]
                    ids.append(transcript_id)
    return ids

# Extract transcript IDs from all three files
lncrna_ids = extract_fasta_ids(lncrna_file)
mrna_ids = extract_fasta_ids(mrna_file)
noorf_ids = extract_fasta_ids(noorf_file)

# Combine non-coding transcripts (lncRNA + noORF)
noncoding_ids = lncrna_ids + noorf_ids

# Read cutoff value from summary file
cutoff = None
if os.path.exists(summary_file):
    with open(summary_file, 'r') as f:
        for line in f:
            if 'cutoff' in line.lower() or 'with_cutoff' in line.lower():
                # Extract numeric value from line like "-With_cutoff:	0.3992"
                match = re.search(r'[:\s]+([\d.]+)', line)
                if match:
                    cutoff = float(match.group(1))
                    break

# Calculate statistics
total = len(lncrna_ids) + len(mrna_ids) + len(noorf_ids)
noncoding = len(noncoding_ids)
coding = len(mrna_ids)

# Write statistics file
stats_file = os.path.join(output_dir, "feelnc_stats.txt")
with open(stats_file, 'w') as out:
    out.write("FEELnc Coding Potential Assessment\\n")
    out.write("=" * 60 + "\\n\\n")
    out.write(f"Reference: ${REFERENCE_GTF}\\n")
    out.write(f"Model trained with known mRNA and lncRNA\\n")
    if cutoff is not None:
        out.write(f"Best cutoff determined by FEELnc: {cutoff:.4f}\\n")
    out.write(f"Total transcripts: {total}\\n")
    noncoding_pct = (noncoding/total*100 if total > 0 else 0)
    coding_pct = (coding/total*100 if total > 0 else 0)
    out.write(f"Non-coding (lncRNA + noORF): {noncoding} ({noncoding_pct:.2f}%)\\n")
    out.write(f"  - lncRNA: {len(lncrna_ids)}\\n")
    out.write(f"  - noORF: {len(noorf_ids)}\\n")
    out.write(f"Coding (mRNA): {coding} ({coding_pct:.2f}%)\\n")

# Write non-coding transcript IDs - ensure clean format
noncoding_file = os.path.join(output_dir, "feelnc_noncoding.txt")

# Clean transcript IDs: remove any commas, spaces, or special characters
clean_ids = []
for transcript_id in noncoding_ids:
    seq_id_str = str(transcript_id).strip()
    # Remove '>' prefix if present
    if seq_id_str.startswith('>'):
        seq_id_str = seq_id_str[1:].strip()
    # Split by any delimiter (comma, tab, or space) and take first field
    seq_id_str = re.split(r'[,\t\s]+', seq_id_str)[0]
    clean_ids.append(seq_id_str)

with open(noncoding_file, 'w') as out:
    for transcript_id in clean_ids:
        out.write(f"{transcript_id}\\n")

print(f"FEELnc completed:")
print(f"  Total: {total}")
print(f"  Non-coding (lncRNA + noORF): {noncoding}")
print(f"    - lncRNA: {len(lncrna_ids)}")
print(f"    - noORF: {len(noorf_ids)}")
print(f"  Coding (mRNA): {coding}")
if cutoff is not None:
    print(f"  Cutoff: {cutoff:.4f}")
print(f"Non-coding IDs written to: {noncoding_file}")
print(f"Statistics written to: {stats_file}")
EOF

printf "\n\nCompleted Job: 07-feelnc\n\n"
pwd
date
