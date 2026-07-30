#!/bin/bash
#SBATCH --job-name=10-coding_consensus
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

###############################################################################
# 10_coding_consensus.sh - Extract lncRNA candidates based on coding potential consensus
# 
# This script uses results from script 09 to identify and extract lncRNA candidates
# based on a user-defined consensus strategy. It creates GTF, FASTA, and metadata
# files for the selected transcripts, and uses TransDecoder to predict ORFs for 
# downstream analysis.
#
# Usage: sbatch 10_coding_consensus.sh
#
# Strategy values, for number of tools used given by num_tools:
#   0 = Strict: All tools predict non-coding (overlap = num_tools)
#   1 = Moderate: num_tools - 1 tools predict non-coding (overlap = num_tools - 1)
#   2 = Liberal: num_tools - 2 tools predict non-coding (overlap = num_tools - 2)
#   3 = Very Liberal: num_tools - 3 tools predict non-coding (overlap = num_tools - 3)
#
# Outputs:
#   ${OUTDIR}/lncrna_candidates/candidate_transcripts.txt
#   ${OUTDIR}/lncrna_candidates/candidate_transcripts.gtf
#   ${OUTDIR}/lncrna_candidates/candidate_transcripts.fasta
#   ${OUTDIR}/lncrna_candidates/candidate_info.tsv
#   ${OUTDIR}/lncrna_candidates/transdecoder/longest_orfs.pep (for Pfam search in script 11)
#   ${OUTDIR}/lncrna_candidates/extraction_stats.txt
###############################################################################

printf "\n\nStarted Job: 10-coding_consensus\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results" # Make sure to provide the same path for all scripts

# Conda environment name (created from environment.yml)
CONDA_ENV="lncrna"


# Consensus strategy (0-3):
#   0 = Strict: All tools must predict non-coding (overlap = num_tools)
#   1 = Moderate: num_tools - 1 tools predict non-coding (overlap = num_tools - 1)
#   2 = Liberal: num_tools - 2 tools predict non-coding (overlap = num_tools - 2)
#   3 = Very Liberal: num_tools - 3 tools predict non-coding (overlap = num_tools - 3)
STRATEGY=1

# Reference genome FASTA (for extracting sequences with gffread)
REFERENCE_FASTA="/storage/users/ahinsu/lncPredict/reference/Bos_taurus.chromsomal.fa"

# TransDecoder parameters
TRANSDECODER_MIN_PROTEIN_LENGTH=100  # Minimum protein length in amino acids

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

# Create output directories
mkdir -p "${OUTDIR}/lncrna_candidates"
mkdir -p "${OUTDIR}/lncrna_candidates/transdecoder"
mkdir -p logs

# Check input files from script 09
DETAILED_OVERLAP="${OUTDIR}/coding_potential/final_stats/detailed_overlap.tsv"
CLASS_CODE_INFO="${OUTDIR}/filtered_transcripts/class_code_info.tsv"
FILTERED_FASTA="${OUTDIR}/filtered_transcripts/filtered.fasta"
ANNOTATED_GTF="${OUTDIR}/filtered_transcripts/filtered.gtf"

# Check required files
for f in "${DETAILED_OVERLAP}" "${CLASS_CODE_INFO}" "${FILTERED_FASTA}" "${ANNOTATED_GTF}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Required file not found: ${f}"
        echo "Please run previous steps first (especially script 09)"
        exit 1
    fi
done

# Check reference FASTA
if [[ ! -f "${REFERENCE_FASTA}" ]]; then
    echo "ERROR: Reference FASTA not found: ${REFERENCE_FASTA}"
    echo "Please provide a valid reference genome FASTA file"
    exit 1
fi

echo ""
echo "Configuration:"
echo "  Output directory: ${OUTDIR}"
echo "  Strategy: ${STRATEGY}"
echo "  Conda environment: ${CONDA_ENV}"
echo "  Reference FASTA: ${REFERENCE_FASTA}"
echo "  TransDecoder min protein length: ${TRANSDECODER_MIN_PROTEIN_LENGTH} aa"
echo ""

# Generate consensus
echo "Selecting transcripts based on strategy ${STRATEGY}..."
python3 << EOF
import pandas as pd
import json

# Parameters
strategy = ${STRATEGY}
outdir = '${OUTDIR}'

print(f"Strategy: {strategy}")

# Read detailed overlap from script 09
overlap_df = pd.read_csv('${DETAILED_OVERLAP}', sep='\t')
print(f"Total transcripts in detailed_overlap: {len(overlap_df)}")

# Determine total number of tools used
# Count columns that end with '_noncoding'
noncoding_cols = [col for col in overlap_df.columns if col.endswith('_noncoding')]
num_tools = len(noncoding_cols)
print(f"Number of tools used: {num_tools}")
print(f"Tools: {[col.replace('_noncoding', '') for col in noncoding_cols]}")

# Validate strategy and calculate minimum votes
# Strategy 0: overlap = num_tools (strict)
# Strategy 1: overlap = num_tools - 1 (moderate)
# Strategy 2: overlap = num_tools - 2 (liberal)
# Strategy 3: overlap = num_tools - 3 (very liberal)

# Check that num_tools > strategy
if num_tools <= strategy:
    print(f"ERROR: num_tools ({num_tools}) must be greater than strategy ({strategy})")
    print(f"With {num_tools} tools, maximum valid strategy is {num_tools - 1}")
    exit(1)

if strategy == 0:
    # All tools must agree
    min_votes = num_tools
    strategy_name = "Strict (all tools agree)"
elif strategy == 1:
    # All but one tool agree
    min_votes = num_tools - 1
    strategy_name = "Moderate (num_tools - 1 tools agree)"
elif strategy == 2:
    # All but two tools agree
    min_votes = num_tools - 2
    strategy_name = "Liberal (num_tools - 2 tools agree)"
elif strategy == 3:
    # All but three tools agree
    min_votes = num_tools - 3
    strategy_name = "Very Liberal (num_tools - 3 tools agree)"
else:
    print(f"ERROR: Invalid strategy {strategy}. Must be 0, 1, 2, or 3.")
    exit(1)

print(f"Strategy name: {strategy_name}")
print(f"Minimum votes required: {min_votes}/{num_tools}")

# Filter transcripts based on total_votes
selected = overlap_df[overlap_df['total_votes'] >= min_votes].copy()
print(f"Transcripts selected: {len(selected)}")

if len(selected) == 0:
    print("WARNING: No transcripts meet the criteria!")
    # Create empty files
    with open(f'{outdir}/lncrna_candidates/candidate_transcripts.txt', 'w') as f:
        pass
else:
    # Save transcript IDs
    selected['transcript_id'].to_csv(
        f'{outdir}/lncrna_candidates/candidate_transcripts.txt',
        index=False,
        header=False
    )

# Read class code info to get additional metadata
class_df = pd.read_csv('${CLASS_CODE_INFO}', sep='\t')

# Merge selected transcripts with class code info
candidate_info = selected.merge(
    class_df,
    on='transcript_id',
    how='left'
)

# Save candidate info table
candidate_info.to_csv(
    f'{outdir}/lncrna_candidates/candidate_info.tsv',
    sep='\t',
    index=False
)

# Generate statistics
total_transcripts = len(overlap_df)
selected_count = len(selected)
pct_selected = (selected_count / total_transcripts * 100) if total_transcripts > 0 else 0

# Class code distribution
class_counts = selected['class_code'].value_counts()

# Count genes in selected GTF entries
# We'll count this after GTF extraction
# For now, save the counts for later update

# Write statistics
with open(f'{outdir}/lncrna_candidates/extraction_stats.txt', 'w') as out:
    out.write("lncRNA Candidate Extraction Statistics\n")
    out.write("=" * 70 + "\n\n")
    out.write(f"Strategy: {strategy} ({strategy_name})\n")
    out.write(f"Minimum votes required: {min_votes}/{num_tools}\n\n")
    
    out.write(f"Total transcripts analyzed: {total_transcripts}\n")
    out.write(f"Transcripts selected: {selected_count} ({pct_selected:.2f}%)\n\n")
    
    # Placeholder for FASTA and GTF counts - will be updated later
    out.write("Note: Detailed extraction counts will be added after file generation\n\n")
    
    out.write("Class Code Distribution:\n")
    out.write("-" * 70 + "\n")
    for class_code, count in class_counts.items():
        pct = (count / selected_count * 100) if selected_count > 0 else 0
        out.write(f"  {class_code}: {count} ({pct:.2f}%)\n")
    out.write("\n")
    
    # Vote distribution
    out.write("Vote Distribution (selected transcripts):\n")
    out.write("-" * 70 + "\n")
    vote_counts = selected['total_votes'].value_counts().sort_index(ascending=False)
    for votes, count in vote_counts.items():
        pct = (count / selected_count * 100) if selected_count > 0 else 0
        out.write(f"  {votes}/{num_tools} tools: {count} ({pct:.2f}%)\n")

print(f"\nExtraction completed:")
print(f"  Selected {selected_count} transcripts ({pct_selected:.2f}%)")
print(f"  Class codes: {len(class_counts)}")
EOF

# Extract consensus FASTA sequences
echo ""
echo "Extracting candidate FASTA sequences..."
if [[ -f "${OUTDIR}/lncrna_candidates/candidate_transcripts.txt" && -s "${OUTDIR}/lncrna_candidates/candidate_transcripts.txt" ]]; then
    seqkit grep -f "${OUTDIR}/lncrna_candidates/candidate_transcripts.txt" \
        "${FILTERED_FASTA}" > "${OUTDIR}/lncrna_candidates/candidate_transcripts.fasta"
    
    # Count transcripts in FASTA
    FASTA_COUNT=$(grep -c "^>" "${OUTDIR}/lncrna_candidates/candidate_transcripts.fasta" || echo "0")
    echo "  Created: candidate_transcripts.fasta (${FASTA_COUNT} sequences)"
else
    echo "  No candidates selected - creating empty FASTA file"
    touch "${OUTDIR}/lncrna_candidates/candidate_transcripts.fasta"
    FASTA_COUNT=0
fi

# Extract GTF for selected transcripts
echo ""
echo "Extracting candidate GTF..."
python3 << EOF
import sys

import pandas as pd
import json

# Parameters
STRATEGY = ${STRATEGY}
OUTDIR = '${OUTDIR}'
ANNOTATED_GTF = '${ANNOTATED_GTF}'

# Read candidate transcript IDs
candidate_ids = set()
try:
    with open('${OUTDIR}/lncrna_candidates/candidate_transcripts.txt', 'r') as f:
        for line in f:
            tid = line.strip()
            if tid:
                candidate_ids.add(tid)
except FileNotFoundError:
    print("No candidates file found - creating empty GTF")
    candidate_ids = set()

print(f"Extracting GTF for {len(candidate_ids)} transcripts...")

# Read annotated GTF and extract matching transcripts
extracted = 0
extracted_genes = set()
extracted_transcripts = set()

with open('${ANNOTATED_GTF}', 'r') as infile, \
     open('${OUTDIR}/lncrna_candidates/candidate_transcripts.gtf', 'w') as outfile:
    
    for line in infile:
        if line.startswith('#'):
            # Keep comments
            outfile.write(line)
            continue
        
        # Parse GTF line to extract transcript_id and gene_id
        if '\ttranscript\t' in line or '\texon\t' in line:
            # Look for transcript_id in attributes
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                attributes = parts[8]
                # Extract transcript_id and gene_id
                transcript_id = None
                gene_id = None
                
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('transcript_id'):
                        transcript_id = attr.split('"')[1] if '"' in attr else attr.split()[1].strip('"')
                    elif attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1] if '"' in attr else attr.split()[1].strip('"')
                
                if transcript_id and transcript_id in candidate_ids:
                    outfile.write(line)
                    extracted += 1
                    extracted_transcripts.add(transcript_id)
                    if gene_id:
                        extracted_genes.add(gene_id)

print(f"Extracted {extracted} GTF entries")
print(f"Unique genes: {len(extracted_genes)}")
print(f"Unique transcripts: {len(extracted_transcripts)}")

# Save counts for statistics update
with open('${OUTDIR}/lncrna_candidates/gtf_counts.txt', 'w') as f:
    f.write(f"{len(extracted_genes)}\n")
    f.write(f"{len(extracted_transcripts)}\n")
EOF

# Prepare ORF sequences using TransDecoder
echo ""
echo "Running TransDecoder for ORF prediction..."
if [[ -f "${OUTDIR}/lncrna_candidates/candidate_transcripts.fasta" && -s "${OUTDIR}/lncrna_candidates/candidate_transcripts.fasta" ]]; then
    cd "${OUTDIR}/lncrna_candidates/transdecoder" || exit 1
    
    # Extract long ORFs using TransDecoder.LongOrfs
    echo " Extracting long ORFs with TransDecoder.LongOrfs..."
    TransDecoder.LongOrfs \
        -t ../candidate_transcripts.fasta \
        -m ${TRANSDECODER_MIN_PROTEIN_LENGTH} \
        --output_dir . \
        || {
            echo "  ERROR: TransDecoder.LongOrfs failed"
            echo "  Make sure TransDecoder is installed in your conda environment"
            exit 1
        }
    
    # Step 2: Predict likely coding regions using TransDecoder.Predict
    echo "  Moving TransDecoder to previous folder"
        DIR=$(ls -d */) && cp -r ${DIR}* . && rm -r ${DIR} && rm -r __checkpoints_longorfs
    
    
    cd - > /dev/null || exit 1
    
    # The longest_orfs.pep file is created by TransDecoder.LongOrfs and will be used in script 11
    if [[ -f "${OUTDIR}/lncrna_candidates/transdecoder/longest_orfs.pep" ]]; then
        echo "  Created: transdecoder/longest_orfs.pep (for Pfam search in script 11)"
    else
        echo "  WARNING: longest_orfs.pep not created"
    fi
    
    # Count ORFs in the longest_orfs.pep file
    if [[ -f "${OUTDIR}/lncrna_candidates/transdecoder/longest_orfs.pep" ]]; then
        ORF_COUNT=$(grep -c "^>" "${OUTDIR}/lncrna_candidates/transdecoder/longest_orfs.pep" || echo "0")
        echo "  Total ORFs predicted: ${ORF_COUNT}"
    fi
else
    echo "  No candidates - skipping TransDecoder"
    mkdir -p "${OUTDIR}/lncrna_candidates/transdecoder"
    touch "${OUTDIR}/lncrna_candidates/transdecoder/longest_orfs.pep"
fi

# Update extraction statistics with detailed counts
echo ""
echo "Updating extraction statistics..."
python3 << EOF
# Read the counts
fasta_count = ${FASTA_COUNT}
gtf_genes = 0
gtf_transcripts = 0

try:
    with open('${OUTDIR}/lncrna_candidates/gtf_counts.txt', 'r') as f:
        gtf_genes = int(f.readline().strip())
        gtf_transcripts = int(f.readline().strip())
except:
    pass

# Read existing stats file and update it
stats_file = '${OUTDIR}/lncrna_candidates/extraction_stats.txt'
with open(stats_file, 'r') as f:
    lines = f.readlines()

# Find the line to insert after
new_lines = []
for i, line in enumerate(lines):
    new_lines.append(line)
    if 'Transcripts selected:' in line:
        # Insert detailed extraction counts after this line
        new_lines.append('\n')
        new_lines.append('Detailed Extraction Counts:\n')
        new_lines.append('-' * 70 + '\n')
        new_lines.append(f'  Transcripts extracted from FASTA: {fasta_count}\n')
        new_lines.append(f'  Genes extracted from GTF: {gtf_genes}\n')
        new_lines.append(f'  Transcripts extracted from GTF: {gtf_transcripts}\n')
    elif 'Note: Detailed extraction counts will be added after file generation' in line:
        # Skip this placeholder line
        continue

# Write updated stats
with open(stats_file, 'w') as f:
    f.writelines(new_lines)

print(f"Statistics updated:")
print(f"  FASTA transcripts: {fasta_count}")
print(f"  GTF genes: {gtf_genes}")
print(f"  GTF transcripts: {gtf_transcripts}")
EOF

# Clean up temporary file
rm -f "${OUTDIR}/lncrna_candidates/gtf_counts.txt"

echo ""
echo "Output files created in ${OUTDIR}/lncrna_candidates/:"
echo "  - candidate_transcripts.txt (list of transcript IDs)"
echo "  - candidate_transcripts.fasta (sequences)"
echo "  - candidate_transcripts.gtf (annotations)"
echo "  - candidate_info.tsv (metadata with class codes and tool votes)"
echo "  - transdecoder/longest_orfs.pep (ORF predictions for Pfam search in script 11)"
echo "  - extraction_stats.txt (summary statistics)"
echo ""

printf "\n\nCompleted Job: 10-coding_consensus\n\n"
pwd
date
