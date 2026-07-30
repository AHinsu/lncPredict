#!/bin/bash
#SBATCH --job-name=12-create_unified_gtf
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

###############################################################################
# 12_create_unified_gtf.sh - Create unified GTF combining reference and novel lncRNAs
# 
# This script creates a unified GTF file by combining the reference annotation
# with novel lncRNA candidates that have been filtered through the pipeline.
# Supports flexible filtering by transcript IDs and/or class codes.
# It produces two merged GTF files: one with all novel candidates (merged_nofilter.gtf)
# and one with novel transcripts filtered to remove same-strand overlaps with
# protein-coding exons (merged_filter.gtf).
#
# Filtering options (applied in order):
#   1. Filter by transcript IDs from a file (one ID per line, no header)
#   2. Filter by class codes (u, x, i, o, j, etc.)
#   3. Filter out transcripts without strand information ('.' in column 7)
#   4. Filter out same-strand overlaps with protein-coding exons
#   - Set CLASS_CODES="all" or empty to skip class code filtering
#
# NOTE: This script (step 12) should be run BEFORE 13a_setup_blast.sh and
#       13_annotation.sh. Script 13 annotates the filtered candidates
#       produced here.
#
# Prerequisites:
#   - Run 11_orf_pfam.sh first to generate final lncRNA set
#
# Usage: sbatch 12_create_unified_gtf.sh
#
# Outputs:
#   ${OUTDIR}/unified_reference/novel_candidates.gtf
#   ${OUTDIR}/unified_reference/novel_candidates_stats.txt
#   ${OUTDIR}/unified_reference/ref_protein_exons.gtf
#   ${OUTDIR}/unified_reference/novel_no_same_strand_overlap.gtf
#   ${OUTDIR}/unified_reference/novel_no_same_strand_overlap_stats.txt
#   ${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf
#   ${OUTDIR}/unified_reference/final_reference/merged_nofilter_stats.txt
#   ${OUTDIR}/unified_reference/final_reference/merged_filter.gtf
#   ${OUTDIR}/unified_reference/final_reference/merged_filter_stats.txt
#   ${OUTDIR}/unified_reference/novel_candidates_table.tsv
#   ${OUTDIR}/unified_reference/novel_no_same_strand_overlap_table.tsv
###############################################################################

printf "\n\nStarted Job: 12-create_unified_gtf\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results" # Make sure to provide the same path for all scripts

# Reference genome annotation GTF file
REFERENCE_GTF="/storage/users/ahinsu/lncPredict/reference/Bos_taurus.cleaned.gtf"

# Conda environment name (created from environment.yml)
CONDA_ENV="lncrna"

# Class codes to include in the unified GTF
# u: intergenic, x: antisense, i: intronic, o: overlap, j: multi-exon junction
# User can modify this to include different class codes
# Set to "all" or leave empty to skip class code filtering
CLASS_CODES="u,x,i"

# Optional: File with transcript IDs to filter (one ID per line, no header)
# Leave empty to skip transcript ID filtering
# Filtering order: transcript IDs first (if provided), then class codes (if not "all")
TRANSCRIPT_ID_FILE="/storage/users/ahinsu/lncPredict/results/lncRNA_keep.txt"

# Number of threads to use (will use SLURM_CPUS_PER_TASK if available)
THREADS=${SLURM_CPUS_PER_TASK:-4}


# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

# Create output directories
mkdir -p "${OUTDIR}/unified_reference"
mkdir -p "${OUTDIR}/unified_reference/final_reference"
mkdir -p logs

# Check input files
FINAL_LNCRNA_GTF="${OUTDIR}/final_lncRNA/final_lncrna.gtf"

if [[ ! -f "${FINAL_LNCRNA_GTF}" ]]; then
    echo "ERROR: Final lncRNA GTF not found: ${FINAL_LNCRNA_GTF}"
    echo "Please run 11_orf_pfam.sh first"
    exit 1
fi

if [[ ! -f "${REFERENCE_GTF}" ]]; then
    echo "ERROR: Reference GTF not found: ${REFERENCE_GTF}"
    echo "Please provide a valid reference annotation GTF file"
    exit 1
fi

echo ""
echo "Configuration:"
echo "  Final lncRNA GTF: ${FINAL_LNCRNA_GTF}"
echo "  Reference GTF: ${REFERENCE_GTF}"
echo "  Transcript ID file: ${TRANSCRIPT_ID_FILE:-none}"
echo "  Class codes to include: ${CLASS_CODES}"
echo "  Threads: ${THREADS}"
echo ""

# ============================================================================
# Step 1: Filter final_lncrna.gtf by transcript IDs, class_codes, and strand
# Output: unified_reference/novel_candidates.gtf
# Filtering order:
#   1. Transcript IDs (if file provided)
#   2. Class codes (if not "all" or empty)
#   3. Strand: remove transcripts with '.' in strand column (column 7)
# ============================================================================
echo ""
echo "Step 1: Filtering final_lncrna.gtf"
echo "  Transcript ID file: ${TRANSCRIPT_ID_FILE:-none}"
echo "  Class codes: ${CLASS_CODES}"
echo "  Strand filter: removing unstranded transcripts ('.' in column 7)"

python3 << EOF
import sys
import re
from collections import defaultdict

# Load transcript IDs from file if provided
transcript_ids_to_keep = None
transcript_id_file = '${TRANSCRIPT_ID_FILE}'.strip()

if transcript_id_file:
    try:
        print(f"Loading transcript IDs from: {transcript_id_file}")
        transcript_ids_to_keep = set()
        with open(transcript_id_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:  # Skip empty lines
                    transcript_ids_to_keep.add(line)
        print(f"  Loaded {len(transcript_ids_to_keep)} transcript IDs")
    except FileNotFoundError:
        print(f"ERROR: Transcript ID file not found: {transcript_id_file}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read transcript ID file: {e}")
        sys.exit(1)
else:
    print("No transcript ID file provided - using all transcripts")

# Parse class codes
class_codes_str = '${CLASS_CODES}'.replace(' ', '').strip()
skip_class_code_filter = False

if not class_codes_str or class_codes_str.lower() == 'all':
    skip_class_code_filter = True
    print("Class code filtering disabled (using all class codes)")
    class_codes = set()
else:
    class_codes = set(class_codes_str.split(','))
    print(f"Target class codes: {class_codes}")

# Helper to add gene_biotype attribute to a GTF line if not already present
def add_gene_biotype(line, biotype="lncRNApredicted"):
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 9:
        return line
    attributes = parts[8].rstrip()
    # Check if gene_biotype is already an attribute key (not a substring of another key)
    attr_keys = [a.strip().split()[0] for a in attributes.split(';') if a.strip()]
    if 'gene_biotype' not in attr_keys:
        attributes = attributes.rstrip('; ').rstrip(';').rstrip()
        attributes += f'; gene_biotype "{biotype}";'
        parts[8] = attributes
        return '\t'.join(parts) + '\n'
    return line

# Filter transcripts
transcripts_written = set()
genes_written = set()
class_code_transcripts = defaultdict(set)
class_code_genes = defaultdict(set)
unstranded_removed = 0

try:
    with open('${FINAL_LNCRNA_GTF}', 'r') as infile, \
         open('${OUTDIR}/unified_reference/novel_candidates.gtf', 'w') as outfile:
        
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature = parts[2]
            attributes = parts[8]
            
            # Extract transcript_id, gene_id, and class_code
            transcript_id = None
            gene_id = None
            class_code = None
            
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    match = re.search(r'transcript_id "([^"]+)"', attr)
                    if match:
                        transcript_id = match.group(1)
                elif attr.startswith('gene_id'):
                    match = re.search(r'gene_id "([^"]+)"', attr)
                    if match:
                        gene_id = match.group(1)
                elif attr.startswith('class_code'):
                    match = re.search(r'class_code "([^"]+)"', attr)
                    if match:
                        class_code = match.group(1)
            
            # Apply filters for transcripts
            if feature == "transcript":
                keep_transcript = True
                
                # First filter: transcript ID (if provided)
                if transcript_ids_to_keep is not None and transcript_id not in transcript_ids_to_keep:
                    keep_transcript = False
                
                # Second filter: class code (if not skipped)
                if keep_transcript and not skip_class_code_filter and class_code not in class_codes:
                    keep_transcript = False
                
                # Third filter: strand - remove transcripts without strand info ('.' in column 7)
                if keep_transcript and parts[6] == '.':
                    keep_transcript = False
                    unstranded_removed += 1
                
                # Write transcript if it passes all filters
                if keep_transcript:
                    outfile.write(add_gene_biotype(line))
                    if transcript_id:
                        transcripts_written.add(transcript_id)
                    if gene_id:
                        genes_written.add(gene_id)
                    if class_code and transcript_id:
                        class_code_transcripts[class_code].add(transcript_id)
                        if gene_id:
                            class_code_genes[class_code].add(gene_id)
            elif feature == "exon":
                # For exons, only write if their transcript was already written
                if transcript_id in transcripts_written:
                    outfile.write(add_gene_biotype(line))

    print(f"Filtered {len(transcripts_written)} transcripts")
    print(f"Filtered {len(genes_written)} genes")
    print(f"Removed {unstranded_removed} unstranded transcripts (strand='.')")
    
    # Write statistics
    with open('${OUTDIR}/unified_reference/novel_candidates_stats.txt', 'w') as out:
        out.write("Novel Candidates Statistics\n")
        out.write("=" * 70 + "\n\n")
        
        # Document filters applied
        out.write("Filters Applied:\n")
        out.write("-" * 70 + "\n")
        if transcript_id_file:
            out.write(f"  Transcript ID file: {transcript_id_file}\n")
            out.write(f"  Transcript IDs loaded: {len(transcript_ids_to_keep) if transcript_ids_to_keep else 0}\n")
        else:
            out.write("  Transcript ID filter: none (all transcripts)\n")
        
        if skip_class_code_filter:
            out.write("  Class code filter: none (all class codes)\n")
        else:
            out.write(f"  Class codes included: ${CLASS_CODES}\n")
        out.write(f"  Strand filter: unstranded transcripts removed (strand='.'): {unstranded_removed}\n")
        out.write("\n")
        
        out.write(f"Total transcripts: {len(transcripts_written)}\n")
        out.write(f"Total genes: {len(genes_written)}\n\n")
        
        if class_code_transcripts:
            out.write("Breakdown by Class Code:\n")
            out.write("-" * 70 + "\n")
            for class_code in sorted(class_code_transcripts.keys()):
                n_trans = len(class_code_transcripts[class_code])
                n_genes = len(class_code_genes[class_code])
                out.write(f"  Class code '{class_code}':\n")
                out.write(f"    Transcripts: {n_trans}\n")
                out.write(f"    Genes: {n_genes}\n")
    
    print(f"Statistics written to: novel_candidates_stats.txt")
    
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to filter novel candidates"
    exit 1
fi

NOVEL_COUNT=$(grep -c $'\ttranscript\t' "${OUTDIR}/unified_reference/novel_candidates.gtf" 2>/dev/null || echo "0")
echo "  Novel candidates: ${NOVEL_COUNT} transcripts"

# ============================================================================
# Step 2: Merge novel_candidates.gtf with reference GTF
# Output: final_reference/merged_nofilter.gtf
# ============================================================================
echo ""
echo "Step 2: Merging novel_candidates.gtf with reference GTF (no filtering)..."

# Concatenate reference and novel
cat "${REFERENCE_GTF}" "${OUTDIR}/unified_reference/novel_candidates.gtf" > "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.tmp"

# Check if temporary file was created successfully
if [[ ! -s "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.tmp" ]]; then
    echo "ERROR: Failed to create merged temporary file or file is empty"
    exit 1
fi

# Sort and validate GTF using gffread
# gffread will validate the GTF format and output sorted by coordinates
echo "  Sorting and validating GTF using gffread..."

# Define temporary directory for sort operations (used by fallback if needed)
SORT_TMPDIR="${OUTDIR}/unified_reference/final_reference/sort_tmp"

# Try using gffread for sorting and validation
if command -v gffread &> /dev/null; then
    # Redirect gffread errors to log for debugging
    gffread "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.tmp" -o "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf" --sort-alpha -T --keep-comments --keep-genes --force-exons -F --keep-exon-attrs 2> "${OUTDIR}/unified_reference/final_reference/gffread_nofilter.log"
    
    # Check if gffread was successful
    if [[ ! -s "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf" ]]; then
        echo "  WARNING: gffread produced empty output, falling back to native sort"
        echo "  See ${OUTDIR}/unified_reference/final_reference/gffread_nofilter.log for details"
        
        # Create temporary directory for sort operations in the output folder
        mkdir -p "${SORT_TMPDIR}"
        
        # Extract and save header lines (if any)
        grep '^#' "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.tmp" | head -20 > "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.header" 2>/dev/null || touch "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.header"
        
        # Sort the non-header lines using custom temp directory
        grep -v '^#' "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.tmp" | TMPDIR="${SORT_TMPDIR}" sort -k1,1 -k4,4n > "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.sorted"
        
        # Check if sorting was successful
        if [[ ! -s "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.sorted" ]]; then
            echo "ERROR: Sorting failed or produced empty file"
            rm -rf "${SORT_TMPDIR}"
            rm -f "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.tmp" "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.sorted"
            exit 1
        fi
        
        # Combine header and sorted content
        cat "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.sorted" > "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf"
        
        # Clean up temporary files and directory
        rm -f "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.sorted"
        rm -rf "${SORT_TMPDIR}"
    else
        echo "  ✓ Successfully sorted with gffread"
    fi
else
    echo "  WARNING: gffread not found, using native sort with custom temp directory"
    
    # Create temporary directory for sort operations in the output folder
    mkdir -p "${SORT_TMPDIR}"
    
    # Extract and save header lines (if any)
    grep '^#' "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.tmp" | head -20 > "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.header" 2>/dev/null || touch "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.header"
    
    # Sort the non-header lines using custom temp directory
    grep -v '^#' "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.tmp" | TMPDIR="${SORT_TMPDIR}" sort -k1,1 -k4,4n > "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.sorted"
    
    # Check if sorting was successful
    if [[ ! -s "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.sorted" ]]; then
        echo "ERROR: Sorting failed or produced empty file"
        rm -rf "${SORT_TMPDIR}"
        rm -f "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.tmp" "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.sorted"
        exit 1
    fi
    
    # Combine header and sorted content
    cat "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.sorted" > "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf"
    
    # Clean up temporary files and directory
    rm -f "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.sorted"
    rm -rf "${SORT_TMPDIR}"
fi

# Clean up temporary file
rm -f "${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf.tmp"

echo "  Merged GTF created: merged_nofilter.gtf"

# Validate and generate stats
echo "  Validating merged_nofilter.gtf..."

python3 << EOF
import sys
import re
from collections import defaultdict

def count_gtf_features(gtf_file):
    """Count genes and transcripts in a GTF file"""
    transcripts = set()
    genes = set()
    
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature = parts[2]
                attributes = parts[8]
                
                if feature != "transcript":
                    continue
                
                transcript_id = None
                gene_id = None
                
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('transcript_id'):
                        match = re.search(r'transcript_id "([^"]+)"', attr)
                        if match:
                            transcript_id = match.group(1)
                    elif attr.startswith('gene_id'):
                        match = re.search(r'gene_id "([^"]+)"', attr)
                        if match:
                            gene_id = match.group(1)
                
                if transcript_id:
                    transcripts.add(transcript_id)
                if gene_id:
                    genes.add(gene_id)
    
    except Exception as e:
        print(f"Error counting features: {e}")
        return 0, 0
    
    return len(transcripts), len(genes)

# Count features in input files
print("Counting features in reference GTF...")
ref_trans, ref_genes = count_gtf_features('${REFERENCE_GTF}')
print(f"  Reference: {ref_trans} transcripts, {ref_genes} genes")

print("Counting features in novel_candidates.gtf...")
novel_trans, novel_genes = count_gtf_features('${OUTDIR}/unified_reference/novel_candidates.gtf')
print(f"  Novel candidates: {novel_trans} transcripts, {novel_genes} genes")

print("Counting features in merged_nofilter.gtf...")
merged_trans, merged_genes = count_gtf_features('${OUTDIR}/unified_reference/final_reference/merged_nofilter.gtf')
print(f"  Merged (no filter): {merged_trans} transcripts, {merged_genes} genes")

# Write statistics
with open('${OUTDIR}/unified_reference/final_reference/merged_nofilter_stats.txt', 'w') as out:
    out.write("Merged GTF Statistics (No Overlap Filtering)\n")
    out.write("=" * 70 + "\n\n")
    
    out.write("Input Files:\n")
    out.write("-" * 70 + "\n")
    out.write(f"Reference GTF:\n")
    out.write(f"  Transcripts: {ref_trans}\n")
    out.write(f"  Genes: {ref_genes}\n\n")
    
    out.write(f"Novel Candidates GTF:\n")
    out.write(f"  Transcripts: {novel_trans}\n")
    out.write(f"  Genes: {novel_genes}\n\n")
    
    out.write("Output File (merged_nofilter.gtf):\n")
    out.write("-" * 70 + "\n")
    out.write(f"  Transcripts: {merged_trans}\n")
    out.write(f"  Genes: {merged_genes}\n\n")
    
    out.write("Validation:\n")
    out.write("-" * 70 + "\n")
    expected = ref_trans + novel_trans
    if merged_trans == expected:
        out.write(f"  PASS: Output contains expected number of transcripts\n")
        out.write(f"  ({ref_trans} + {novel_trans} = {expected})\n")
    else:
        out.write(f"  WARNING: Unexpected transcript count\n")
        out.write(f"  Expected: {expected}, Found: {merged_trans}\n")

print("Statistics written to: merged_nofilter_stats.txt")

# Exit with error if no transcripts found in merged file
if merged_trans == 0:
    print("ERROR: No transcripts found in merged_nofilter.gtf!")
    print("This indicates a problem with the merge or sort process.")
    sys.exit(1)

# Exit with error if transcript count doesn't match expected
if merged_trans != expected:
    print(f"ERROR: Transcript count validation failed!")
    print(f"  Expected: {expected} (ref: {ref_trans} + novel: {novel_trans})")
    print(f"  Found: {merged_trans}")
    print("This indicates transcripts were lost during merge or sort.")
    sys.exit(1)

print(f"✓ Validation PASSED: {merged_trans} transcripts found as expected")
EOF

# ============================================================================
# Step 3: Extract reference protein coding exons for filtering
# ============================================================================
echo ""
echo "Step 3: Extracting reference protein coding exons..."

python3 << EOF
import re

exon_count = 0

try:
    with open('${REFERENCE_GTF}', 'r') as infile, \
         open('${OUTDIR}/unified_reference/ref_protein_exons.gtf', 'w') as outfile:
        
        for line in infile:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature = parts[2]
            attributes = parts[8]
            
            # Check for protein_coding gene_biotype or gene_type
            if feature == "exon":
                if 'gene_biotype "protein_coding"' in attributes or \
                   'gene_type "protein_coding"' in attributes:
                    outfile.write(line)
                    exon_count += 1

    print(f"Extracted {exon_count} protein-coding exons from reference")
    
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()
EOF

# ============================================================================
# Step 4: Filter novel_candidates.gtf to remove same-strand overlaps
# Output: unified_reference/novel_no_same_strand_overlap.gtf
# ============================================================================
echo ""
echo "Step 4: Filtering novel transcripts overlapping protein-coding exons (same strand)..."

if [[ ! -s "${OUTDIR}/unified_reference/novel_candidates.gtf" ]]; then
    echo "  WARNING: No novel candidates found. Creating empty output."
    touch "${OUTDIR}/unified_reference/novel_no_same_strand_overlap.gtf"
else
    # 1. Identify transcripts that overlap protein-coding exons
    bedtools intersect -s -a "${OUTDIR}/unified_reference/novel_candidates.gtf" \
        -b "${OUTDIR}/unified_reference/ref_protein_exons.gtf" \
        -wa -u \
        | awk '$3=="transcript"' \
        | sed 's/.*transcript_id "//;s/".*//' \
        | sort -u \
        > "${OUTDIR}/unified_reference/overlapping_transcripts.txt"
        
    # 2. Remove those transcripts AND all their exons
    awk '
        BEGIN {
            while ((getline < "'"${OUTDIR}/unified_reference/overlapping_transcripts.txt"'" ) > 0) {
                remove[$1] = 1
            }
        }
        {
            # extract transcript_id
            match($0, /transcript_id "([^"]+)"/, m)
            tid = m[1]

            # keep only if transcript_id is NOT in remove list
            if (!(tid in remove)) print $0
        }
    ' "${OUTDIR}/unified_reference/novel_candidates.gtf" > "${OUTDIR}/unified_reference/novel_no_same_strand_overlap.gtf"
    
    # Count filtered transcripts
    BEFORE=$(grep -c $'\ttranscript\t' "${OUTDIR}/unified_reference/novel_candidates.gtf" || echo "0")
    AFTER=$(grep -c $'\ttranscript\t' "${OUTDIR}/unified_reference/novel_no_same_strand_overlap.gtf" || echo "0")
    REMOVED=$((BEFORE - AFTER))
    
    echo "  Transcripts before filtering: ${BEFORE}"
    echo "  Transcripts after filtering: ${AFTER}"
    echo "  Transcripts removed: ${REMOVED}"
fi

# Generate stats for filtered novel transcripts
python3 << EOF
import re
from collections import defaultdict

def count_gtf_features_with_class(gtf_file):
    """Count genes, transcripts, and class codes in a GTF file"""
    transcripts = set()
    genes = set()
    class_code_transcripts = defaultdict(set)
    class_code_genes = defaultdict(set)
    
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature = parts[2]
                attributes = parts[8]
                
                if feature != "transcript":
                    continue
                
                transcript_id = None
                gene_id = None
                class_code = None
                
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('transcript_id'):
                        match = re.search(r'transcript_id "([^"]+)"', attr)
                        if match:
                            transcript_id = match.group(1)
                    elif attr.startswith('gene_id'):
                        match = re.search(r'gene_id "([^"]+)"', attr)
                        if match:
                            gene_id = match.group(1)
                    elif attr.startswith('class_code'):
                        match = re.search(r'class_code "([^"]+)"', attr)
                        if match:
                            class_code = match.group(1)
                
                if transcript_id:
                    transcripts.add(transcript_id)
                if gene_id:
                    genes.add(gene_id)
                if class_code and transcript_id:
                    class_code_transcripts[class_code].add(transcript_id)
                    if gene_id:
                        class_code_genes[class_code].add(gene_id)
    
    except Exception as e:
        print(f"Error counting features: {e}")
        return 0, 0, {}, {}
    
    return len(transcripts), len(genes), class_code_transcripts, class_code_genes

# Count before and after filtering
print("Counting features before filtering...")
before_trans, before_genes, before_cc_trans, before_cc_genes = count_gtf_features_with_class('${OUTDIR}/unified_reference/novel_candidates.gtf')

print("Counting features after filtering...")
after_trans, after_genes, after_cc_trans, after_cc_genes = count_gtf_features_with_class('${OUTDIR}/unified_reference/novel_no_same_strand_overlap.gtf')

# Write statistics
with open('${OUTDIR}/unified_reference/novel_no_same_strand_overlap_stats.txt', 'w') as out:
    out.write("Novel Transcripts - Same-Strand Overlap Filtering Statistics\n")
    out.write("=" * 70 + "\n\n")
    
    out.write("Before Filtering (novel_candidates.gtf):\n")
    out.write("-" * 70 + "\n")
    out.write(f"  Transcripts: {before_trans}\n")
    out.write(f"  Genes: {before_genes}\n\n")
    
    out.write("After Filtering (novel_no_same_strand_overlap.gtf):\n")
    out.write("-" * 70 + "\n")
    out.write(f"  Transcripts: {after_trans}\n")
    out.write(f"  Genes: {after_genes}\n")
    
    if after_cc_trans:
        out.write("\n  Breakdown by Class Code:\n")
        for class_code in sorted(after_cc_trans.keys()):
            n_trans = len(after_cc_trans[class_code])
            n_genes = len(after_cc_genes[class_code])
            out.write(f"    Class code '{class_code}': {n_trans} transcripts, {n_genes} genes\n")
    
    out.write("\n")
    out.write("Filtering Summary:\n")
    out.write("-" * 70 + "\n")
    out.write(f"  Transcripts removed: {before_trans - after_trans}\n")
    out.write(f"  Genes removed: {before_genes - after_genes}\n")
    if before_trans > 0:
        out.write(f"  Retention rate: {after_trans/before_trans*100:.2f}%\n")
    out.write("\n")
    out.write("Note: Filtered transcripts overlapped protein-coding exons on the same strand.\n")

print(f"Statistics written to: novel_no_same_strand_overlap_stats.txt")
print(f"  Before: {before_trans} transcripts, {before_genes} genes")
print(f"  After: {after_trans} transcripts, {after_genes} genes")
print(f"  Removed: {before_trans - after_trans} transcripts")
EOF

# ============================================================================
# Step 5: Create combined GTF with filtered novel transcripts
# Output: final_reference/merged_filter.gtf
# ============================================================================
echo ""
echo "Step 5: Creating combined GTF with filtered novel transcripts..."

# Concatenate reference and filtered novel
cat "${REFERENCE_GTF}" "${OUTDIR}/unified_reference/novel_no_same_strand_overlap.gtf" > "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.tmp"

# Check if temporary file was created successfully
if [[ ! -s "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.tmp" ]]; then
    echo "ERROR: Failed to create merged temporary file or file is empty"
    exit 1
fi

# Sort and validate GTF using gffread
# gffread will validate the GTF format and output sorted by coordinates
echo "  Sorting and validating GTF using gffread..."

# Define temporary directory for sort operations (used by fallback if needed)
SORT_TMPDIR="${OUTDIR}/unified_reference/final_reference/sort_tmp"

# Try using gffread for sorting and validation
if command -v gffread &> /dev/null; then
    # Redirect gffread errors to log for debugging
    gffread "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.tmp" -o "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf" --sort-alpha -T --keep-comments --keep-genes --force-exons -F --keep-exon-attrs 2> "${OUTDIR}/unified_reference/final_reference/gffread_filter.log"
    
    # Check if gffread was successful
    if [[ ! -s "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf" ]]; then
        echo "  WARNING: gffread produced empty output, falling back to native sort"
        echo "  See ${OUTDIR}/unified_reference/final_reference/gffread_filter.log for details"
        
        # Create temporary directory for sort operations in the output folder
        mkdir -p "${SORT_TMPDIR}"
        
        # Extract and save header lines (if any)
        grep '^#' "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.tmp" | head -20 > "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.header" 2>/dev/null || touch "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.header"
        
        # Sort the non-header lines using custom temp directory
        grep -v '^#' "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.tmp" | TMPDIR="${SORT_TMPDIR}" sort -k1,1 -k4,4n > "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.sorted"
        
        # Check if sorting was successful
        if [[ ! -s "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.sorted" ]]; then
            echo "ERROR: Sorting failed or produced empty file"
            rm -rf "${SORT_TMPDIR}"
            rm -f "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.tmp" "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.sorted"
            exit 1
        fi
        
        # Combine header and sorted content
        cat "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.sorted" > "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf"
        
        # Clean up temporary files and directory
        rm -f "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.sorted"
        rm -rf "${SORT_TMPDIR}"
    else
        echo "  ✓ Successfully sorted with gffread"
    fi
else
    echo "  WARNING: gffread not found, using native sort with custom temp directory"
    
    # Create temporary directory for sort operations in the output folder
    mkdir -p "${SORT_TMPDIR}"
    
    # Extract and save header lines (if any)
    grep '^#' "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.tmp" | head -20 > "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.header" 2>/dev/null || touch "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.header"
    
    # Sort the non-header lines using custom temp directory
    grep -v '^#' "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.tmp" | TMPDIR="${SORT_TMPDIR}" sort -k1,1 -k4,4n > "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.sorted"
    
    # Check if sorting was successful
    if [[ ! -s "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.sorted" ]]; then
        echo "ERROR: Sorting failed or produced empty file"
        rm -rf "${SORT_TMPDIR}"
        rm -f "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.tmp" "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.sorted"
        exit 1
    fi
    
    # Combine header and sorted content
    cat "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.sorted" > "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf"
    
    # Clean up temporary files and directory
    rm -f "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.header" "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.sorted"
    rm -rf "${SORT_TMPDIR}"
fi

# Clean up temporary file
rm -f "${OUTDIR}/unified_reference/final_reference/merged_filter.gtf.tmp"

echo "  Merged GTF created: merged_filter.gtf"

# Validate and generate stats
echo "  Validating merged_filter.gtf..."

python3 << EOF
import sys
import re

def count_gtf_features(gtf_file):
    """Count genes and transcripts in a GTF file"""
    transcripts = set()
    genes = set()
    
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature = parts[2]
                attributes = parts[8]
                
                if feature != "transcript":
                    continue
                
                transcript_id = None
                gene_id = None
                
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('transcript_id'):
                        match = re.search(r'transcript_id "([^"]+)"', attr)
                        if match:
                            transcript_id = match.group(1)
                    elif attr.startswith('gene_id'):
                        match = re.search(r'gene_id "([^"]+)"', attr)
                        if match:
                            gene_id = match.group(1)
                
                if transcript_id:
                    transcripts.add(transcript_id)
                if gene_id:
                    genes.add(gene_id)
    
    except Exception as e:
        print(f"Error counting features: {e}")
        return 0, 0
    
    return len(transcripts), len(genes)

# Count features in input files
print("Counting features in reference GTF...")
ref_trans, ref_genes = count_gtf_features('${REFERENCE_GTF}')
print(f"  Reference: {ref_trans} transcripts, {ref_genes} genes")

print("Counting features in novel_no_same_strand_overlap.gtf...")
novel_trans, novel_genes = count_gtf_features('${OUTDIR}/unified_reference/novel_no_same_strand_overlap.gtf')
print(f"  Novel (filtered): {novel_trans} transcripts, {novel_genes} genes")

print("Counting features in merged_filter.gtf...")
merged_trans, merged_genes = count_gtf_features('${OUTDIR}/unified_reference/final_reference/merged_filter.gtf')
print(f"  Merged (with filter): {merged_trans} transcripts, {merged_genes} genes")

# Write statistics
with open('${OUTDIR}/unified_reference/final_reference/merged_filter_stats.txt', 'w') as out:
    out.write("Merged GTF Statistics (With Overlap Filtering)\n")
    out.write("=" * 70 + "\n\n")
    
    out.write("Input Files:\n")
    out.write("-" * 70 + "\n")
    out.write(f"Reference GTF:\n")
    out.write(f"  Transcripts: {ref_trans}\n")
    out.write(f"  Genes: {ref_genes}\n\n")
    
    out.write(f"Novel Transcripts (after same-strand overlap filtering):\n")
    out.write(f"  Transcripts: {novel_trans}\n")
    out.write(f"  Genes: {novel_genes}\n\n")
    
    out.write("Output File (merged_filter.gtf):\n")
    out.write("-" * 70 + "\n")
    out.write(f"  Transcripts: {merged_trans}\n")
    out.write(f"  Genes: {merged_genes}\n\n")
    
    out.write("Validation:\n")
    out.write("-" * 70 + "\n")
    expected = ref_trans + novel_trans
    if merged_trans == expected:
        out.write(f"  PASS: Output contains expected number of transcripts\n")
        out.write(f"  ({ref_trans} + {novel_trans} = {expected})\n")
    else:
        out.write(f"  WARNING: Unexpected transcript count\n")
        out.write(f"  Expected: {expected}, Found: {merged_trans}\n")

print("Statistics written to: merged_filter_stats.txt")

# Exit with error if no transcripts found in merged file
if merged_trans == 0:
    print("ERROR: No transcripts found in merged_filter.gtf!")
    print("This indicates a problem with the merge or sort process.")
    sys.exit(1)

# Exit with error if transcript count doesn't match expected
if merged_trans != expected:
    print(f"ERROR: Transcript count validation failed!")
    print(f"  Expected: {expected} (ref: {ref_trans} + novel: {novel_trans})")
    print(f"  Found: {merged_trans}")
    print("This indicates transcripts were lost during merge or sort.")
    sys.exit(1)

print(f"✓ Validation PASSED: {merged_trans} transcripts found as expected")
EOF

# ============================================================================
# Step 6: Create gene-level TSV tables for novel_candidates.gtf and
#         novel_no_same_strand_overlap.gtf
# Columns: gene_id, chromosome, start, end, class_code, transcript_count,
#          exon_count, xloc, cmp_ref
# ============================================================================
echo ""
echo "Step 6: Creating gene-level TSV tables..."

python3 << EOF
import re
import sys
from collections import defaultdict

def gtf_to_gene_table(gtf_file, output_file):
    """Parse a GTF file and write a gene-level TSV summary table."""
    gene_data = {}

    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                chrom   = parts[0]
                feature = parts[2]
                start   = int(parts[3])
                end     = int(parts[4])
                attrs   = parts[8]

                gene_id       = None
                transcript_id = None
                class_code    = None
                xloc          = None
                cmp_ref       = None

                for attr in attrs.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        m = re.search(r'gene_id "([^"]+)"', attr)
                        if m: gene_id = m.group(1)
                    elif attr.startswith('transcript_id'):
                        m = re.search(r'transcript_id "([^"]+)"', attr)
                        if m: transcript_id = m.group(1)
                    elif attr.startswith('class_code'):
                        m = re.search(r'class_code "([^"]+)"', attr)
                        if m: class_code = m.group(1)
                    elif attr.startswith('xloc'):
                        m = re.search(r'xloc "([^"]+)"', attr)
                        if m: xloc = m.group(1)
                    elif attr.startswith('cmp_ref'):
                        m = re.search(r'cmp_ref "([^"]+)"', attr)
                        if m: cmp_ref = m.group(1)

                if not gene_id:
                    continue

                # Initialize gene record on first encounter
                if gene_id not in gene_data:
                    gene_data[gene_id] = {
                        'chrom':       chrom,
                        'start':       start,
                        'end':         end,
                        'class_code':  '.',
                        'transcripts': set(),
                        'exon_count':  0,
                        'xloc':        '.',
                        'cmp_ref':     '.',
                    }

                g = gene_data[gene_id]
                if start < g['start']:
                    g['start'] = start
                if end > g['end']:
                    g['end'] = end

                if feature == 'transcript':
                    if transcript_id:
                        g['transcripts'].add(transcript_id)
                    if class_code and g['class_code'] == '.':
                        g['class_code'] = class_code
                    if xloc and g['xloc'] == '.':
                        g['xloc'] = xloc
                    if cmp_ref and g['cmp_ref'] == '.':
                        g['cmp_ref'] = cmp_ref
                elif feature == 'exon':
                    g['exon_count'] += 1

    except FileNotFoundError:
        print(f"ERROR: File not found: {gtf_file}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR reading {gtf_file}: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    try:
        with open(output_file, 'w') as out:
            out.write('gene_id\tchromosome\tstart\tend\tclass_code\ttranscript_count\texon_count\txloc\tcmp_ref\n')
            for gene_id in sorted(gene_data.keys()):
                g = gene_data[gene_id]
                out.write(
                    f"{gene_id}\t{g['chrom']}\t{g['start']}\t{g['end']}\t"
                    f"{g['class_code']}\t{len(g['transcripts'])}\t{g['exon_count']}\t"
                    f"{g['xloc']}\t{g['cmp_ref']}\n"
                )
    except Exception as e:
        print(f"ERROR writing {output_file}: {e}")
        sys.exit(1)

    return len(gene_data)

# Process novel_candidates.gtf
print("Creating table for novel_candidates.gtf...")
n = gtf_to_gene_table(
    '${OUTDIR}/unified_reference/novel_candidates.gtf',
    '${OUTDIR}/unified_reference/novel_candidates_table.tsv'
)
print(f"  Wrote {n} gene rows to novel_candidates_table.tsv")

# Process novel_no_same_strand_overlap.gtf
print("Creating table for novel_no_same_strand_overlap.gtf...")
n = gtf_to_gene_table(
    '${OUTDIR}/unified_reference/novel_no_same_strand_overlap.gtf',
    '${OUTDIR}/unified_reference/novel_no_same_strand_overlap_table.tsv'
)
print(f"  Wrote {n} gene rows to novel_no_same_strand_overlap_table.tsv")
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create gene-level TSV tables"
    exit 1
fi

echo ""
echo "Unified GTF creation completed successfully!"
echo ""
echo "Output files:"
echo ""
echo "In ${OUTDIR}/unified_reference/:"
echo "  - novel_candidates.gtf (filtered by class codes)"
echo "  - novel_candidates_stats.txt"
echo "  - novel_candidates_table.tsv (gene-level TSV table)"
echo "  - ref_protein_exons.gtf (protein-coding exons from reference)"
echo "  - novel_no_same_strand_overlap.gtf (after same-strand overlap filtering)"
echo "  - novel_no_same_strand_overlap_stats.txt"
echo "  - novel_no_same_strand_overlap_table.tsv (gene-level TSV table)"
echo ""
echo "In ${OUTDIR}/unified_reference/final_reference/:"
echo "  - merged_nofilter.gtf (reference + all novel candidates)"
echo "  - merged_nofilter_stats.txt"
echo "  - merged_filter.gtf (reference + filtered novel transcripts)"
echo "  - merged_filter_stats.txt"
echo ""

printf "\n\nCompleted Job: 12-create_unified_gtf\n\n"
pwd
date
