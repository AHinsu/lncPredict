#!/bin/bash
#SBATCH --job-name=13-annotation
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

###############################################################################
# 13_annotation.sh - NONCODE annotation and comprehensive statistics
# 
# This script annotates filtered lncRNAs (after unified GTF creation in script 12)
# against the NONCODE database using BLAST. It generates comprehensive statistics
# about matched and novel (non-matched) lncRNAs, both overall and per class-code.
#
# NOTE: Run this script AFTER 12_create_unified_gtf.sh. Script 12 performs the
#       necessary filtering steps (class codes, strand, overlap) before annotation.
#
# Prerequisites:
#   - Run 11_orf_pfam.sh first to generate final lncRNA files
#   - Run 12_create_unified_gtf.sh to generate the filtered candidate set
#   - Run 13a_setup_blast.sh to download NONCODE database (on node with internet)
#
# Usage: sbatch 13_annotation.sh
#
# Outputs:
#   ${OUTDIR}/unified_reference/annotation/filtered_lncrna.fasta
#   ${OUTDIR}/unified_reference/annotation/noncode_blast_results.txt
#   ${OUTDIR}/unified_reference/annotation/noncode_annotation_stats.txt
#   ${OUTDIR}/unified_reference/annotation/matched_lncrnas.txt
#   ${OUTDIR}/unified_reference/annotation/novel_lncrnas.txt
#   ${OUTDIR}/unified_reference/annotation/stats_by_class_code.txt
###############################################################################

printf "\n\nStarted Job: 13-annotation\n\n"
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

# Number of threads to use (will use SLURM_CPUS_PER_TASK if available)
THREADS=${SLURM_CPUS_PER_TASK:-12}


# Path to NONCODE database for BLAST (without extension)
# Set after running 13a_setup_blast.sh
LNCRNA_DATABASE="/storage/users/ahinsu/lncPredict/databases/lncrna_blast/lncrna_database"

# BLAST E-value threshold
BLAST_EVALUE=1e-5

# Maximum BLAST hits per query
BLAST_MAX_TARGETS=5

# Minimum alignment length for BLAST hits (nucleotides)
MIN_ALIGNMENT_LENGTH=200

# Minimum percent identity for a BLAST hit to be considered a match
IDENTITY_THRESHOLD=90

# Minimum query coverage (%) for a BLAST hit to be considered a match
MIN_QUERY_COVERAGE=80

# Minimum subject coverage (%) for a BLAST hit to be considered a match
MIN_SUBJECT_COVERAGE=70

# Filtered candidates from script 12 to be taken as input
FILTERED_GTF="${OUTDIR}/unified_reference/novel_no_same_strand_overlap.gtf" # Change this to novel_candidates.gtf if same strand overlapped file is not to be used further.


# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

# Create output directories
mkdir -p "${OUTDIR}/unified_reference/annotation"
mkdir -p logs

# Check input files from script 11
FULL_FASTA="${OUTDIR}/final_lncRNA/final_lncrna.fasta"
INPUT_INFO="${OUTDIR}/final_lncRNA/final_lncrna_info.tsv"

if [[ ! -f "${FULL_FASTA}" ]]; then
    echo "ERROR: Final lncRNA FASTA not found: ${FULL_FASTA}"
    echo "Please run 11_orf_pfam.sh first"
    exit 1
fi

if [[ ! -f "${INPUT_INFO}" ]]; then
    echo "ERROR: Final lncRNA info table not found: ${INPUT_INFO}"
    echo "Please run 11_orf_pfam.sh first"
    exit 1
fi

# Check filtered candidates from script 12
FILTERED_GTF="${OUTDIR}/unified_reference/novel_no_same_strand_overlap.gtf"

if [[ ! -f "${FILTERED_GTF}" ]]; then
    echo "ERROR: Filtered candidates GTF not found: ${FILTERED_GTF}"
    echo "Please run 12_create_unified_gtf.sh before running this script"
    exit 1
fi

# Extract transcript IDs from the filtered GTF and create a filtered FASTA
# for BLAST annotation (only annotate lncRNAs that passed all filtering steps)
INPUT_FASTA="${OUTDIR}/unified_reference/annotation/filtered_lncrna.fasta"

echo ""
echo "Extracting transcript IDs from filtered candidates GTF..."
python3 << PYEOF
import re
import sys

# Get transcript IDs from filtered GTF
filtered_ids = set()
transcript_id_pattern = re.compile(r'transcript_id "([^"]+)"')
try:
    with open('${FILTERED_GTF}', 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'transcript':
                continue
            for attr in parts[8].split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    m = transcript_id_pattern.search(attr)
                    if m:
                        filtered_ids.add(m.group(1))
                        break
except Exception as e:
    print(f"ERROR reading filtered GTF: {e}")
    sys.exit(1)

print(f"  Found {len(filtered_ids)} transcript IDs in filtered GTF")

# Filter the full FASTA to keep only filtered transcripts
kept = 0
skipped = 0
write_seq = False
try:
    with open('${FULL_FASTA}', 'r') as fin, \
         open('${INPUT_FASTA}', 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                write_seq = seq_id in filtered_ids
                if write_seq:
                    kept += 1
                    fout.write(line)
                else:
                    skipped += 1
            elif write_seq:
                fout.write(line)
except Exception as e:
    print(f"ERROR filtering FASTA: {e}")
    sys.exit(1)

print(f"  Kept {kept} sequences, skipped {skipped} sequences")
PYEOF

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create filtered FASTA for annotation"
    exit 1
fi

# Check if NONCODE database is available
if [[ -z "${LNCRNA_DATABASE:-}" ]]; then
    echo "ERROR: LNCRNA_DATABASE variable is not set"
    echo "Please set LNCRNA_DATABASE to the path of the NONCODE BLAST database (without extension)"
    echo "Or run 13a_setup_blast.sh first to download the database"
    exit 1
fi

# Check for NONCODE FASTA file
BLAST_FASTA="${LNCRNA_DATABASE}.fasta"
if [[ ! -f "${BLAST_FASTA}" ]]; then
    echo "ERROR: NONCODE BLAST database FASTA not found: ${BLAST_FASTA}"
    echo "Please run 13a_setup_blast.sh first to download the database"
    exit 1
fi

# Check if database is indexed, if not, index it
BLAST_INDEX="${LNCRNA_DATABASE}.nhr"

if [[ ! -f "${BLAST_INDEX}" ]]; then
    echo ""
    echo "BLAST database is not indexed. Indexing now..."
    echo "This may take several minutes depending on database size..."
    
    makeblastdb \
        -in "${BLAST_FASTA}" \
        -dbtype nucl \
        -out "${LNCRNA_DATABASE}" \
        -parse_seqids \
        -title "NONCODE_database" || {
            echo "ERROR: Failed to create BLAST database index"
            echo "Make sure BLAST+ tools are installed in your environment"
            exit 1
        }
    
    echo "Indexing completed successfully"
fi

# Check if input file is empty
if [[ ! -s "${INPUT_FASTA}" ]]; then
    echo "WARNING: Input FASTA is empty. No lncRNAs to annotate."
    echo "Creating placeholder output files..."
    
    echo "No filtered lncRNAs from script 12" > "${OUTDIR}/unified_reference/annotation/noncode_annotation_stats.txt"
    
    echo "Job completed (no lncRNAs to annotate)"
    exit 0
fi

# ============================================================================
# BLAST Annotation against NONCODE
# ============================================================================
echo ""
echo "Running BLAST against NONCODE database..."
echo "  Database: ${LNCRNA_DATABASE}"
echo "  E-value threshold: ${BLAST_EVALUE}"
echo "  Max target sequences: ${BLAST_MAX_TARGETS}"
echo "  Min alignment length: ${MIN_ALIGNMENT_LENGTH} nt"
echo "  Identity threshold: ${IDENTITY_THRESHOLD}%"
echo "  Min query coverage: ${MIN_QUERY_COVERAGE}%"
echo "  Min subject coverage: ${MIN_SUBJECT_COVERAGE}%"
echo "  Threads: ${THREADS}"

blastn \
    -query "${INPUT_FASTA}" \
    -db "${LNCRNA_DATABASE}" \
    -out "${OUTDIR}/unified_reference/annotation/noncode_blast_results_temp.txt" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue "${BLAST_EVALUE}" \
    -max_target_seqs "${BLAST_MAX_TARGETS}" \
    -num_threads ${THREADS} || {
        echo "ERROR: BLAST search failed"
        rm -f "${OUTDIR}/unified_reference/annotation/noncode_blast_results_temp.txt"
        exit 1
    }

# Add column headers to BLAST output file
echo "  Adding column headers to BLAST output..."
{
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen"
    cat "${OUTDIR}/unified_reference/annotation/noncode_blast_results_temp.txt"
} > "${OUTDIR}/unified_reference/annotation/noncode_blast_results.txt" || {
    echo "ERROR: Failed to add headers to BLAST output"
    rm -f "${OUTDIR}/unified_reference/annotation/noncode_blast_results_temp.txt"
    exit 1
}
rm -f "${OUTDIR}/unified_reference/annotation/noncode_blast_results_temp.txt"
    
# Analyze BLAST results with comprehensive statistics
echo "  Analyzing BLAST results and generating statistics..."
python3 << EOF
import re
import pandas as pd
from collections import defaultdict

# Read final lncRNA info to get class codes
lncrna_info = pd.read_csv('${INPUT_INFO}', sep='\t')

# Filter lncrna_info to only transcripts present in the filtered FASTA
# (those that passed all script 12 filtering steps)
filtered_ids = set()
try:
    with open('${INPUT_FASTA}', 'r') as f:
        for line in f:
            if line.startswith('>'):
                filtered_ids.add(line[1:].split()[0])
except Exception as e:
    print(f"Warning: Could not read filtered FASTA for ID list: {e}")

if filtered_ids:
    lncrna_info = lncrna_info[lncrna_info['transcript_id'].isin(filtered_ids)]
    print(f"Restricted annotation to {len(lncrna_info)} filtered lncRNAs (from script 12)")

# Handle class_code column (may be class_code, class_code_x, or class_code_y due to merge)
# In script 10, pandas creates class_code_x and class_code_y when merging dataframes.
# Script 11 should fix this, but for backward compatibility we handle all cases here.
# Preference: class_code > class_code_x > class_code_y
if 'class_code' in lncrna_info.columns:
    transcript_to_class = dict(zip(lncrna_info['transcript_id'], lncrna_info['class_code']))
elif 'class_code_x' in lncrna_info.columns:
    # Use class_code_x if class_code doesn't exist (from old format before script 11 fix)
    print("Warning: Found 'class_code_x' instead of 'class_code'. Using class_code_x.")
    transcript_to_class = dict(zip(lncrna_info['transcript_id'], lncrna_info['class_code_x']))
elif 'class_code_y' in lncrna_info.columns:
    # Use class_code_y if neither class_code nor class_code_x exist
    print("Warning: Found 'class_code_y' instead of 'class_code'. Using class_code_y.")
    transcript_to_class = dict(zip(lncrna_info['transcript_id'], lncrna_info['class_code_y']))
else:
    # If no class_code column exists at all, use 'unknown' as default
    print("Warning: No 'class_code' column found in input file. Using 'unknown' as default.")
    transcript_to_class = {tid: 'unknown' for tid in lncrna_info['transcript_id']}

print(f"Total filtered lncRNAs to annotate: {len(lncrna_info)}")

# Read BLAST results
hits = []
queries_with_hits = set()
best_hits = {}  # Best hit per query

min_length = ${MIN_ALIGNMENT_LENGTH}
identity_threshold = ${IDENTITY_THRESHOLD}
min_query_cov = ${MIN_QUERY_COVERAGE}
min_subject_cov = ${MIN_SUBJECT_COVERAGE}

try:
    with open('${OUTDIR}/unified_reference/annotation/noncode_blast_results.txt', 'r') as f:
        # Skip the header line
        header = f.readline()
        for line in f:
            if line.strip():
                fields = line.strip().split('\t')
                if len(fields) >= 14:
                    qseqid = fields[0]
                    sseqid = fields[1]
                    pident = float(fields[2])
                    length = int(fields[3])
                    evalue = fields[10]
                    bitscore = float(fields[11])
                    qstart = int(fields[6])
                    qend = int(fields[7])
                    sstart = int(fields[8])
                    send = int(fields[9])
                    qlen = int(fields[12])
                    slen = int(fields[13])
                    query_cov = abs(qend - qstart + 1) / qlen * 100
                    subject_cov = abs(send - sstart + 1) / slen * 100
                    
                    # Apply all filters: alignment length, identity, query coverage, subject coverage
                    if (length >= min_length and pident > identity_threshold
                            and query_cov > min_query_cov and subject_cov > min_subject_cov):
                        hits.append(fields)
                        queries_with_hits.add(qseqid)
                        
                        # Track best hit per query (by bitscore)
                        if qseqid not in best_hits or bitscore > best_hits[qseqid]['bitscore']:
                            best_hits[qseqid] = {
                                'sseqid': sseqid,
                                'pident': pident,
                                'length': length,
                                'evalue': evalue,
                                'bitscore': bitscore
                            }
except FileNotFoundError:
    print("No BLAST results file found")

total_lncrnas = len(lncrna_info)
matched_lncrnas = len(queries_with_hits)
novel_lncrnas = total_lncrnas - matched_lncrnas

print(f"Matched lncRNAs (with NONCODE hits): {matched_lncrnas}")
print(f"Novel lncRNAs (no NONCODE hits): {novel_lncrnas}")

# Collect identities of matched lncRNAs
if best_hits:
    identities = [hit['pident'] for hit in best_hits.values()]
else:
    identities = []

# Statistics by class code
class_code_stats = defaultdict(lambda: {
    'total': 0,
    'matched': 0,
    'novel': 0
})

for transcript_id in lncrna_info['transcript_id']:
    class_code = transcript_to_class.get(transcript_id, 'unknown')
    class_code_stats[class_code]['total'] += 1
    
    if transcript_id in queries_with_hits:
        class_code_stats[class_code]['matched'] += 1
    else:
        class_code_stats[class_code]['novel'] += 1

# Write matched lncRNAs list
with open('${OUTDIR}/unified_reference/annotation/matched_lncrnas.txt', 'w') as out:
    out.write("transcript_id\tclass_code\tbest_hit\tpercent_identity\talignment_length\tevalue\n")
    for transcript_id in sorted(queries_with_hits):
        class_code = transcript_to_class.get(transcript_id, 'unknown')
        if transcript_id in best_hits:
            hit = best_hits[transcript_id]
            out.write(f"{transcript_id}\t{class_code}\t{hit['sseqid']}\t{hit['pident']:.2f}\t{hit['length']}\t{hit['evalue']}\n")

# Write novel lncRNAs list
novel_transcript_ids = set(lncrna_info['transcript_id']) - queries_with_hits
with open('${OUTDIR}/unified_reference/annotation/novel_lncrnas.txt', 'w') as out:
    out.write("transcript_id\tclass_code\n")
    for transcript_id in sorted(novel_transcript_ids):
        class_code = transcript_to_class.get(transcript_id, 'unknown')
        out.write(f"{transcript_id}\t{class_code}\n")

# Write overall statistics
with open('${OUTDIR}/unified_reference/annotation/noncode_annotation_stats.txt', 'w') as out:
    out.write("NONCODE Annotation Statistics\n")
    out.write("=" * 70 + "\n\n")
    out.write(f"Database: ${LNCRNA_DATABASE}\n")
    out.write(f"E-value threshold: ${BLAST_EVALUE}\n")
    out.write(f"Min alignment length: {min_length} nt\n")
    out.write(f"Identity threshold: >{identity_threshold}%\n")
    out.write(f"Min query coverage: >{min_query_cov}%\n")
    out.write(f"Min subject coverage: >{min_subject_cov}%\n\n")
    
    out.write("Overall Statistics:\n")
    out.write("-" * 70 + "\n")
    out.write(f"Total final lncRNAs: {total_lncrnas}\n")
    out.write(f"Matched lncRNAs (with NONCODE hits): {matched_lncrnas} ({matched_lncrnas/total_lncrnas*100:.2f}%)\n")
    out.write(f"Novel lncRNAs (no NONCODE hits): {novel_lncrnas} ({novel_lncrnas/total_lncrnas*100:.2f}%)\n\n")
    
    if matched_lncrnas > 0:
        out.write("Identity Statistics (matched lncRNAs):\n")
        out.write("-" * 70 + "\n")
        avg_identity = sum(identities) / len(identities)
        out.write(f"Average identity: {avg_identity:.2f}%\n\n")
    
    out.write("Interpretation:\n")
    out.write("-" * 70 + "\n")
    out.write("- Matched lncRNAs have homologs in NONCODE (known lncRNAs)\n")
    out.write("- Novel lncRNAs have no significant NONCODE hits (potentially novel)\n")
    out.write("- All matched hits have >{}% identity, >{}% query coverage, and >{}% subject coverage\n".format(identity_threshold, min_query_cov, min_subject_cov))

# Write statistics by class code
with open('${OUTDIR}/unified_reference/annotation/stats_by_class_code.txt', 'w') as out:
    out.write("NONCODE Annotation Statistics by Class Code\n")
    out.write("=" * 70 + "\n\n")
    
    for class_code in sorted(class_code_stats.keys()):
        stats = class_code_stats[class_code]
        total = stats['total']
        matched = stats['matched']
        novel = stats['novel']
        
        out.write(f"Class Code: {class_code}\n")
        out.write("-" * 70 + "\n")
        out.write(f"Total lncRNAs: {total}\n")
        out.write(f"Matched (with NONCODE hits): {matched} ({matched/total*100:.2f}%)\n")
        out.write(f"Novel (no NONCODE hits): {novel} ({novel/total*100:.2f}%)\n")
        
        out.write("\n")

print("\nAnnotation analysis completed successfully")
EOF

echo ""
echo "Annotation files created in ${OUTDIR}/unified_reference/annotation/:"
echo "  - noncode_blast_results.txt (raw BLAST output)"
echo "  - matched_lncrnas.txt (lncRNAs with NONCODE hits)"
echo "  - novel_lncrnas.txt (lncRNAs without NONCODE hits)"
echo "  - noncode_annotation_stats.txt (overall statistics)"
echo "  - stats_by_class_code.txt (statistics per class code)"
echo ""

printf "\n\nCompleted Job: 13-annotation\n\n"
pwd
date
