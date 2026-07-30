#!/bin/bash
#SBATCH --job-name=11-orf_pfam
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20


###############################################################################
# 11_orf_pfam.sh - Pfam domain scan and filtering of protein-coding transcripts
# 
# This script uses TransDecoder-predicted ORFs from script 10 to scan for protein
# domains using Pfam. Transcripts with Pfam domains are identified as likely
# protein-coding and filtered out. The remaining transcripts are moved to the
# final_lncRNA folder.
#
# Prerequisites:
#   - Run 10_coding_consensus.sh first to generate candidate transcripts and ORFs
#   - Run 11a_setup_pfam.sh to download Pfam database (on node with internet)
#
# Usage: sbatch 11_orf_pfam.sh
#
# Outputs:
#   ${OUTDIR}/lncrna_candidates/pfam_analysis/pfamscan_output.txt
#   ${OUTDIR}/lncrna_candidates/pfam_analysis/pfam_hits_by_transcript.txt
#   ${OUTDIR}/lncrna_candidates/pfam_analysis/pfam_stats.txt
#   ${OUTDIR}/lncrna_candidates/pfam_analysis/transcripts_with_domains.txt
#   ${OUTDIR}/final_lncRNA/final_lncrna.gtf
#   ${OUTDIR}/final_lncRNA/final_lncrna.fasta
#   ${OUTDIR}/final_lncRNA/final_lncrna_ids.txt
#   ${OUTDIR}/final_lncRNA/final_lncrna_info.tsv
#   ${OUTDIR}/final_lncRNA/filtering_stats.txt
###############################################################################

printf "\n\nStarted Job: 11-orf_pfam\n\n"
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

# Path to Pfam database (Pfam-A.hmm) - set after running 11a_setup_pfam.sh
PFAM_DATABASE="/storage/users/ahinsu/lncPredict/databases/pfam/Pfam-A.hmm"

# Pfam E-value threshold
PFAM_EVALUE=0.001

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
mkdir -p "${OUTDIR}/lncrna_candidates/pfam_analysis"
mkdir -p "${OUTDIR}/final_lncRNA"
mkdir -p logs

# Check input files from script 10
INPUT_FASTA="${OUTDIR}/lncrna_candidates/candidate_transcripts.fasta"
INPUT_GTF="${OUTDIR}/lncrna_candidates/candidate_transcripts.gtf"
INPUT_INFO="${OUTDIR}/lncrna_candidates/candidate_info.tsv"
TRANSDECODER_PEP="${OUTDIR}/lncrna_candidates/transdecoder/longest_orfs.pep"

if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "ERROR: Candidate FASTA not found: ${INPUT_FASTA}"
    echo "Please run 10_coding_consensus.sh first"
    exit 1
fi

if [[ ! -f "${INPUT_GTF}" ]]; then
    echo "ERROR: Candidate GTF not found: ${INPUT_GTF}"
    echo "Please run 10_coding_consensus.sh first"
    exit 1
fi

if [[ ! -f "${TRANSDECODER_PEP}" ]]; then
    echo "ERROR: TransDecoder ORF predictions not found: ${TRANSDECODER_PEP}"
    echo "Please run 10_coding_consensus.sh first"
    exit 1
fi

# Check if Pfam database is available
if [[ -z "${PFAM_DATABASE:-}" ]]; then
    echo "ERROR: PFAM_DATABASE variable is not set"
    echo "Please set PFAM_DATABASE to the path of Pfam-A.hmm"
    echo "Or run 11a_setup_pfam.sh first to download the database"
    exit 1
fi

if [[ ! -f "${PFAM_DATABASE}" ]]; then
    echo "ERROR: Pfam database not found: ${PFAM_DATABASE}"
    echo "Please run 11a_setup_pfam.sh first to download the database"
    exit 1
fi

# Check if database is indexed, if not, index it
PFAM_DIR=$(dirname "${PFAM_DATABASE}")
PFAM_BASE=$(basename "${PFAM_DATABASE}")
PFAM_INDEX="${PFAM_DIR}/${PFAM_BASE}.h3m"

if [[ ! -f "${PFAM_INDEX}" ]]; then
    echo ""
    echo "Pfam database is not indexed. Indexing now..."
    echo "This may take 10-20 minutes..."
    
    cd "${PFAM_DIR}" || exit 1
    hmmpress "${PFAM_BASE}" || {
        echo "ERROR: Failed to index Pfam database"
        echo "Make sure HMMER3 (hmmpress) is installed in your environment"
        exit 1
    }
    cd - > /dev/null || exit 1
    
    echo "Indexing completed successfully"
fi

# Check if input file is empty
if [[ ! -s "${INPUT_FASTA}" ]]; then
    echo "WARNING: Input FASTA is empty. No candidates to analyze."
    echo "Creating placeholder output files..."
    
    echo "No candidates selected in script 10" > "${OUTDIR}/lncrna_candidates/pfam_analysis/pfam_stats.txt"
    
    # Create empty final files
    touch "${OUTDIR}/final_lncRNA/final_lncrna.gtf"
    touch "${OUTDIR}/final_lncRNA/final_lncrna.fasta"
    touch "${OUTDIR}/final_lncRNA/final_lncrna_ids.txt"
    touch "${OUTDIR}/final_lncRNA/final_lncrna_info.tsv"
    
    echo "Job completed (no candidates to analyze)"
    exit 0
fi

# ============================================================================
# Pfam-scan Analysis using TransDecoder ORF predictions
# ============================================================================
echo ""
echo "Running Pfam-scan on TransDecoder ORF predictions..."
echo "  Database: ${PFAM_DATABASE}"
echo "  E-value threshold: ${PFAM_EVALUE}"
echo "  Threads: ${THREADS}"
echo "  Input ORFs: ${TRANSDECODER_PEP}"

PFAM_DIR=$(dirname "${PFAM_DATABASE}")

if [[ -f "${TRANSDECODER_PEP}" && -s "${TRANSDECODER_PEP}" ]]; then
    pfam_scan.pl \
        -fasta "${TRANSDECODER_PEP}" \
        -dir "${PFAM_DIR}" \
        -e_seq "${PFAM_EVALUE}" \
        -cpu ${THREADS} \
        -outfile "${OUTDIR}/lncrna_candidates/pfam_analysis/pfamscan_output.txt" 2>&1 || {
            echo "  pfam_scan.pl failed (this is common if no domains are found)"
        }
        
        # Analyze Pfam results
        echo "  Analyzing Pfam results..."
        python3 << EOF
import re
from collections import defaultdict

# Parse Pfam-scan output
pfam_hits = defaultdict(list)
total_hits = 0

try:
    with open('${OUTDIR}/lncrna_candidates/pfam_analysis/pfamscan_output.txt', 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            # Parse Pfam output
            # Format: <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>
            fields = line.strip().split()
            if len(fields) >= 13:
                seq_id = fields[0]
                pfam_acc = fields[5]
                pfam_name = fields[6]
                evalue = fields[12]
                
                # Extract transcript ID from TransDecoder ORF ID
                # TransDecoder format: TRANSCRIPT_ID.p# (e.g., MSTRG.1.1.p1)
                # Extract everything before .p
                if '.p' in seq_id:
                    transcript_id = seq_id.rsplit('.p', 1)[0]
                else:
                    transcript_id = seq_id
                
                pfam_hits[transcript_id].append({
                    'pfam_acc': pfam_acc,
                    'pfam_name': pfam_name,
                    'evalue': evalue,
                    'orf_id': seq_id
                })
                total_hits += 1
except FileNotFoundError:
    print("No Pfam output file found")
except Exception as e:
    print(f"Error parsing Pfam output: {e}")

transcripts_with_domains = len(pfam_hits)

# Write detailed Pfam hits
if pfam_hits:
    with open('${OUTDIR}/lncrna_candidates/pfam_analysis/pfam_hits_by_transcript.txt', 'w') as out:
        out.write("transcript_id\tpfam_acc\tpfam_name\tevalue\torf_id\n")
        for transcript_id in sorted(pfam_hits.keys()):
            for hit in pfam_hits[transcript_id]:
                out.write(f"{transcript_id}\t{hit['pfam_acc']}\t{hit['pfam_name']}\t{hit['evalue']}\t{hit['orf_id']}\n")
    
    # Write list of transcripts with Pfam domains (to be filtered out)
    with open('${OUTDIR}/lncrna_candidates/pfam_analysis/transcripts_with_domains.txt', 'w') as out:
        for transcript_id in sorted(pfam_hits.keys()):
            domain_count = len(pfam_hits[transcript_id])
            out.write(f"{transcript_id}\t{domain_count}\n")

# Write statistics
with open('${OUTDIR}/lncrna_candidates/pfam_analysis/pfam_stats.txt', 'w') as out:
    out.write("Pfam-scan Analysis\n")
    out.write("=" * 70 + "\n\n")
    out.write(f"Database: ${PFAM_DATABASE}\n")
    out.write(f"E-value threshold: ${PFAM_EVALUE}\n\n")
    
    out.write(f"Total Pfam domain hits: {total_hits}\n")
    out.write(f"Transcripts with Pfam domains: {transcripts_with_domains}\n\n")
    
    out.write("Note: Transcripts with Pfam domains are likely protein-coding\n")
    out.write("and will be filtered out from the final lncRNA set.\n")

print(f"Pfam-scan completed:")
print(f"  Domain hits: {total_hits}")
print(f"  Transcripts with domains: {transcripts_with_domains}")
EOF
else
    echo "  No ORF predictions found or file is empty"
    echo "No ORFs to scan" > "${OUTDIR}/lncrna_candidates/pfam_analysis/pfam_stats.txt"
fi

# ============================================================================
# Filter transcripts and create final lncRNA files
# ============================================================================
echo ""
echo "Filtering transcripts with Pfam domains and creating final lncRNA set..."

python3 << EOF
import pandas as pd

# Read candidate transcript IDs
candidate_ids = set()
with open('${OUTDIR}/lncrna_candidates/candidate_transcripts.txt', 'r') as f:
    for line in f:
        tid = line.strip()
        if tid:
            candidate_ids.add(tid)

print(f"Total candidate transcripts: {len(candidate_ids)}")

# Read transcripts with Pfam domains (to be filtered out)
pfam_transcripts = set()
pfam_file = '${OUTDIR}/lncrna_candidates/pfam_analysis/transcripts_with_domains.txt'
try:
    with open(pfam_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields:
                pfam_transcripts.add(fields[0])
except FileNotFoundError:
    print("No Pfam domains file found - no transcripts to filter")

print(f"Transcripts with Pfam domains (protein-coding): {len(pfam_transcripts)}")

# Filter out transcripts with Pfam domains
final_lncrna_ids = candidate_ids - pfam_transcripts
print(f"Final lncRNA transcripts (after filtering): {len(final_lncrna_ids)}")

# Save final lncRNA IDs
with open('${OUTDIR}/final_lncRNA/final_lncrna_ids.txt', 'w') as out:
    for tid in sorted(final_lncrna_ids):
        out.write(f"{tid}\n")

# Read candidate info and filter
final_info = None
try:
    candidate_info = pd.read_csv('${INPUT_INFO}', sep='\t')
    final_info = candidate_info[candidate_info['transcript_id'].isin(final_lncrna_ids)]
    
    # Fix class_code column if it was duplicated during merge (class_code_x and class_code_y)
    # In script 10, both detailed_overlap and class_code_info contain class_code, causing
    # pandas to create class_code_x and class_code_y. Both columns should have identical
    # values, but we use class_code_x as it comes from the primary working dataset.
    if 'class_code_x' in final_info.columns and 'class_code_y' in final_info.columns:
        # Use class_code_x and drop class_code_y (both should be identical)
        final_info['class_code'] = final_info['class_code_x']
        final_info = final_info.drop(columns=['class_code_x', 'class_code_y'])
    elif 'class_code_x' in final_info.columns:
        final_info['class_code'] = final_info['class_code_x']
        final_info = final_info.drop(columns=['class_code_x'])
    elif 'class_code_y' in final_info.columns:
        final_info['class_code'] = final_info['class_code_y']
        final_info = final_info.drop(columns=['class_code_y'])
    
    final_info.to_csv('${OUTDIR}/final_lncRNA/final_lncrna_info.tsv', sep='\t', index=False)
    print(f"Saved final lncRNA info table: {len(final_info)} transcripts")
except Exception as e:
    print(f"Warning: Could not process candidate info: {e}")

# Note: Gene counting and statistics will be generated after GTF extraction
# where gene_id information is reliably available from the GTF file itself.
print("\nFiltering completed successfully")
print(f"Final lncRNA transcripts: {len(final_lncrna_ids)}")
EOF

# Extract final FASTA sequences
echo ""
echo "Extracting final lncRNA FASTA sequences..."
seqkit grep -f "${OUTDIR}/final_lncRNA/final_lncrna_ids.txt" \
    "${INPUT_FASTA}" > "${OUTDIR}/final_lncRNA/final_lncrna.fasta"

FINAL_FASTA_COUNT=$(grep -c "^>" "${OUTDIR}/final_lncRNA/final_lncrna.fasta" || echo "0")
echo "  Created: final_lncrna.fasta (${FINAL_FASTA_COUNT} sequences)"

# Extract final GTF
echo ""
echo "Extracting final lncRNA GTF..."
python3 << EOF
import re
from collections import defaultdict

# Read final lncRNA IDs
final_ids = set()
with open('${OUTDIR}/final_lncRNA/final_lncrna_ids.txt', 'r') as f:
    for line in f:
        tid = line.strip()
        if tid:
            final_ids.add(tid)

print(f"Extracting GTF for {len(final_ids)} final lncRNA transcripts...")

# Extract from candidate GTF and collect class_code information
extracted = 0
extracted_genes = set()
extracted_transcripts = set()
class_code_transcripts = defaultdict(set)
class_code_genes = defaultdict(set)

with open('${INPUT_GTF}', 'r') as infile, \
     open('${OUTDIR}/final_lncRNA/final_lncrna.gtf', 'w') as outfile:
    
    for line in infile:
        if line.startswith('#'):
            outfile.write(line)
            continue
        
        if '\ttranscript\t' in line or '\texon\t' in line:
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                feature = parts[2]
                attributes = parts[8]
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
                
                if transcript_id and transcript_id in final_ids:
                    outfile.write(line)
                    extracted += 1
                    extracted_transcripts.add(transcript_id)
                    if gene_id:
                        extracted_genes.add(gene_id)
                    # Collect class_code info from transcript lines only
                    if feature == "transcript" and class_code:
                        class_code_transcripts[class_code].add(transcript_id)
                        if gene_id:
                            class_code_genes[class_code].add(gene_id)

print(f"Extracted {extracted} GTF entries")
print(f"Unique genes: {len(extracted_genes)}")
print(f"Unique transcripts: {len(extracted_transcripts)}")

# Save gene and transcript counts for statistics generation
with open('${OUTDIR}/final_lncRNA/gtf_counts.txt', 'w') as f:
    f.write(f"genes:{len(extracted_genes)}\n")
    f.write(f"transcripts:{len(extracted_transcripts)}\n")
    # Save class_code statistics
    for class_code in sorted(class_code_transcripts.keys()):
        n_transcripts = len(class_code_transcripts[class_code])
        n_genes = len(class_code_genes[class_code])
        f.write(f"class_code_{class_code}:{n_transcripts},{n_genes}\n")
EOF

# ============================================================================
# Generate filtering statistics with accurate gene counts from GTF
# ============================================================================
echo ""
echo "Generating filtering statistics with gene counts from GTF..."

python3 << EOF
import pandas as pd

# Read counts from GTF extraction
gtf_genes = 0
gtf_transcripts = 0
class_code_stats = {}

try:
    with open('${OUTDIR}/final_lncRNA/gtf_counts.txt', 'r') as f:
        for line in f:
            if line.startswith('genes:'):
                gtf_genes = int(line.strip().split(':')[1])
            elif line.startswith('transcripts:'):
                gtf_transcripts = int(line.strip().split(':')[1])
            elif line.startswith('class_code_'):
                # Format: class_code_X:transcripts,genes
                parts = line.strip().split(':')
                class_code = parts[0].replace('class_code_', '')
                counts = parts[1].split(',')
                class_code_stats[class_code] = {
                    'transcripts': int(counts[0]),
                    'genes': int(counts[1])
                }
except FileNotFoundError:
    print("Warning: GTF counts file not found")

# Read candidate IDs for initial count
candidate_ids = set()
with open('${OUTDIR}/lncrna_candidates/candidate_transcripts.txt', 'r') as f:
    for line in f:
        tid = line.strip()
        if tid:
            candidate_ids.add(tid)

# Read Pfam-filtered transcripts
pfam_transcripts = set()
pfam_file = '${OUTDIR}/lncrna_candidates/pfam_analysis/transcripts_with_domains.txt'
try:
    with open(pfam_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields:
                pfam_transcripts.add(fields[0])
except FileNotFoundError:
    pass

# Generate filtering statistics
with open('${OUTDIR}/final_lncRNA/filtering_stats.txt', 'w') as out:
    out.write("lncRNA Filtering Statistics (Script 11)\n")
    out.write("=" * 70 + "\n\n")
    out.write(f"Initial candidate transcripts: {len(candidate_ids)}\n")
    out.write(f"Transcripts with Pfam domains (removed): {len(pfam_transcripts)}\n")
    out.write(f"Final lncRNA transcripts: {gtf_transcripts}\n")
    if len(candidate_ids) > 0:
        out.write(f"Retention rate: {gtf_transcripts/len(candidate_ids)*100:.2f}%\n\n")
    else:
        out.write("Retention rate: N/A\n\n")
    
    # Add final numbers section
    out.write("Final lncRNA Set:\n")
    out.write("-" * 70 + "\n")
    out.write(f"Total transcripts: {gtf_transcripts}\n")
    out.write(f"Total genes: {gtf_genes}\n\n")
    
    # Add per-class-code statistics
    if class_code_stats:
        out.write("Breakdown by Class Code:\n")
        out.write("-" * 70 + "\n")
        for class_code in sorted(class_code_stats.keys()):
            stats = class_code_stats[class_code]
            out.write(f"  Class code '{class_code}':\n")
            out.write(f"    Transcripts: {stats['transcripts']}\n")
            out.write(f"    Genes: {stats['genes']}\n")
        out.write("\n")
    
    out.write("Note: Transcripts with Pfam protein domains were identified as\n")
    out.write("likely protein-coding and removed from the final lncRNA set.\n")

print(f"Statistics generated: {gtf_transcripts} transcripts, {gtf_genes} genes")
if class_code_stats:
    print(f"Class codes found: {', '.join(sorted(class_code_stats.keys()))}")
EOF

echo ""
echo "Final lncRNA files created in ${OUTDIR}/final_lncRNA/:"
echo "  - final_lncrna_ids.txt (list of transcript IDs)"
echo "  - final_lncrna.fasta (sequences)"
echo "  - final_lncrna.gtf (annotations)"
echo "  - final_lncrna_info.tsv (metadata table)"
echo "  - filtering_stats.txt (filtering statistics)"
echo ""
echo "Pfam analysis results in ${OUTDIR}/lncrna_candidates/pfam_analysis/:"
echo "  - pfamscan_output.txt (raw Pfam-scan output)"
echo "  - pfam_hits_by_transcript.txt (detailed hits)"
echo "  - transcripts_with_domains.txt (filtered transcripts)"
echo "  - pfam_stats.txt (summary statistics)"
echo ""

printf "\n\nCompleted Job: 11-orf_pfam\n\n"
pwd
date
