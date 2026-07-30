#!/bin/bash
#SBATCH --job-name=05a-cpat
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

###############################################################################
# 05a_cpat.sh - CPAT coding potential assessment with custom model building
# 
# This script builds CPAT models from reference FASTA and GTF files, then
# runs CPAT (Coding-Potential Assessment Tool) to identify non-coding
# transcripts based on sequence features.
#
# Use this script when:
#   - You need to build custom CPAT models for your organism
#   - Pre-built models are not available for your species
#   - You want to use organism-specific reference data from Ensembl
#
# Requirements:
#   - CPAT (installed via conda)
#   - gffread (for sequence extraction) - conda install -c bioconda gffread
#   - R with ROCR package (for automatic threshold calculation):
#       conda install -c conda-forge r-base r-rocr
#
# Usage: sbatch 05a_cpat.sh
#
# Outputs:
#   ${OUTDIR}/coding_potential/cpat/models/hexamer.tsv
#   ${OUTDIR}/coding_potential/cpat/models/logitModel.logit.RData
#   ${OUTDIR}/coding_potential/cpat/models/logitModel.feature.xls
#   ${OUTDIR}/coding_potential/cpat/models/logitModel.make_logitModel.r
#   ${OUTDIR}/coding_potential/cpat/models/cpat_threshold.txt (if auto)
#   ${OUTDIR}/coding_potential/cpat/models/cpat_performance_plots.png (if auto)
#   ${OUTDIR}/coding_potential/cpat/models/cpat_stats.txt (if auto)
#   ${OUTDIR}/coding_potential/cpat/cpat_output.txt
#   ${OUTDIR}/coding_potential/cpat/cpat_stats.txt
#   ${OUTDIR}/coding_potential/cpat/cpat_noncoding.txt
###############################################################################

printf "\n\nStarted Job: 05a-cpat\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results" # Make sure to provide the same path for all scripts

# Reference files for model building (from Ensembl or other source)
# These should contain known coding and non-coding sequences
REFERENCE_FASTA="/storage/users/ahinsu/lncPredict/reference/Bos_taurus.chromsomal.fa"
REFERENCE_GTF="/storage/users/ahinsu/lncPredict/reference/Bos_taurus.cleaned.gtf"

# CPAT threshold calculation method
# Options:
#   "auto" - Automatically calculate optimal threshold using ROC analysis
#   <number> - Use a specific threshold value (e.g., 0.364 for human)
# 
# Pre-defined thresholds for reference:
#   Human: 0.364, Mouse: 0.44, Fly: 0.39, Zebrafish: 0.38
# 
# For custom models, "auto" is recommended to determine the optimal threshold
# based on your specific training data
CPAT_THRESHOLD="auto"

# Training set sizes (optional, leave empty for defaults)
# Number of coding genes to use for training (default: 1000)
CODING_TRAINING_SIZE=1000

# Number of non-coding genes to use for training (default: 1000)
NONCODING_TRAINING_SIZE=1000

# Conda environment name (created from environment.yml)
CONDA_ENV="lncrna"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

# Create output directories
mkdir -p "${OUTDIR}/coding_potential/cpat/models"
mkdir -p "${OUTDIR}/coding_potential/cpat/training_data"
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

# Check input file
INPUT_FASTA="${OUTDIR}/filtered_transcripts/filtered.fasta"
if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "ERROR: Filtered FASTA not found: ${INPUT_FASTA}"
    echo "Please run 04_filter_transcripts.sh first"
    exit 1
fi

echo "============================================================================"
echo "Step 1: Extracting training sequences from reference"
echo "============================================================================"
echo ""

# Extract coding sequences from reference GTF
echo "Extracting coding sequences..."
python3 - "${REFERENCE_GTF}" "${OUTDIR}/coding_potential/cpat/training_data/coding_gene_ids.txt" ${CODING_TRAINING_SIZE} << 'EOF1'
import sys
import re

gtf_file = sys.argv[1]
output_file = sys.argv[2]
max_genes = int(sys.argv[3]) if len(sys.argv) > 3 else 1000

coding_genes = set()

# Read GTF and extract protein-coding genes
with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        feature = fields[2]
        attributes = fields[8]
        
        # Look for protein_coding genes
        if 'gene_biotype "protein_coding"' in attributes or \
           'transcript_biotype "protein_coding"' in attributes:
            # Extract gene_id
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            if gene_id_match:
                coding_genes.add(gene_id_match.group(1))
                if len(coding_genes) >= max_genes:
                    break

print(f"Found {len(coding_genes)} coding genes")

# Write gene IDs to file
with open(output_file, 'w') as out:
    for gene_id in sorted(coding_genes):
        out.write(f"{gene_id}\n")

EOF1

# Extract non-coding sequences from reference GTF
echo ""
echo "Extracting lncRNA sequences (lncRNA biotype)..."
python3 - "${REFERENCE_GTF}" "${OUTDIR}/coding_potential/cpat/training_data/noncoding_gene_ids.txt" ${NONCODING_TRAINING_SIZE} << 'EOF2'
import sys
import re

gtf_file = sys.argv[1]
output_file = sys.argv[2]
max_genes = int(sys.argv[3]) if len(sys.argv) > 3 else 1000

noncoding_genes = set()

# Read GTF and extract lncRNA genes
with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        feature = fields[2]
        attributes = fields[8]
        
        # Look for lncRNA biotype specifically
        if 'gene_biotype "lncRNA"' in attributes or \
           'transcript_biotype "lncRNA"' in attributes:
            # Extract gene_id
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            if gene_id_match:
                noncoding_genes.add(gene_id_match.group(1))
                if len(noncoding_genes) >= max_genes:
                    break

print(f"Found {len(noncoding_genes)} lncRNA genes")

# Write gene IDs to file
with open(output_file, 'w') as out:
    for gene_id in sorted(noncoding_genes):
        out.write(f"{gene_id}\n")

EOF2

echo ""
echo "============================================================================"
echo "Step 2: Extracting coding CDS and non-coding sequences"
echo "============================================================================"
echo ""

# Extract coding CDS sequences (without UTR) for hexamer table
echo "Extracting coding CDS sequences (without UTR)..."
CODING_CDS_FASTA="${OUTDIR}/coding_potential/cpat/training_data/coding_cds.fasta"
CODING_CDS_GTF="${OUTDIR}/coding_potential/cpat/training_data/coding_cds.gtf"

# Filter GTF to get only CDS features for protein-coding genes
python3 - "${REFERENCE_GTF}" "${OUTDIR}/coding_potential/cpat/training_data/coding_gene_ids.txt" "${CODING_CDS_GTF}" << 'EOF_EXTRACT_CDS'
import sys
import re

gtf_file = sys.argv[1]
gene_ids_file = sys.argv[2]
output_gtf = sys.argv[3]

# Read target gene IDs
with open(gene_ids_file, 'r') as f:
    target_genes = set(line.strip() for line in f)

print(f"Target coding genes: {len(target_genes)}")

# Filter GTF to only include CDS features of target protein-coding genes
count = 0
with open(gtf_file, 'r') as f_in, open(output_gtf, 'w') as f_out:
    for line in f_in:
        if line.startswith('#'):
            f_out.write(line)
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        feature = fields[2]
        attributes = fields[8]
        
        # Only include CDS features (not exons, not full transcript)
        if feature != 'CDS':
            continue
        
        # Extract gene_id
        gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
        
        if gene_id_match:
            gene_id = gene_id_match.group(1)
            
            if gene_id in target_genes:
                f_out.write(line)
                count += 1

print(f"Wrote {count} CDS entries for coding genes")

EOF_EXTRACT_CDS

# Extract CDS sequences using gffread
if command -v gffread &> /dev/null; then
    echo "Extracting CDS sequences with gffread..."
    # Use -x to extract only CDS (coding sequence without UTR)
    gffread -x "${CODING_CDS_FASTA}" -g "${REFERENCE_FASTA}" "${CODING_CDS_GTF}"
else
    echo "ERROR: gffread not found. Please install it to extract CDS sequences."
    echo "Install via: conda install -c bioconda gffread"
    exit 1
fi

# Extract non-coding sequences
echo ""
echo "Extracting non-coding sequences (full transcripts)..."
NONCODING_FASTA="${OUTDIR}/coding_potential/cpat/training_data/noncoding.fasta"

# Create a temporary GTF with only the target non-coding genes
NONCODING_GTF="${OUTDIR}/coding_potential/cpat/training_data/noncoding.gtf"

python3 - "${REFERENCE_GTF}" "${OUTDIR}/coding_potential/cpat/training_data/noncoding_gene_ids.txt" "${NONCODING_GTF}" << 'EOF_EXTRACT_NONCODING_GTF_FULL'
import sys
import re

gtf_file = sys.argv[1]
gene_ids_file = sys.argv[2]
output_gtf = sys.argv[3]

# Read target gene IDs
with open(gene_ids_file, 'r') as f:
    target_genes = set(line.strip() for line in f)

print(f"Target non-coding genes: {len(target_genes)}")

# Filter GTF to only include target genes
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
        
        # Extract gene_id
        gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
        
        if gene_id_match:
            gene_id = gene_id_match.group(1)
            
            if gene_id in target_genes:
                f_out.write(line)
                count += 1

print(f"Wrote {count} GTF entries for non-coding genes")

EOF_EXTRACT_NONCODING_GTF_FULL

# Use gffread to extract sequences
if command -v gffread &> /dev/null; then
    gffread -w "${NONCODING_FASTA}" -g "${REFERENCE_FASTA}" "${NONCODING_GTF}"
else
    echo "ERROR: gffread not found. Please install it to extract sequences."
    echo "Install via: conda install -c bioconda gffread"
    exit 1
fi

echo ""
echo "============================================================================"
echo "Step 3: Building hexamer frequency table"
echo "============================================================================"
echo ""

# Build hexamer table using make_hexamer_tab.py
# This creates a table of hexamer frequencies from the reference sequences
HEXAMER_TABLE="${OUTDIR}/coding_potential/cpat/models/hexamer.tsv"

echo "Building hexamer table using CPAT's make_hexamer_tab.py..."
echo "This may take several minutes..."

# Check if make_hexamer_tab.py is available
if ! command -v make_hexamer_tab.py &> /dev/null; then
    echo "ERROR: make_hexamer_tab.py not found. Please install CPAT."
    echo "Install via: conda install -c bioconda cpat"
    exit 1
fi

# Run make_hexamer_tab.py
echo "Command: make_hexamer_tab.py -c \"${CODING_CDS_FASTA}\" -n \"${NONCODING_FASTA}\" > ${HEXAMER_TABLE}"
if ! make_hexamer_tab.py -c "${CODING_CDS_FASTA}" -n "${NONCODING_FASTA}" > "${HEXAMER_TABLE}"; then
    echo "ERROR: make_hexamer_tab.py failed"
    exit 1
fi

echo "Hexamer table created: ${HEXAMER_TABLE}"
if [[ -f "${HEXAMER_TABLE}" ]]; then
    echo "Table size: $(wc -l < ${HEXAMER_TABLE}) lines"
else
    echo "ERROR: Hexamer table file not created"
    exit 1
fi

echo ""
echo "============================================================================"
echo "Step 4: Building logit model"
echo "============================================================================"
echo ""

# Build logit model using make_logitModel.py
LOGIT_MODEL_PREFIX="${OUTDIR}/coding_potential/cpat/models/logitModel"
LOGIT_MODEL="${LOGIT_MODEL_PREFIX}.logit.RData"

echo "Building logit regression model using CPAT's make_logitModel.py..."
echo "This may take several minutes..."

# Check if make_logitModel.py is available
if ! command -v make_logitModel.py &> /dev/null; then
    echo "ERROR: make_logitModel.py not found. Please install CPAT."
    echo "Install via: conda install -c bioconda cpat"
    exit 1
fi

# Run make_logitModel.py
# Note: make_logitModel.py requires the sequences in FASTA format
# It will generate:
#   - ${LOGIT_MODEL_PREFIX}.feature.xls (features from training data)
#   - ${LOGIT_MODEL_PREFIX}.logit.RData (R logit model)
#   - ${LOGIT_MODEL_PREFIX}.make_logitModel.r (R script used to build model)

echo "Command: make_logitModel.py -x \"${HEXAMER_TABLE}\" -c \"${CODING_CDS_FASTA}\" -n \"${NONCODING_FASTA}\" -o \"${LOGIT_MODEL_PREFIX}\""
if ! make_logitModel.py -x "${HEXAMER_TABLE}" -c "${CODING_CDS_FASTA}" -n "${NONCODING_FASTA}" -o "${LOGIT_MODEL_PREFIX}"; then
    echo "ERROR: make_logitModel.py failed"
    exit 1
fi

echo ""
echo "Logit model created:"
echo "  Model file: ${LOGIT_MODEL}"
echo "  Features file: ${LOGIT_MODEL_PREFIX}.feature.xls"
echo "  R script: ${LOGIT_MODEL_PREFIX}.make_logitModel.r"

# Check if model was created
if [[ ! -f "${LOGIT_MODEL}" ]]; then
    echo "WARNING: Logit model file not found. Check if R is installed and working."
    echo "The model file should be at: ${LOGIT_MODEL}"
    echo "You may need to run the R script manually:"
    echo "  Rscript ${LOGIT_MODEL_PREFIX}.make_logitModel.r"
fi

echo ""
echo "============================================================================"
echo "Step 5: Calculate optimal threshold using 10-fold cross-validation"
echo "============================================================================"
echo ""

# Determine threshold to use
FINAL_THRESHOLD="${CPAT_THRESHOLD}"

if [[ "${CPAT_THRESHOLD}" == "auto" ]]; then
    echo "Automatic threshold calculation requested"
    echo ""
    
    # Use the feature file created by make_logitModel.py
    # This file contains: ID, mRNA, ORF, Fickett, Hexamer, Label
    # Label: 1 = coding, 0 = noncoding
    FEATURE_FILE="${LOGIT_MODEL_PREFIX}.feature.xls"
    
    if [[ ! -f "${FEATURE_FILE}" ]]; then
        echo "ERROR: Feature file not found: ${FEATURE_FILE}"
        echo "The feature file should have been created by make_logitModel.py"
        exit 1
    fi
    
    echo "Using feature file for 10-fold cross-validation:"
    echo "  Feature file: ${FEATURE_FILE}"
    echo ""
    echo "This file contains extracted features (mRNA, ORF, Fickett, Hexamer) and labels"
    echo "from both coding and non-coding training sequences."
    echo ""
    
    # Count sequences in feature file
    NUM_SEQUENCES=$(tail -n +2 "${FEATURE_FILE}" | wc -l)
    echo "  Total training sequences: ${NUM_SEQUENCES}"
    echo ""
    
    echo "Performing 10-fold cross-validation to calculate optimal threshold..."
    echo "  This matches the CPAT original publication methodology:"
    echo "    - Split data into 10 folds"
    echo "    - For each fold: train logistic regression on 90%, predict on 10%"
    echo "    - Combine all predictions for ROC analysis"
    echo "    - Find threshold where sensitivity = specificity"
    echo ""
    
    # Run R script to calculate threshold using 10-fold cross-validation
    THRESHOLD_PREFIX="${OUTDIR}/coding_potential/cpat/models/cpat"
    Rscript scripts/calculate_cpat_threshold.R "${FEATURE_FILE}" "${THRESHOLD_PREFIX}"
    
    # Read calculated threshold
    THRESHOLD_FILE="${THRESHOLD_PREFIX}_threshold.txt"
    if [[ -f "${THRESHOLD_FILE}" ]]; then
        FINAL_THRESHOLD=$(cat "${THRESHOLD_FILE}")
        echo ""
        echo "Calculated optimal threshold: ${FINAL_THRESHOLD}"
        echo "See ${THRESHOLD_PREFIX}_stats.txt for detailed statistics"
        echo "See ${THRESHOLD_PREFIX}_performance_plots.png for performance visualizations"
    else
        echo "ERROR: Threshold calculation failed"
        echo "Using default threshold of 0.5 (midpoint, conservative)"
        echo "Consider using a species-specific threshold (Human: 0.364, Mouse: 0.44, etc.)"
        FINAL_THRESHOLD=0.5
    fi
else
    echo "Using manually specified threshold: ${CPAT_THRESHOLD}"
    FINAL_THRESHOLD="${CPAT_THRESHOLD}"
fi

echo ""
echo "============================================================================"
echo "Step 6: Running CPAT with custom models"
echo "============================================================================"
echo ""

echo "CPAT threshold: ${FINAL_THRESHOLD}"
echo "Hexamer table: ${HEXAMER_TABLE}"
echo "Logit model: ${LOGIT_MODEL}"
echo ""

# Build CPAT command with custom models
CPAT_CMD="cpat.py -g ${INPUT_FASTA} --min-orf=50 --top-orf=50 -o ${OUTDIR}/coding_potential/cpat/cpat_output.txt"

# Add hexamer table
if [[ -f "${HEXAMER_TABLE}" ]]; then
    CPAT_CMD="${CPAT_CMD} -x ${HEXAMER_TABLE}"
fi

# Add logit model
if [[ -f "${LOGIT_MODEL}" ]]; then
    CPAT_CMD="${CPAT_CMD} -d ${LOGIT_MODEL}"
fi

# Run CPAT
echo "Running CPAT with custom models..."
echo "Command: ${CPAT_CMD}"
echo ""

eval "${CPAT_CMD}"

# Generate statistics
echo ""
echo "Generating CPAT statistics..."
python3 << EOF
import pandas as pd
import re

threshold = ${FINAL_THRESHOLD}

# Read CPAT results (using the best ORF probability file)
df = pd.read_csv('${OUTDIR}/coding_potential/cpat/cpat_output.txt.ORF_prob.best.tsv', sep='\t')

# Handle column name variations between CPAT versions
# The file has both 'seq_ID' (transcript ID) and 'ID' (ORF ID) columns
# We want to use 'seq_ID' which contains the transcript IDs (e.g., MSTRG.1.1)
# NOT 'ID' which contains ORF IDs (e.g., MSTRG.1.1_ORF_2)
if 'Coding_prob' in df.columns:
    df.rename(columns={'Coding_prob': 'coding_prob'}, inplace=True)
# Only rename 'ID' to 'seq_ID' if 'seq_ID' doesn't already exist
if 'ID' in df.columns and 'seq_ID' not in df.columns:
    df.rename(columns={'ID': 'seq_ID'}, inplace=True)

# Classify transcripts
df['classification'] = df['coding_prob'].apply(
    lambda x: 'coding' if x >= threshold else 'noncoding'
)

# Statistics
total = len(df)
noncoding = (df['classification'] == 'noncoding').sum()
coding = (df['classification'] == 'coding').sum()

# Write statistics
with open('${OUTDIR}/coding_potential/cpat/cpat_stats.txt', 'w') as out:
    out.write("CPAT Coding Potential Assessment (Custom Models)\n")
    out.write("=" * 60 + "\n\n")
    out.write(f"Models built from: ${REFERENCE_GTF}\n")
    out.write(f"Threshold: {threshold}\n")
    out.write(f"Total transcripts: {total}\n")
    out.write(f"Non-coding (< {threshold}): {noncoding} ({noncoding/total*100:.2f}%)\n")
    out.write(f"Coding (>= {threshold}): {coding} ({coding/total*100:.2f}%)\n\n")
    
    out.write("Score Distribution:\n")
    out.write(f"  Min:  {df['coding_prob'].min():.4f}\n")
    out.write(f"  Max:  {df['coding_prob'].max():.4f}\n")
    out.write(f"  Mean: {df['coding_prob'].mean():.4f}\n")
    out.write(f"  Median: {df['coding_prob'].median():.4f}\n")

# Write non-coding transcript IDs - ensure clean format
id_col = 'seq_ID' if 'seq_ID' in df.columns else df.columns[0]
noncoding_df = df[df['classification'] == 'noncoding']

# Clean transcript IDs: remove any commas, spaces, or special characters
# Extract only the transcript ID (first field if comma/tab/space-separated)
clean_ids = []
for seq_id in noncoding_df[id_col].values:
    seq_id_str = str(seq_id).strip()
    # Remove '>' prefix if present
    if seq_id_str.startswith('>'):
        seq_id_str = seq_id_str[1:].strip()
    # Split by any delimiter (comma, tab, or space) and take first field
    seq_id_str = re.split(r'[,\t\s]+', seq_id_str)[0]
    clean_ids.append(seq_id_str)

# Write clean IDs one per line
with open('${OUTDIR}/coding_potential/cpat/cpat_noncoding.txt', 'w') as out:
    for transcript_id in clean_ids:
        out.write(f"{transcript_id}\n")

print(f"CPAT completed: {noncoding} non-coding, {coding} coding (threshold: {threshold})")
EOF

printf "\n\nCompleted Job: 05a-cpat\n\n"
pwd
date
