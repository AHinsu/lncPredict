#!/bin/bash
#SBATCH --job-name=03-gffcompare
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

###############################################################################
# 03_gffcompare.sh - Compare merged assembly with reference annotation
# 
# This script uses gffcompare to compare the merged transcriptome with
# the reference annotation and assign class codes to identify novel
# transcripts for lncRNA detection.
#
# Usage: sbatch 03_gffcompare.sh
#
# Outputs:
#   ${OUTDIR}/gffcompare/gffcmp.annotated.gtf
#   ${OUTDIR}/gffcompare/gffcmp.tracking
#   ${OUTDIR}/gffcompare/gffcmp.loci
#   ${OUTDIR}/gffcompare/gffcmp.stats
#   ${OUTDIR}/gffcompare/class_code_stats.txt
#   ${OUTDIR}/gffcompare/class_code_summary.json
###############################################################################

printf "\n\nStarted Job: 03-gffcompare\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Class codes to keep for lncRNA detection
# u: intergenic, x: antisense, i: intronic, o: overlap, j: multi-exon junction
CLASS_CODES="u,x,i,o,j"

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results" # Make sure to provide the same path for all scripts

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
mkdir -p "${OUTDIR}/gffcompare"
mkdir -p logs

# Check input file
MERGED_GTF="${OUTDIR}/stringtie_merge/merged.gtf"
if [[ ! -f "${MERGED_GTF}" ]]; then
    echo "ERROR: Merged GTF not found: ${MERGED_GTF}"
    echo "Please run 02_stringtie_merge.sh first"
    exit 1
fi

# Run gffcompare
echo "Running gffcompare..."
cd "${OUTDIR}/gffcompare"
gffcompare \
    -r "${REFERENCE_GTF}" \
    -o gffcmp \
    "${MERGED_GTF}"

echo "Generating class code statistics..."

# Generate class code statistics
python3 << EOF
import json
from collections import defaultdict

# Parse the annotated GTF to get class codes
class_counts = defaultdict(int)
transcript_classes = {}

with open('gffcmp.annotated.gtf', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9 or fields[2] != 'transcript':
            continue
        
        attributes = fields[8]
        class_code = None
        transcript_id = None
        
        for attr in attributes.split(';'):
            attr = attr.strip()
            if attr.startswith('class_code'):
                class_code = attr.split('"')[1]
            elif attr.startswith('transcript_id'):
                transcript_id = attr.split('"')[1]
        
        if class_code and transcript_id:
            class_counts[class_code] += 1
            transcript_classes[transcript_id] = class_code

# Class code descriptions
class_descriptions = {
    '=': 'Complete intron chain match',
    'c': 'Contained in reference',
    'j': 'Multi-exon with at least one junction match',
    'e': 'Single exon transfrag overlapping a reference exon',
    'i': 'Fully contained within reference intron',
    'o': 'Generic exonic overlap with reference',
    'p': 'Possible polymerase run-on',
    'r': 'Repeat (at least 50% overlap)',
    's': 'Antisense (intron match on opposite strand)',
    'u': 'Intergenic (novel lncRNA candidate)',
    'x': 'Antisense exonic overlap (novel lncRNA candidate)',
    'y': 'Exon overlap on opposite strand',
    '.': 'Tracking file only'
}

# Target class codes for lncRNA detection
target_classes = set('${CLASS_CODES}'.split(','))

# Write statistics
with open('class_code_stats.txt', 'w') as out:
    out.write("GffCompare Class Code Statistics\n")
    out.write("=" * 70 + "\n\n")
    out.write(f"Total transcripts: {sum(class_counts.values())}\n\n")
    
    out.write("Class Code Distribution:\n")
    out.write("-" * 70 + "\n")
    
    for code in sorted(class_counts.keys()):
        count = class_counts[code]
        desc = class_descriptions.get(code, 'Unknown')
        marker = " *TARGET*" if code in target_classes else ""
        out.write(f"{code:3s} : {count:7d} - {desc}{marker}\n")
    
    out.write("\n")
    out.write(f"Target class codes for lncRNA: {', '.join(sorted(target_classes))}\n")
    out.write(f"Total target transcripts: {sum(class_counts[c] for c in target_classes if c in class_counts)}\n")

# Write JSON summary
summary = {
    'total_transcripts': sum(class_counts.values()),
    'class_counts': dict(class_counts),
    'target_classes': list(target_classes),
    'target_count': sum(class_counts[c] for c in target_classes if c in class_counts),
    'class_descriptions': class_descriptions
}

with open('class_code_summary.json', 'w') as out:
    json.dump(summary, out, indent=2)

print(f"Class code statistics generated: {len(class_counts)} unique codes found")
print(f"Target transcripts (class codes {', '.join(sorted(target_classes))}): {summary['target_count']}")
EOF

printf "\n\nCompleted Job: 03-gffcompare\n\n"
pwd
date
