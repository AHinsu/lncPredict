#!/bin/bash
#SBATCH --job-name=04-filter_transcripts
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

###############################################################################
# 04_filter_transcripts.sh - Filter candidate lncRNA transcripts
# 
# This script filters transcripts based on class code, length, and exon
# count criteria to identify lncRNA candidates.
#
# Usage: sbatch 04_filter_transcripts.sh
#
# Outputs:
#   ${OUTDIR}/filtered_transcripts/filtered.gtf
#   ${OUTDIR}/filtered_transcripts/filtered.fasta
#   ${OUTDIR}/filtered_transcripts/filter_stats.txt
#   ${OUTDIR}/filtered_transcripts/filter_summary.json
#   ${OUTDIR}/filtered_transcripts/class_code_info.tsv
###############################################################################

printf "\n\nStarted Job: 04-filter_transcripts\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Reference genome FASTA file
REFERENCE_FASTA="/storage/users/ahinsu/lncPredict/reference/Bos_taurus.chromsomal.fa"

# Minimum transcript length (nucleotides)
MIN_LENGTH=200

# Minimum exon count(s): one value for all class codes (e.g., "2"),
# or one per CLASS_CODES entry in the same order (e.g., "2,1,1,1,2")
MIN_EXONS="1"

# u: intergenic, x: antisense, i: intronic, o: overlap, j: multi-exon junction
CLASS_CODES="u,x,i,o,j"

# Output directory
OUTDIR="/storage/users/ahinsu/lncPredict/results" # Make sure to provide the same path for all scripts

# Conda environment name (created from environment.yml)
CONDA_ENV="lncrna"

# Reference genome annotation GTF file
REFERENCE_GTF="/storage/users/ahinsu/lncPredict/reference/Bos_taurus.cleaned.gtf"



# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

# Create output directory
mkdir -p "${OUTDIR}/filtered_transcripts"
mkdir -p logs

# Check input files
ANNOTATED_GTF="${OUTDIR}/gffcompare/gffcmp.annotated.gtf"
if [[ ! -f "${ANNOTATED_GTF}" ]]; then
    echo "ERROR: Annotated GTF not found: ${ANNOTATED_GTF}"
    echo "Please run 03_gffcompare.sh first"
    exit 1
fi

echo "Filter criteria:"
echo "  Min length: ${MIN_LENGTH} nt"
echo "  Min exons: ${MIN_EXONS}"
echo "  Target class codes: ${CLASS_CODES}"
echo ""

# Filter transcripts
python3 << EOF
import json
from collections import defaultdict
import re

# Parameters
min_length = ${MIN_LENGTH}
min_exons_raw = '${MIN_EXONS}'
target_class_list = [code.strip() for code in '${CLASS_CODES}'.split(',') if code.strip()]

if not target_class_list:
    raise ValueError("CLASS_CODES cannot be empty")

if len(target_class_list) != len(set(target_class_list)):
    raise ValueError("CLASS_CODES contains duplicates; use unique class codes when mapping MIN_EXONS")

min_exons_values = [x.strip() for x in min_exons_raw.split(',') if x.strip()]
if not min_exons_values:
    raise ValueError("MIN_EXONS cannot be empty")

try:
    min_exons_values = [int(x) for x in min_exons_values]
except ValueError:
    raise ValueError(f"MIN_EXONS must contain integer value(s), got: {min_exons_raw}")

if any(x < 1 for x in min_exons_values):
    raise ValueError(f"MIN_EXONS values must be >= 1, got: {min_exons_raw}")

if len(min_exons_values) == 1:
    class_min_exons = {cc: min_exons_values[0] for cc in target_class_list}
elif len(min_exons_values) == len(target_class_list):
    class_min_exons = dict(zip(target_class_list, min_exons_values))
else:
    raise ValueError(
        "MIN_EXONS must have either 1 value or the same number of values as CLASS_CODES "
        f"(got {len(min_exons_values)} MIN_EXONS values and {len(target_class_list)} CLASS_CODES values)"
    )

target_classes = set(target_class_list)

print(f"Filtering with parameters:")
print(f"  Min length: {min_length} nt")
print(f"  Min exons (raw): {min_exons_raw}")
print(f"  Target class codes: {', '.join(sorted(target_classes))}")
print("  Min exons by class code:")
for class_code in target_class_list:
    print(f"    {class_code}: >= {class_min_exons[class_code]} exons")

# Parse GTF and filter transcripts
transcripts = {}
transcript_info = {}
exon_counts = defaultdict(int)

with open('${ANNOTATED_GTF}', 'r') as f:
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
        strand = fields[6]
        attributes = fields[8]
        
        # Parse attributes
        attr_dict = {}
        for attr in attributes.split(';'):
            attr = attr.strip()
            if not attr:
                continue
            match = re.match(r'(\S+)\s+"([^"]+)"', attr)
            if match:
                key, value = match.groups()
                attr_dict[key] = value
        
        transcript_id = attr_dict.get('transcript_id')
        gene_id = attr_dict.get('gene_id')
        class_code = attr_dict.get('class_code', '.')
        
        if not transcript_id:
            continue
        
        # Store transcript info
        if feature == 'transcript':
            length = end - start + 1
            transcript_info[transcript_id] = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'length': length,
                'class_code': class_code,
                'gene_id': gene_id,
                'line': line
            }
        
        # Count exons
        if feature == 'exon':
            exon_counts[transcript_id] += 1
            
            if transcript_id not in transcripts:
                transcripts[transcript_id] = []
            transcripts[transcript_id].append({
                'start': start,
                'end': end
            })

# Filter transcripts
filtered_transcripts = {}
filter_reasons = defaultdict(int)
class_code_counts = {'before_filter': defaultdict(int), 'after_filter': defaultdict(int)}

for tid, info in transcript_info.items():
    class_code_counts['before_filter'][info['class_code']] += 1
    
    # Check class code
    if info['class_code'] not in target_classes:
        filter_reasons['wrong_class_code'] += 1
        continue
    
    # Check length
    if info['length'] < min_length:
        filter_reasons['too_short'] += 1
        continue
    
    # Check exon count
    required_min_exons = class_min_exons[info['class_code']]
    if exon_counts[tid] < required_min_exons:
        filter_reasons['too_few_exons'] += 1
        continue
    
    # Passed all filters
    filtered_transcripts[tid] = info
    class_code_counts['after_filter'][info['class_code']] += 1

print(f"\nFiltering results:")
print(f"  Total transcripts: {len(transcript_info)}")
print(f"  Filtered transcripts: {len(filtered_transcripts)}")
print(f"  Removed: {len(transcript_info) - len(filtered_transcripts)}")

# Write filtered GTF
with open('${OUTDIR}/filtered_transcripts/filtered.gtf', 'w') as out:
    out.write("# Filtered transcripts for lncRNA detection\n")
    exons_by_class = ",".join(f"{cc}:{class_min_exons[cc]}" for cc in target_class_list)
    out.write(f"# Filters: length >= {min_length}, exons by class {{{exons_by_class}}}, class codes: {','.join(sorted(target_classes))}\n")
    
    with open('${ANNOTATED_GTF}', 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            attributes = fields[8]
            transcript_id = None
            
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    transcript_id = attr.split('"')[1]
                    break
            
            if transcript_id and transcript_id in filtered_transcripts:
                out.write(line)

# Write class code info
with open('${OUTDIR}/filtered_transcripts/class_code_info.tsv', 'w') as out:
    out.write("transcript_id\tclass_code\tlength\texon_count\tchrom\tstart\tend\tstrand\n")
    for tid, info in filtered_transcripts.items():
        out.write(f"{tid}\t{info['class_code']}\t{info['length']}\t{exon_counts[tid]}\t")
        out.write(f"{info['chrom']}\t{info['start']}\t{info['end']}\t{info['strand']}\n")

# Write statistics
with open('${OUTDIR}/filtered_transcripts/filter_stats.txt', 'w') as out:
    out.write("Transcript Filtering Statistics\n")
    out.write("=" * 70 + "\n\n")
    out.write(f"Input transcripts: {len(transcript_info)}\n")
    out.write(f"Filtered transcripts: {len(filtered_transcripts)}\n")
    out.write(f"Removal rate: {((len(transcript_info) - len(filtered_transcripts)) / len(transcript_info) * 100):.2f}%\n\n")
    
    out.write("Filter Criteria:\n")
    out.write(f"  Minimum length: {min_length} nt\n")
    out.write(f"  Minimum exons (raw): {min_exons_raw}\n")
    out.write("  Minimum exons by class code:\n")
    for class_code in target_class_list:
        out.write(f"    {class_code}: {class_min_exons[class_code]}\n")
    out.write(f"  Target class codes: {', '.join(sorted(target_classes))}\n\n")
    
    out.write("Removal Reasons:\n")
    for reason, count in sorted(filter_reasons.items()):
        out.write(f"  {reason}: {count}\n")
    
    out.write("\nClass Code Distribution (Before Filter):\n")
    for code, count in sorted(class_code_counts['before_filter'].items()):
        out.write(f"  {code}: {count}\n")
    
    out.write("\nClass Code Distribution (After Filter):\n")
    for code, count in sorted(class_code_counts['after_filter'].items()):
        out.write(f"  {code}: {count}\n")

# Write JSON summary
summary = {
    'input_transcripts': len(transcript_info),
    'filtered_transcripts': len(filtered_transcripts),
    'removal_rate': ((len(transcript_info) - len(filtered_transcripts)) / len(transcript_info) * 100) if transcript_info else 0,
    'filter_criteria': {
        'min_length': min_length,
        'min_exons_raw': min_exons_raw,
        'min_exons_by_class': class_min_exons,
        'target_classes': list(target_classes)
    },
    'removal_reasons': dict(filter_reasons),
    'class_code_distribution': {
        'before': dict(class_code_counts['before_filter']),
        'after': dict(class_code_counts['after_filter'])
    }
}

with open('${OUTDIR}/filtered_transcripts/filter_summary.json', 'w') as out:
    json.dump(summary, out, indent=2)

print("\nFiltering completed successfully")
EOF

# Extract sequences from reference genome
echo ""
echo "Extracting sequences from reference genome..."
gffread "${OUTDIR}/filtered_transcripts/filtered.gtf" \
    -g "${REFERENCE_FASTA}" \
    -w "${OUTDIR}/filtered_transcripts/filtered.fasta"

# Count sequences
SEQ_COUNT=$(grep -c "^>" "${OUTDIR}/filtered_transcripts/filtered.fasta" || echo "0")
echo "Extracted ${SEQ_COUNT} sequences"

printf "\n\nCompleted Job: 04-filter_transcripts\n\n"
pwd
date
