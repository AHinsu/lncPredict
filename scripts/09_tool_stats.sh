#!/bin/bash
#SBATCH --job-name=09-tool_stats
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

###############################################################################
# 09_tool_stats.sh - Combined and overlapping statistics for coding potential tools
# 
# This script collates results from all coding potential assessment tools
# (CPAT, CPC2, FEELnc, PLEKv2) and generates comprehensive statistics:
#   - Per-tool statistics overall and by class code
#   - Combined statistics showing tool agreement and overlap
#   - Venn diagrams (2-3 tools) or UpSet plots (4+ tools) for visual comparison
#
# Statistics are based on sequences present in the FASTA file (not GTF),
# as tools use FASTA as input and not all transcripts have exons.
#
# This script helps decide which transcripts to take forward by providing
# detailed analysis of tool consensus before the final consensus step (script 10).
#
# Usage: sbatch 09_tool_stats.sh
#
# Outputs:
#   ${OUTDIR}/coding_potential/final_stats/overall_stats.txt
#   ${OUTDIR}/coding_potential/final_stats/per_tool_stats.txt
#   ${OUTDIR}/coding_potential/final_stats/per_class_code_stats.txt
#   ${OUTDIR}/coding_potential/final_stats/overlap_matrix.tsv
#   ${OUTDIR}/coding_potential/final_stats/detailed_overlap.tsv
#   ${OUTDIR}/coding_potential/final_stats/summary_report.txt
#   ${OUTDIR}/coding_potential/final_stats/venn_diagrams/venn_overall.png
#   ${OUTDIR}/coding_potential/final_stats/venn_diagrams/venn_by_class/venn_*.png
###############################################################################

printf "\n\nStarted Job: 09-tool_stats\n\n"
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

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

# Create output directories
mkdir -p "${OUTDIR}/coding_potential/final_stats/venn_diagrams/venn_by_class"
mkdir -p logs

# ============================================================================
# Check input files
# ============================================================================

echo "Checking input files..."

# Class code information from script 4
CLASS_CODE_INFO="${OUTDIR}/filtered_transcripts/class_code_info.tsv"
if [[ ! -f "${CLASS_CODE_INFO}" ]]; then
    echo "ERROR: Class code info not found: ${CLASS_CODE_INFO}"
    echo "Please run 04_filter_transcripts.sh first"
    exit 1
fi

# Tool output files
CPAT_NONCODING="${OUTDIR}/coding_potential/cpat/cpat_noncoding.txt"
CPC2_NONCODING="${OUTDIR}/coding_potential/cpc2/cpc2_noncoding.txt"
FEELNC_NONCODING="${OUTDIR}/coding_potential/feelnc/feelnc_noncoding.txt"
PLEK2_NONCODING="${OUTDIR}/coding_potential/plek2/plek2_noncoding.txt"

# Check which tools have completed
TOOLS_AVAILABLE=()
TOOLS_MISSING=()

if [[ -f "${CPAT_NONCODING}" ]]; then
    TOOLS_AVAILABLE+=("CPAT")
    echo "✓ CPAT results found"
else
    TOOLS_MISSING+=("CPAT")
    echo "✗ CPAT results not found: ${CPAT_NONCODING}"
fi

if [[ -f "${CPC2_NONCODING}" ]]; then
    TOOLS_AVAILABLE+=("CPC2")
    echo "✓ CPC2 results found"
else
    TOOLS_MISSING+=("CPC2")
    echo "✗ CPC2 results not found: ${CPC2_NONCODING}"
fi

if [[ -f "${FEELNC_NONCODING}" ]]; then
    TOOLS_AVAILABLE+=("FEELnc")
    echo "✓ FEELnc results found"
else
    TOOLS_MISSING+=("FEELnc")
    echo "✗ FEELnc results not found: ${FEELNC_NONCODING}"
fi

if [[ -f "${PLEK2_NONCODING}" ]]; then
    TOOLS_AVAILABLE+=("PLEK2")
    echo "✓ PLEK2 results found"
else
    TOOLS_MISSING+=("PLEK2")
    echo "✗ PLEK2 results not found: ${PLEK2_NONCODING}"
fi

# Check if we have at least 2 tools
NUM_TOOLS=${#TOOLS_AVAILABLE[@]}
if [[ ${NUM_TOOLS} -lt 2 ]]; then
    echo ""
    echo "ERROR: At least 2 coding potential tools are required for meaningful statistics"
    echo "Available tools: ${TOOLS_AVAILABLE[*]}"
    echo "Missing tools: ${TOOLS_MISSING[*]}"
    echo ""
    echo "Please run at least 2 of the following scripts:"
    echo "  - 05a_cpat.sh or 05b_cpat.sh (CPAT)"
    echo "  - 06_cpc2.sh (CPC2)"
    echo "  - 07_feelnc.sh (FEELnc)"
    echo "  - 08_plek.sh (PLEK2)"
    exit 1
fi

echo ""
echo "Proceeding with ${NUM_TOOLS} tools: ${TOOLS_AVAILABLE[*]}"
echo ""

# ============================================================================
# Generate comprehensive statistics
# ============================================================================

echo "Generating comprehensive statistics..."
python3 << EOF
import os
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles

# Try to import upsetplot for 4-way comparisons
try:
    from upsetplot import UpSet, from_memberships
    UPSET_AVAILABLE = True
except ImportError:
    UPSET_AVAILABLE = False
    print("Note: upsetplot not available. Will use alternative visualization for 4 tools.")

# Configuration - using bash variable substitution
outdir = '${OUTDIR}'
class_code_file = os.path.join(outdir, 'filtered_transcripts', 'class_code_info.tsv')
fasta_file = os.path.join(outdir, 'filtered_transcripts', 'filtered.fasta')

# Tool files
tool_files = {
    'CPAT': os.path.join(outdir, 'coding_potential/cpat/cpat_noncoding.txt'),
    'CPC2': os.path.join(outdir, 'coding_potential/cpc2/cpc2_noncoding.txt'),
    'FEELnc': os.path.join(outdir, 'coding_potential/feelnc/feelnc_noncoding.txt'),
    'PLEK2': os.path.join(outdir, 'coding_potential/plek2/plek2_noncoding.txt')
}

# Load sequences from FASTA file (these are the actual sequences used by tools)
print("Loading sequences from FASTA file...")
fasta_sequences = set()
if os.path.exists(fasta_file):
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract transcript ID (remove '>' and any trailing description)
                seq_id = line[1:].strip().split()[0]
                fasta_sequences.add(seq_id)
    print(f"Total sequences in FASTA file: {len(fasta_sequences)}")
else:
    print(f"ERROR: FASTA file not found: {fasta_file}")
    exit(1)

# Load class code information
print("Loading class code information...")
class_df = pd.read_csv(class_code_file, sep='\t')
print(f"Total transcripts in class_code_info.tsv: {len(class_df)}")

# Filter class_df to only include transcripts present in FASTA file
class_df = class_df[class_df['transcript_id'].isin(fasta_sequences)]
print(f"Transcripts with sequences in FASTA (after filtering): {len(class_df)}")

# Create a mapping of transcript_id to class_code
transcript_to_class = dict(zip(class_df['transcript_id'], class_df['class_code']))

# Load non-coding transcripts from each tool
tool_noncoding = {}
for tool, filepath in tool_files.items():
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            # Read transcript IDs and remove any whitespace
            transcripts = set(line.strip() for line in f if line.strip())
        tool_noncoding[tool] = transcripts
        print(f"{tool}: {len(transcripts)} non-coding transcripts")
    else:
        print(f"{tool}: Not available (file not found)")

# Check if we have at least 2 tools
available_tools = list(tool_noncoding.keys())
num_tools = len(available_tools)

if num_tools < 2:
    print(f"\nERROR: Only {num_tools} tool(s) available. Need at least 2 for overlap analysis.")
    exit(1)

print(f"\nProceeding with {num_tools} tools: {', '.join(available_tools)}")

# ============================================================================
# Generate Overall Statistics
# ============================================================================

print("\n" + "="*70)
print("OVERALL STATISTICS")
print("="*70)

overall_stats_file = os.path.join(outdir, 'coding_potential/final_stats/overall_stats.txt')
with open(overall_stats_file, 'w') as out:
    out.write("Overall Statistics for Coding Potential Assessment Tools\n")
    out.write("="*70 + "\n\n")
    
    # Total transcripts
    total_transcripts = len(class_df)
    out.write(f"Total filtered transcripts (with sequences in FASTA): {total_transcripts}\n")
    out.write(f"Note: Statistics based on sequences in filtered.fasta\n\n")
    
    # Per-tool summary
    out.write("Non-coding Transcripts per Tool:\n")
    out.write("-"*70 + "\n")
    for tool in available_tools:
        count = len(tool_noncoding[tool])
        pct = (count / total_transcripts * 100) if total_transcripts > 0 else 0
        out.write(f"  {tool:10s}: {count:6d} ({pct:5.2f}%)\n")
    
    # Union and intersection
    all_noncoding = set()
    for transcripts in tool_noncoding.values():
        all_noncoding.update(transcripts)
    
    # Start with first tool, then intersect with all others
    common_noncoding = tool_noncoding[available_tools[0]].copy()
    for tool in available_tools[1:]:
        common_noncoding &= tool_noncoding[tool]
    
    out.write("\nCombined Statistics:\n")
    out.write("-"*70 + "\n")
    union_count = len(all_noncoding)
    union_pct = (union_count / total_transcripts * 100) if total_transcripts > 0 else 0
    out.write(f"  Union (any tool):      {union_count:6d} ({union_pct:5.2f}%)\n")
    
    intersect_count = len(common_noncoding)
    intersect_pct = (intersect_count / total_transcripts * 100) if total_transcripts > 0 else 0
    out.write(f"  Intersection (all):    {intersect_count:6d} ({intersect_pct:5.2f}%)\n")
    
    # Agreement statistics (for all tools)
    out.write("\nAgreement Analysis:\n")
    out.write("-"*70 + "\n")
    
    # Count how many tools agree on each transcript
    vote_counts = defaultdict(int)
    for transcript in all_noncoding:
        votes = sum(1 for tool_set in tool_noncoding.values() if transcript in tool_set)
        vote_counts[votes] += 1
    
    for i in range(num_tools, 0, -1):
        count = vote_counts.get(i, 0)
        pct = (count / union_count * 100) if union_count > 0 else 0
        out.write(f"  {i}/{num_tools} tools agree:     {count:6d} ({pct:5.2f}% of union)\n")

print(f"Overall statistics written to: {overall_stats_file}")

# ============================================================================
# Generate Per-Tool Statistics
# ============================================================================

print("\n" + "="*70)
print("PER-TOOL STATISTICS")
print("="*70)

per_tool_stats_file = os.path.join(outdir, 'coding_potential/final_stats/per_tool_stats.txt')
with open(per_tool_stats_file, 'w') as out:
    out.write("Per-Tool Detailed Statistics\n")
    out.write("="*70 + "\n\n")
    
    for tool in available_tools:
        out.write(f"{tool}\n")
        out.write("-"*70 + "\n")
        
        transcripts = tool_noncoding[tool]
        total = len(transcripts)
        
        # Class code distribution for this tool
        class_counts = Counter()
        for tid in transcripts:
            class_code = transcript_to_class.get(tid, 'unknown')
            class_counts[class_code] += 1
        
        out.write(f"Total non-coding transcripts: {total}\n\n")
        out.write("Distribution by class code:\n")
        
        for class_code in sorted(class_counts.keys()):
            count = class_counts[class_code]
            pct = (count / total * 100) if total > 0 else 0
            out.write(f"  {class_code:3s}: {count:6d} ({pct:5.2f}%)\n")
        
        out.write("\n")

print(f"Per-tool statistics written to: {per_tool_stats_file}")

# ============================================================================
# Generate Per-Class-Code Statistics
# ============================================================================

print("\n" + "="*70)
print("PER-CLASS-CODE STATISTICS")
print("="*70)

# Get all unique class codes
all_class_codes = sorted(class_df['class_code'].unique())

per_class_stats_file = os.path.join(outdir, 'coding_potential/final_stats/per_class_code_stats.txt')
with open(per_class_stats_file, 'w') as out:
    out.write("Per-Class-Code Statistics Across Tools\n")
    out.write("="*70 + "\n\n")
    
    for class_code in all_class_codes:
        out.write(f"Class Code: {class_code}\n")
        out.write("-"*70 + "\n")
        
        # Get all transcripts with this class code
        class_transcripts = set(class_df[class_df['class_code'] == class_code]['transcript_id'])
        total_in_class = len(class_transcripts)
        
        out.write(f"Total transcripts: {total_in_class}\n\n")
        
        # Per-tool counts for this class code
        out.write("Non-coding predictions per tool:\n")
        
        tool_class_counts = {}
        for tool in available_tools:
            # Intersection of tool's non-coding with this class
            noncoding_in_class = tool_noncoding[tool] & class_transcripts
            count = len(noncoding_in_class)
            pct = (count / total_in_class * 100) if total_in_class > 0 else 0
            tool_class_counts[tool] = noncoding_in_class
            out.write(f"  {tool:10s}: {count:6d} ({pct:5.2f}%)\n")
        
        # Union and intersection for this class
        class_union = set()
        for transcripts in tool_class_counts.values():
            class_union.update(transcripts)
        
        class_intersect = tool_class_counts[available_tools[0]].copy()
        for tool in available_tools[1:]:
            class_intersect &= tool_class_counts[tool]
        
        out.write("\nCombined (for this class code):\n")
        union_count = len(class_union)
        union_pct = (union_count / total_in_class * 100) if total_in_class > 0 else 0
        out.write(f"  Union (any tool):      {union_count:6d} ({union_pct:5.2f}%)\n")
        
        intersect_count = len(class_intersect)
        intersect_pct = (intersect_count / total_in_class * 100) if total_in_class > 0 else 0
        out.write(f"  Intersection (all):    {intersect_count:6d} ({intersect_pct:5.2f}%)\n")
        
        out.write("\n" + "="*70 + "\n\n")

print(f"Per-class-code statistics written to: {per_class_stats_file}")

# ============================================================================
# Generate Overlap Matrix
# ============================================================================

print("\n" + "="*70)
print("GENERATING OVERLAP MATRIX")
print("="*70)

# Create pairwise overlap matrix
overlap_matrix = pd.DataFrame(index=available_tools, columns=available_tools, dtype=int)

for tool1 in available_tools:
    for tool2 in available_tools:
        overlap = len(tool_noncoding[tool1] & tool_noncoding[tool2])
        overlap_matrix.loc[tool1, tool2] = overlap

overlap_matrix_file = os.path.join(outdir, 'coding_potential/final_stats/overlap_matrix.tsv')
overlap_matrix.to_csv(overlap_matrix_file, sep='\t')
print(f"Overlap matrix written to: {overlap_matrix_file}")

# ============================================================================
# Generate Detailed Overlap Table
# ============================================================================

print("\n" + "="*70)
print("GENERATING DETAILED OVERLAP TABLE")
print("="*70)

# Create a detailed table with all transcripts and tool assignments
detailed_data = []

for transcript in all_noncoding:
    row = {
        'transcript_id': transcript,
        'class_code': transcript_to_class.get(transcript, 'unknown')
    }
    
    # Add column for each tool
    for tool in available_tools:
        row[f'{tool}_noncoding'] = 'Yes' if transcript in tool_noncoding[tool] else 'No'
    
    # Count total votes
    votes = sum(1 for tool in available_tools if transcript in tool_noncoding[tool])
    row['total_votes'] = votes
    row['agreement_level'] = f'{votes}/{num_tools}'
    
    detailed_data.append(row)

detailed_df = pd.DataFrame(detailed_data)

# Sort by total votes (descending) and then by transcript_id
detailed_df = detailed_df.sort_values(['total_votes', 'transcript_id'], ascending=[False, True])

detailed_overlap_file = os.path.join(outdir, 'coding_potential/final_stats/detailed_overlap.tsv')
detailed_df.to_csv(detailed_overlap_file, sep='\t', index=False)
print(f"Detailed overlap table written to: {detailed_overlap_file}")

# ============================================================================
# Generate Venn Diagrams - Overall
# ============================================================================

print("\n" + "="*70)
print("GENERATING VENN DIAGRAMS - OVERALL")
print("="*70)

def create_venn_or_upset_diagram(sets_dict, title, output_file):
    """Create Venn diagram for 2-3 sets, or UpSet plot for 4+ sets"""
    tools = list(sets_dict.keys())
    num_sets = len(tools)
    
    if num_sets == 2:
        plt.figure(figsize=(10, 8))
        venn = venn2([sets_dict[tools[0]], sets_dict[tools[1]]], 
                     set_labels=tools)
        venn2_circles([sets_dict[tools[0]], sets_dict[tools[1]]])
        plt.title(title, fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Venn diagram (2-way) saved: {output_file}")
        
    elif num_sets == 3:
        plt.figure(figsize=(10, 8))
        venn = venn3([sets_dict[tools[0]], sets_dict[tools[1]], sets_dict[tools[2]]], 
                     set_labels=tools)
        venn3_circles([sets_dict[tools[0]], sets_dict[tools[1]], sets_dict[tools[2]]])
        plt.title(title, fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Venn diagram (3-way) saved: {output_file}")
        
    elif num_sets >= 4:
        # For 4+ sets, use UpSet plot if available, otherwise create a bar chart
        if UPSET_AVAILABLE:
            # Create membership data for UpSet plot
            memberships = []
            for transcript in set.union(*sets_dict.values()):
                membership = tuple(transcript in sets_dict[tool] for tool in tools)
                memberships.append(membership)
            
            # Create UpSet plot
            upset_data = from_memberships(memberships)
            
            fig = plt.figure(figsize=(12, 8))
            upset = UpSet(upset_data, subset_size='count', show_counts=True,
                         element_size=40, intersection_plot_elements=6)
            # Highlight the complete intersection (all tools agree)
            upset.style_subsets(present=tools, facecolor="steelblue", label=f"All {num_sets} agree")
            upset.plot(fig=fig)
            plt.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  UpSet plot ({num_sets}-way) saved: {output_file}")
        else:
            # Fallback: Create a combination bar chart showing all intersections
            print(f"  Note: UpSet plot not available. Creating combination bar chart for {num_sets} tools.")
            
            # Calculate all possible intersections
            from itertools import combinations
            intersection_data = []
            
            # For each possible combination of tools
            for r in range(num_sets, 0, -1):
                for combo in combinations(tools, r):
                    # Start with the first tool's set
                    intersection = sets_dict[combo[0]].copy()
                    # Intersect with all other tools in combo
                    for tool in combo[1:]:
                        intersection &= sets_dict[tool]
                    # Subtract elements that appear in other tools not in combo
                    for tool in tools:
                        if tool not in combo:
                            intersection -= sets_dict[tool]
                    
                    if len(intersection) > 0:
                        intersection_data.append({
                            'combination': ' & '.join(combo),
                            'count': len(intersection),
                            'num_tools': len(combo)
                        })
            
            # Sort by number of tools (descending), then by count (descending)
            intersection_data.sort(key=lambda x: (x['num_tools'], x['count']), reverse=True)
            
            # Create bar chart
            fig, ax = plt.subplots(figsize=(14, 8))
            combos = [d['combination'] for d in intersection_data]
            counts = [d['count'] for d in intersection_data]
            colors = plt.cm.viridis([d['num_tools']/num_sets for d in intersection_data])
            
            bars = ax.barh(range(len(combos)), counts, color=colors)
            ax.set_yticks(range(len(combos)))
            ax.set_yticklabels(combos, fontsize=9)
            ax.set_xlabel('Number of Transcripts', fontsize=12)
            ax.set_title(title, fontsize=14, fontweight='bold')
            ax.invert_yaxis()
            
            # Add count labels on bars
            for i, (bar, count) in enumerate(zip(bars, counts)):
                ax.text(count, i, f' {count}', va='center', fontsize=9)
            
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  Combination bar chart ({num_sets}-way) saved: {output_file}")

# Overall Venn diagram or UpSet plot
venn_overall_file = os.path.join(outdir, 'coding_potential/final_stats/venn_diagrams/venn_overall.png')
create_venn_or_upset_diagram(
    tool_noncoding,
    f"Non-coding Transcript Overlap - All Tools (n={len(all_noncoding)})",
    venn_overall_file
)

# ============================================================================
# Generate Venn Diagrams - Per Class Code
# ============================================================================

print("\n" + "="*70)
print("GENERATING VENN DIAGRAMS - PER CLASS CODE")
print("="*70)

for class_code in all_class_codes:
    # Get all transcripts with this class code
    class_transcripts = set(class_df[class_df['class_code'] == class_code]['transcript_id'])
    
    # Get non-coding transcripts per tool for this class code
    tool_class_sets = {}
    for tool in available_tools:
        noncoding_in_class = tool_noncoding[tool] & class_transcripts
        tool_class_sets[tool] = noncoding_in_class
    
    # Only create Venn diagram if there are non-coding transcripts
    total_noncoding = len(set.union(*tool_class_sets.values())) if tool_class_sets else 0
    
    if total_noncoding > 0:
        venn_class_file = os.path.join(
            outdir, 
            f'coding_potential/final_stats/venn_diagrams/venn_by_class/venn_class_{class_code}.png'
        )
        
        create_venn_or_upset_diagram(
            tool_class_sets,
            f"Non-coding Overlap - Class {class_code} (n={total_noncoding})",
            venn_class_file
        )
    else:
        print(f"  Skipping class {class_code} - no non-coding transcripts")

# ============================================================================
# Generate Summary Report
# ============================================================================

print("\n" + "="*70)
print("GENERATING SUMMARY REPORT")
print("="*70)

summary_file = os.path.join(outdir, 'coding_potential/final_stats/summary_report.txt')
with open(summary_file, 'w') as out:
    out.write("="*70 + "\n")
    out.write("CODING POTENTIAL TOOL STATISTICS - SUMMARY REPORT\n")
    out.write("="*70 + "\n\n")
    
    out.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    out.write(f"Tools analyzed: {', '.join(available_tools)}\n")
    out.write(f"Total filtered transcripts (with sequences in FASTA): {len(class_df)}\n")
    out.write(f"Note: Statistics based on sequences in filtered.fasta\n\n")
    
    out.write("="*70 + "\n")
    out.write("KEY FINDINGS\n")
    out.write("="*70 + "\n\n")
    
    # Tool agreement
    out.write("1. Tool Agreement on Non-coding Classification:\n\n")
    for i in range(num_tools, 0, -1):
        count = vote_counts.get(i, 0)
        pct = (count / len(all_noncoding) * 100) if len(all_noncoding) > 0 else 0
        out.write(f"   {i}/{num_tools} tools agree: {count:6d} transcripts ({pct:5.2f}%)\n")
    
    out.write("\n")
    
    # Class code distribution
    out.write("2. Class Code Distribution (all non-coding by union):\n\n")
    union_class_counts = Counter()
    for tid in all_noncoding:
        class_code = transcript_to_class.get(tid, 'unknown')
        union_class_counts[class_code] += 1
    
    for class_code in sorted(union_class_counts.keys()):
        count = union_class_counts[class_code]
        pct = (count / len(all_noncoding) * 100) if len(all_noncoding) > 0 else 0
        out.write(f"   {class_code:3s}: {count:6d} ({pct:5.2f}%)\n")
    
    out.write("\n")
    
    # Recommendations
    out.write("="*70 + "\n")
    out.write("RECOMMENDATIONS\n")
    out.write("="*70 + "\n\n")
    
    out.write(f"Based on {num_tools} tools, consider the following consensus strategies:\n\n")
    
    # Calculate transcripts for different strategies
    all_agree = vote_counts.get(num_tools, 0)
    # Majority means more than half: for 4 tools, that's 3 or 4; for 3 tools, that's 2 or 3
    majority_threshold = (num_tools // 2) + 1
    majority = sum(vote_counts.get(i, 0) for i in range(majority_threshold, num_tools + 1))
    any_two = sum(vote_counts.get(i, 0) for i in range(2, num_tools + 1))
    
    out.write(f"  - Conservative (all {num_tools} tools agree):  {all_agree:6d} transcripts\n")
    out.write(f"  - Moderate (majority agree):      {majority:6d} transcripts\n")
    out.write(f"  - Liberal (any 2+ tools agree):   {any_two:6d} transcripts\n")
    out.write(f"  - Very liberal (union, any tool): {len(all_noncoding):6d} transcripts\n\n")
    
    out.write("Use script 10 (10_coding_consensus.sh) to apply your chosen strategy.\n\n")
    
    out.write("="*70 + "\n")
    out.write("OUTPUT FILES\n")
    out.write("="*70 + "\n\n")
    
    out.write(f"  - Overall statistics:     {overall_stats_file}\n")
    out.write(f"  - Per-tool statistics:    {per_tool_stats_file}\n")
    out.write(f"  - Per-class statistics:   {per_class_stats_file}\n")
    out.write(f"  - Overlap matrix:         {overlap_matrix_file}\n")
    out.write(f"  - Detailed overlap table: {detailed_overlap_file}\n")
    out.write(f"  - Venn/UpSet plot (overall): {venn_overall_file}\n")
    out.write(f"  - Venn/UpSet plots (by class): {outdir}/coding_potential/final_stats/venn_diagrams/venn_by_class/\n")
    out.write("\n")

print(f"\nSummary report written to: {summary_file}")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print(f"\nAll output files written to:")
print(f"  {outdir}/coding_potential/final_stats/")
print("\n")

EOF

echo ""
echo "Statistics generation completed successfully!"
echo ""
echo "Output files:"
echo "  - ${OUTDIR}/coding_potential/final_stats/overall_stats.txt"
echo "  - ${OUTDIR}/coding_potential/final_stats/per_tool_stats.txt"
echo "  - ${OUTDIR}/coding_potential/final_stats/per_class_code_stats.txt"
echo "  - ${OUTDIR}/coding_potential/final_stats/overlap_matrix.tsv"
echo "  - ${OUTDIR}/coding_potential/final_stats/detailed_overlap.tsv"
echo "  - ${OUTDIR}/coding_potential/final_stats/summary_report.txt"
echo "  - ${OUTDIR}/coding_potential/final_stats/venn_diagrams/venn_overall.png"
echo "  - ${OUTDIR}/coding_potential/final_stats/venn_diagrams/venn_by_class/*.png"
echo ""

printf "\n\nCompleted Job: 09-tool_stats\n\n"
pwd
date
