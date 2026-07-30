#!/bin/bash

###############################################################################
# 13a_setup_blast.sh - Download and extract lncRNA BLAST database
# 
# This script downloads and extracts lncRNA databases for BLAST annotation.
# Run this script on a node with internet access BEFORE running 13_annotation.sh.
#
# Note: This script only downloads and extracts. Script 13 will handle indexing.
#
# Databases used:
#   - NONCODE: Comprehensive database of non-coding RNAs
#
# Other databases could be integrated as well, but many databases have sequences combined from all organisms and not species-specific. Therefore, it could be too many sequences to BLAST against taking very long time.
#
# Usage: bash 13a_setup_blast.sh
#
# Outputs:
#   ${BLAST_DB_DIR}/lncrna_database.fasta
###############################################################################

printf "\n\nStarted: 13a-setup_blast\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Directory to store BLAST database
BLAST_DB_DIR="/storage/users/ahinsu/lncPredict/databases/lncrna_blast"

# Species specific downlaod link given at http://www.noncode.org/download.php
DB_URL="http://www.noncode.org/datadownload/NONCODEv6_cow.fa.gz"


# ============================================================================
# Create database directory
# ============================================================================

mkdir -p "${BLAST_DB_DIR}"

echo ""
echo "lncRNA BLAST Database Setup"
echo "============================"
echo "Target directory: ${BLAST_DB_DIR}"
echo ""

# Check if database already exists
if [[ -f "${BLAST_DB_DIR}/lncrna_database.fasta" ]]; then
    echo "BLAST database FASTA file already exists."
    echo "If you want to re-download, delete ${BLAST_DB_DIR}/lncrna_database.fasta and run this script again."
    exit 0
fi

# ============================================================================
# Download lncRNA database based on source
# ============================================================================

cd "${BLAST_DB_DIR}" || exit 1

echo "Downloading NONCODE database..."
echo ""

echo "Downloading from: ${DB_URL}"
wget -c "${DB_URL}" -O noncode.fa.gz || {
    echo "ERROR: Failed to download NONCODE database"
    echo "The URL might have changed. Please check: http://www.noncode.org/download.php"
    exit 1
}

gunzip -f noncode.fa.gz
mv noncode.fa lncrna_database.fasta


# ============================================================================
# Verify download
# ============================================================================

echo ""
echo "Verifying download..."
echo ""

if [[ -f "lncrna_database.fasta" ]]; then
    SIZE=$(du -h "lncrna_database.fasta" | cut -f1)
    SEQ_COUNT=$(grep -c "^>" lncrna_database.fasta || echo "0")
    
    echo "  ✓ lncrna_database.fasta (${SIZE})"
    echo ""
    echo "SUCCESS: BLAST database download complete!"
    echo ""
    echo "Database: ${DATABASE_SOURCE}"
    echo "Sequences: ${SEQ_COUNT}"
    echo "Location: ${BLAST_DB_DIR}/lncrna_database.fasta"
    echo ""
    echo "Next step: Run 13_annotation.sh which will index the database if needed."
    echo ""
else
    echo "  ✗ lncrna_database.fasta - MISSING"
    echo ""
    echo "ERROR: Database file is missing"
    echo "Download incomplete"
    exit 1
fi

printf "\n\nCompleted: 13a-setup_blast\n\n"
pwd
date
