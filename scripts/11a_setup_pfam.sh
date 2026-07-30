#!/bin/bash

###############################################################################
# 11a_setup_pfam.sh - Download and extract Pfam database
# 
# This script downloads and extracts the Pfam database for pfam_scan.
# Run this script on a node with internet access BEFORE running 11_orf_pfam.sh.
#
# The Pfam database is large (~300-500 MB compressed, ~15 GB uncompressed) and
# requires internet access to download. The other option is UniProtKB but it is
# large database as well (60G), and sequences are provided for all organisms, 
# and not species-specific. 
#
# Note: This script only downloads and extracts. Script 11 will handle indexing.
#
# Usage: bash 11a_setup_pfam.sh
#
# Outputs:
#   ${PFAM_DB_DIR}/Pfam-A.hmm
#   ${PFAM_DB_DIR}/Pfam-A.hmm.dat
###############################################################################

printf "\n\nStarted: 11a-setup_pfam\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Directory to store Pfam database
PFAM_DB_DIR="/storage/users/ahinsu/lncPredict/databases/pfam"

# Pfam latest release will be downloaded. Modify the URL below for previous versions.

# ============================================================================
# Create database directory
# ============================================================================

mkdir -p "${PFAM_DB_DIR}"

echo ""
echo "Pfam Database Setup"
echo "===================="
echo "Target directory: ${PFAM_DB_DIR}"
echo ""

echo ""
echo "Checking for existing Pfam database..."
echo ""

# Required files
REQUIRED_FILES=("Pfam-A.hmm" "Pfam-A.hmm.dat")

MISSING=0

# Check each required file
for FILE in "${REQUIRED_FILES[@]}"; do
    if [[ -f "${PFAM_DB_DIR}/${FILE}" ]]; then
        SIZE=$(du -h "${PFAM_DB_DIR}/${FILE}" | cut -f1)
        echo "  ✓ ${FILE} already exists (${SIZE})"
    else
        echo "  ✗ ${FILE} is missing"
        MISSING=1
    fi
done

echo ""

# If everything exists, exit cleanly
if [[ $MISSING -eq 0 ]]; then
    echo "Pfam database already exists."
    echo "If you want to re-download, delete the files in:"
    echo "  ${PFAM_DB_DIR}/"
    echo "and run this script again."
    exit 0
fi

echo "Some files are missing — proceeding to download..."
echo ""

# ============================================================================
# Download Pfam-A HMM database
# ============================================================================

echo "Downloading Pfam-A HMM database..."
echo "This may take a while (file is ~400-500 MB compressed)..."
echo ""

cd "${PFAM_DB_DIR}" || exit 1

# ================================
# Pfam database download locations. THis is for latest release
# ================================
PFAM_URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_DAT_URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"

echo "Downloading Pfam HMM database..."
echo "  HMM: ${PFAM_URL}"
echo "  DAT: ${PFAM_DAT_URL}"

# ================================
# Download both files
# ================================
wget -c "${PFAM_URL}" || {
    echo "ERROR: Failed to download Pfam-A.hmm.gz"
    exit 1
}

wget -c "${PFAM_DAT_URL}" || {
    echo "ERROR: Failed to download Pfam-A.hmm.dat.gz"
    exit 1
}

# ================================
# Extract both files
# ================================
echo "Extracting Pfam-A.hmm.gz..."
gunzip -f Pfam-A.hmm.gz || {
    echo "ERROR: Failed to extract Pfam-A.hmm.gz"
    exit 1
}

echo "Extracting Pfam-A.hmm.dat.gz..."
gunzip -f Pfam-A.hmm.dat.gz || {
    echo "ERROR: Failed to extract Pfam-A.hmm.dat.gz"
    exit 1
}

# ============================================================================
# Verify download
# ============================================================================

echo ""
echo "Verifying download..."
echo ""

echo ""
echo "Verifying download..."
echo ""

# List of required files
REQUIRED_FILES=("Pfam-A.hmm" "Pfam-A.hmm.dat")

MISSING=0

for FILE in "${REQUIRED_FILES[@]}"; do
    if [[ -f "$FILE" ]]; then
        SIZE=$(du -h "$FILE" | cut -f1)
        echo "  ✓ $FILE ($SIZE)"
    else
        echo "  ✗ $FILE - MISSING"
        MISSING=1
    fi
done

echo ""

if [[ $MISSING -eq 0 ]]; then
    echo "SUCCESS: Pfam database download complete!"
    echo ""
    echo "Database location: ${PFAM_DB_DIR}/"
    echo ""
    echo "Next step: Run script 11_orf_pfam.sh which will index the database if needed."
    echo ""
else
    echo "ERROR: One or more required Pfam files are missing."
    echo "Download incomplete."
    exit 1
fi


printf "\n\nCompleted: 11a-setup_pfam\n\n"
pwd
date
