#!/bin/bash

# This script downloads and installs hh-suite 
# Please visit https://mmseqs.com/hhsuite/ and choose the appropriate compiled 
# version that is compatible with your system. 
# 
# Usage:
# cd /path/to/VirusPPIScreen/
# chmod +x scripts/download_dependencies.sh
# ./scripts/download_dependencies.sh <file name of precompiled binary on https://mmseqs.com/hhsuite/>
# Example:
# ./scripts/download_dependencies.sh hhsuite-linux-sse2.tar.gz

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <file name to hhsuite precompiled binary .tar.gz file 
    on https://mmseqs.com/hhsuite/ compatible with your system. example: 'hhsuite-linux-sse2.tar.gz'"
    exit 1
fi

echo "https://mmseqs.com/hhsuite/${1}"
wget -O "$1" "https://mmseqs.com/hhsuite/${1}"

if [ $? -eq 0 ]; then
    echo "Download completed successfully."
else
    echo "Download failed."
    exit 1
fi 

tar -xzf $1

if [ $? -eq 0 ]; then
    echo "Extraction completed successfully."
else
    echo "Extraction failed."
    exit 1
fi 

rm $1
echo -e "${1} was extracted to $(pwd)/hhsuite/ \nPlease set
hhsuite_directory=$(pwd)/hhsuite/ in ./scripts/pmsa_gen.py"
