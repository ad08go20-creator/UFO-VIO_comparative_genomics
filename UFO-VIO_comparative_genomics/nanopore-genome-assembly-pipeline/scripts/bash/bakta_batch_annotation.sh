#!/bin/bash
# Batch-run bakta genome annotation over every .fna file in a directory.

set -euo pipefail

outputdir="bakta"
fnadir="/path/to/fna/files"
extension=".fna"
cpus=4
dbdir="/path/to/bakta/db"

mkdir -p "$outputdir"

for file in "${fnadir}"/*"${extension}"; do
    bname=$(basename "$file" "${extension}")
    bakta --db "$dbdir" --prefix "$bname" --output "$outputdir/$bname" \
        --threads "$cpus" --verbose "$file"
done