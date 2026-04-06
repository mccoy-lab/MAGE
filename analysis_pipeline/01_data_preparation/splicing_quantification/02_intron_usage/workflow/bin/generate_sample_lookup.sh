#!/usr/bin/env bash

LEAFCUTTER_DIR=$1
METADATA_FILE=$2

for FILE in ${LEAFCUTTER_DIR}/*.sorted.gz; do
    INTERNAL_ID=$(basename $FILE .sorted.gz | cut -d "." -f 1) # Extract internal library ID from filename 
    TGP_ID=$(grep $INTERNAL_ID $METADATA_FILE | awk -F"\t" '{print $5}') # Extract 1000 Genomes ID from metadata file using internal library ID
    printf "${INTERNAL_ID}\t${TGP_ID}\n" >> ${LEAFCUTTER_DIR}/sample_lookup.txt # Append internal library ID and 1000 Genomes ID to sample lookup file
done
