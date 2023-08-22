#!/bin/bash

output_file="combined_output.csv"

echo "SampleName,Position,PPE,Gene,Coverages" > $output_file

for file in *outcov.txt; do
    prefix=$(basename "$file" _outcov.txt)

    awk -v prefix="$prefix" '{OFS=","; print prefix, $0}' "$file" >> $output_file
done

output_file="combined_output.csv"

echo "SampleName,Position,PPE,Gene,Coverages" > $output_file

for file in *outcov.txt; do
    prefix=$(basename "$file" _outcov.txt)

    awk -v prefix="$prefix"_ONT 'BEGIN {FS="\t"; OFS=","} {print prefix, $1, $2, $3, $4}' "$file" >> $output_file
done
