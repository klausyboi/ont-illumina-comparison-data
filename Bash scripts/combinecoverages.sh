#!/bin/bash
##define output name
output_file="combined_output.csv"
##write the header to the output file
echo "SampleName,Position,PPE,Gene,Coverages" > $output_file
##iterate over all files in the curr dir that end with outcov.txt
for file in *outcov.txt; do
    ##extract base name from the filename
    prefix=$(basename "$file" _outcov.txt)
    
    awk -v prefix="$prefix" '{OFS=","; print prefix, $0}' "$file" >> $output_file
done

output_file="combined_output.csv"

echo "SampleName,Position,PPE,Gene,Coverages" > $output_file

for file in *outcov.txt; do
    prefix=$(basename "$file" _outcov.txt)
    ##field separator as tab and changed to , print the prefix and 1-4 columns
    awk -v prefix="$prefix"_ONT 'BEGIN {FS="\t"; OFS=","} {print prefix, $1, $2, $3, $4}' "$file" >> $output_file
done
