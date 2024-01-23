#!/bin/bash
##create readcounts file
> read_counts.txt
##define an array to skip
skip_values=(6 7 9 10 18)
##delcare and array to store counts
declare -a counts
##loop through 1-72
for i in {1..72}; do
##checki fi cur index is in the list to skip
  if [[ " ${skip_values[@]} " =~ " ${i} " ]]; then
    continue
  fi
  ##calculate the read count for each sample and add it to the readcounts file
  ##read count is calculated as (number of lines in the FASTQ file / 4) * 2
  ##each read in a FASTQ file spans four lines, and each sample has two files (hence the multiplication by 2)
  count=$(($(wc -l < S${i}_ILL.1.fq) / 4 * 2))
  echo "S${i}_ILL.1.fq $count" >> read_counts.txt
  counts+=($count)
done

sum=0
##calculate sum of all read counts
for value in "${counts[@]}"; do
  sum=$(($sum + $value))
done
##calculate mean
mean=$(echo "scale=2; $sum / ${#counts[@]}" | bc)

##min and max calculations
min=${counts[0]}
max=${counts[0]}

for value in "${counts[@]:1}"; do
  if (( value < min )); then
    min=$value
  fi
  if (( value > max )); then
    max=$value
  fi
done
##write final outputs to readcounts
echo "Mean value: $mean" >> read_counts.txt
echo "Minimum value: $min" >> read_counts.txt
echo "Maximum value: $max" >> read_counts.txt
