#!/bin/bash

> read_counts.txt

skip_values=(6 7 9 10 18)

declare -a counts

for i in {1..72}; do
  if [[ " ${skip_values[@]} " =~ " ${i} " ]]; then
    continue
  fi
  count=$(($(wc -l < S${i}_ILL.1.fq) / 4 * 2))
  echo "S${i}_ILL.1.fq $count" >> read_counts.txt
  counts+=($count)
done

sum=0
for value in "${counts[@]}"; do
  sum=$(($sum + $value))
done

mean=$(echo "scale=2; $sum / ${#counts[@]}" | bc)


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

echo "Mean value: $mean" >> read_counts.txt
echo "Minimum value: $min" >> read_counts.txt
echo "Maximum value: $max" >> read_counts.txt
