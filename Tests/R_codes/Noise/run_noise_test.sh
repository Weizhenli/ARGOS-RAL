#!/bin/bash

polys=(3 4 5)
diffs=(3 4 5)
sds=(0.1 1 10)

for p in "${polys[@]}"; do
    for sd in "${sds[@]}"; do
            sbatch --export=ALL,d="$p",p="$p",sd="$sd" random_noise_test.sh
    done
done
