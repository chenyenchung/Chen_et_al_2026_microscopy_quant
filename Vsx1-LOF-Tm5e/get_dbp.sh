#!/bin/bash
set -e

for i in 1 2 3; do
  echo "Processing mCherry_${i} for batch 1..."
  micromamba run -n tiff python ../get_positive.py \
    -i batch1/mCherryi_${i}_D_Ey_Kn.tif \
    -m batch1/mCherryi_${i}_filtered.tif \
    -t 1000 625 -c 1 2 \
    -o batch1/mCherryi_${i}_eykn.tif
done
for i in 1 2 3 4; do
  echo "Processing VKD_${i} for batch 1..."
  micromamba run -n tiff python ../get_positive.py \
    -i batch1/VKD_${i}_D_Ey_Kn.tif \
    -m batch1/VKD_${i}_filtered.tif \
    -t 1000 625 -c 1 2 \
    -o batch1/VKD_${i}_eykn.tif
done
for i in 1 2 3; do
  echo "Processing mCherry_${i} for batch 2..."
  micromamba run -n tiff python ../get_positive.py \
    -i batch2/mCherryi_${i}_D_Ey_Kn.tif \
    -m batch2/mCherryi_${i}_filtered.tif \
    -t 1000 625 -c 1 2 \
    -o batch2/mCherryi_${i}_eykn.tif
done
for i in 1 2 3; do
  echo "Processing VKD_${i} for batch 2..."
  micromamba run -n tiff python ../get_positive.py \
    -i batch2/VKD_${i}_D_Ey_Kn.tif \
    -m batch2/VKD_${i}_filtered.tif \
    -t 1000 625 -c 1 2 \
    -o batch2/VKD_${i}_eykn.tif
done
