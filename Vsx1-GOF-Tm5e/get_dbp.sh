#!/bin/bash
set -e

EY_THRES=800
KN_THRES=800

for i in 1 2; do
  echo "Processing LacZ_${i} for batch 1..."
  micromamba run -n tiff python ../get_positive.py \
    -i batch1/LacZ_${i}_D_Ey_Kn.tif \
    -m batch1/LacZ_${i}_filtered.tif \
    -t ${EY_THRES} ${KN_THRES} -c 1 2 \
    -o batch1/LacZ_${i}_eykn.tif
done
for i in 1 2 3; do
  echo "Processing Vsx1_${i} for batch 1..."
  micromamba run -n tiff python ../get_positive.py \
    -i batch1/Vsx1_${i}_D_Ey_Kn.tif \
    -m batch1/Vsx1_${i}_filtered.tif \
    -t ${EY_THRES} ${KN_THRES} -c 1 2 \
    -o batch1/Vsx1_${i}_eykn.tif
done

for i in 1 2 3 4; do
  echo "Processing LacZ_${i} for batch 2..."
  micromamba run -n tiff python ../get_positive.py \
    -i batch2/LacZ_${i}_D_Ey_Kn.tif \
    -m batch2/LacZ_${i}_filtered.tif \
    -t ${EY_THRES} ${KN_THRES} -c 1 2 \
    -o batch2/LacZ_${i}_eykn.tif
done

for i in 1 2 3; do
  echo "Processing Vsx1_${i} for batch 2..."
  micromamba run -n tiff python ../get_positive.py \
    -i batch2/Vsx1_${i}_D_Ey_Kn.tif \
    -m batch2/Vsx1_${i}_filtered.tif \
    -t ${EY_THRES} ${KN_THRES} -c 1 2 \
    -o batch3/Vsx1_${i}_eykn.tif
done
