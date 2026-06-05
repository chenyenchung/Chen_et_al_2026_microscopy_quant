#!/bin/bash
set -e

BSH_THRES=800

# for i in 1 2 3 4; do
#   echo "Processing LacZ_${i} for batch 1..."
#   micromamba run -n tiff python ../get_positive.py \
#     -i batch1/LacZ_${i}_Dimm_Bsh_Fs.tif \
#     -m batch1/LacZ_${i}_filtered.tif \
#     -t ${BSH_THRES} -c 1 \
#     -o batch1/LacZ_${i}_bsh.tif
# done
# for i in 1 2 3; do
#   echo "Processing Bi_${i} for batch 1..."
#   micromamba run -n tiff python ../get_positive.py \
#     -i batch1/Bi_${i}_Dimm_Bsh_Fs.tif \
#     -m batch1/Bi_${i}_filtered.tif \
#     -t ${BSH_THRES} -c 1 \
#     -o batch1/Bi_${i}_bsh.tif
# done
# 
# for i in 1 2 3; do
#   echo "Processing LacZ_${i} for batch 2..."
#   micromamba run -n tiff python ../get_positive.py \
#     -i batch2/LacZ_${i}_Dimm_Bsh_Fs.tif \
#     -m batch2/LacZ_${i}_filtered.tif \
#     -t ${BSH_THRES} -c 1 \
#     -o batch2/LacZ_${i}_bsh.tif
# done
# for i in 1 2 3; do
#   echo "Processing Bi_${i} for batch 2..."
#   micromamba run -n tiff python ../get_positive.py \
#     -i batch2/Bi_${i}_Dimm_Bsh_Fs.tif \
#     -m batch2/Bi_${i}_filtered.tif \
#     -t ${BSH_THRES} -c 1 \
#     -o batch2/Bi_${i}_bsh.tif
# done

for i in 2 3; do
  echo "Processing LacZ_${i} for batch 3..."
  micromamba run -n tiff python ../get_positive.py \
    -i batch3/LacZ_${i}_Dimm_Bsh_Fs.tif \
    -m batch3/LacZ_${i}_filtered.tif \
    -t ${BSH_THRES} -c 1 \
    -o batch3/LacZ_${i}_bsh.tif
done
for i in 1 2 3 4; do
  echo "Processing Bi_${i} for batch 3..."
  micromamba run -n tiff python ../get_positive.py \
    -i batch3/Bi_${i}_Dimm_Bsh_Fs.tif \
    -m batch3/Bi_${i}_filtered.tif \
    -t ${BSH_THRES} -c 1 \
    -o batch3/Bi_${i}_bsh.tif
done
