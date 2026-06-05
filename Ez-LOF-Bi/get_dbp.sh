#!/bin/bash
set -e

BI_THRES=1.25
for i in 4; do
  echo "Processing Ezi_${i} for batch 1..."
  micromamba run -n tiff python ../get_positive.py \
    -i batch1/Ezi_${i}_highlaser_bi.tif \
    -m batch1/Ezi_${i}_highlaser_filtered.tif \
    -t ${BI_THRES} -c 0 --znorm \
    -o batch1/Ezi_${i}_bi_mask.tif
done
#for i in 2 3 4; do
#  echo "Processing Ezi_${i} for batch 1..."
#  micromamba run -n tiff python ../get_positive.py \
#    -i batch1/Ezi_${i}_highlaser_bi.tif \
#    -m batch1/Ezi_${i}_highlaser_filtered.tif \
#    -t ${BI_THRES} -c 0 --znorm \
#    -o batch1/Ezi_${i}_bi_mask.tif
#done

# for i in 1 2 3 4; do
#   echo "Processing mCherryi_${i} for batch 2..."
#   micromamba run -n tiff python ../get_positive.py \
#     -i batch2/mCherryi_${i}_bi.tif \
#     -m batch2/mCherryi_${i}_filtered.tif \
#     -t ${BI_THRES} -c 0 --znorm \
#     -o batch2/mCherryi_${i}_bi_mask.tif
# done
# for i in 1 2; do
#   echo "Processing Ezi_${i} for batch 2..."
#   micromamba run -n tiff python ../get_positive.py \
#     -i batch2/Ezi_${i}_bi.tif \
#     -m batch2/Ezi_${i}_filtered.tif \
#     -t ${BI_THRES} -c 0 --znorm \
#     -o batch2/Ezi_${i}_bi_mask.tif
# done
