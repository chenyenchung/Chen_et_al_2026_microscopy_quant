#!/usr/bin/env python
import tifffile
import pandas as pd
import numpy as np
from skimage.measure import regionprops
import warnings
import argparse 
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--image', type=str, required=True)
parser.add_argument('-m', '--mask', type=str, required=True)
parser.add_argument('-t', '--thres', type=float, required=True, nargs='+')
parser.add_argument('-c', '--channels', type=int, required=True, nargs='+')
parser.add_argument('-o', '--output', type=str)
parser.add_argument('-s', '--suffix', type=str, default='filtered')

args = parser.parse_args()

if len(args.channels) != len(args.thres):
    raise ValueError("The number of channels and thresholds must be the same")

image_path = Path(args.image)
mask_path = Path(args.mask)

if image_path.is_file():
    img = tifffile.imread(image_path)
    if img.ndim > 3:
      img = np.moveaxis(img, 1, -1)
else:
    raise FileNotFoundError(f"{image_path} is not found")

logging.info(f"Image shape: {img.shape}")

if mask_path.is_file():
    roi = tifffile.imread(mask_path)
else:
    raise FileNotFoundError(f"{mask_path} is not found")

stats = regionprops(roi, img)

out = {'label': [], 'centroid_x': [], 'centroid_y': [], 'centroid_z': [], 'area': [], 'aspherity': []}

if img.ndim == 4:
    for i in range(img.shape[3]):
        out[f'int_{i}'] = []
else:
    out['int'] = []

for i in stats:
    out['label'].append(i.label)
    out['centroid_x'].append(i.centroid[0])
    out['centroid_y'].append(i.centroid[1])
    out['centroid_z'].append(i.centroid[2])

    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      out['area'].append(i.area_convex)
        
    if img.ndim == 3:
      out['int'].append(i.intensity_mean)
    else:
      for j in range(img.shape[3]):
        out[f'int_{j}'].append(i.intensity_mean[j])

    out['aspherity'].append(i.axis_major_length/i.equivalent_diameter_area)

outdf = pd.DataFrame(out)

if len(args.channels) == 1:
    discarddf = outdf[outdf[f'int_{args.channels[0]}'] < args.thres[0]] 
else:
    to_discard = np.zeros(outdf.shape[0], dtype=bool)
    for i in range(len(args.channels)):
        to_discard = to_discard | (outdf[f'int_{args.channels[i]}'] < args.thres[i])
    discarddf = outdf[to_discard] 

roi[np.isin(roi, discarddf['label'])] = 0

if args.output is None:
    fn = f"{mask_path.stem}_{args.suffix}.tif"
    tifffile.imwrite(fn, roi)
else:
    fn = Path(args.output)
    fn_ori = fn
    i = 1
    while fn.is_file():
        fn = Path(f"{fn.parent}/{fn.stem}_{i}.tif")
        i += 1

    warnings.warn(f"File {fn_ori} already exists, saving to {fn}")
    tifffile.imwrite(fn, roi)
