#!/usr/bin/env python
import argparse
import logging
import warnings
from pathlib import Path

import numpy as np
import tifffile
from skimage.measure import regionprops_table

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


def validate_channels(img, channels):
    if img.ndim == 3:
        invalid = [channel for channel in channels if channel != 0]
        if invalid:
            raise ValueError("Single-channel images only support channel 0")
        return

    if img.ndim == 4:
        n_channels = img.shape[3]
        invalid = [channel for channel in channels if channel < 0 or channel >= n_channels]
        if invalid:
            raise ValueError(f"Channel(s) {invalid} are outside the image channel range 0-{n_channels - 1}")
        return

    raise ValueError(f"Expected a 3D or 4D image after channel-axis normalization, got shape {img.shape}")


def intensity_means_for_channels(stats, img, channels):
    if img.ndim == 3:
        return np.column_stack([stats["intensity_mean"] for _ in channels])

    return np.column_stack([stats[f"intensity_mean-{channel}"] for channel in channels])


image_path = Path(args.image)
mask_path = Path(args.mask)

if image_path.is_file():
    img = tifffile.imread(image_path)
    if img.ndim > 3:
        img = np.moveaxis(img, 1, -1)
else:
    raise FileNotFoundError(f"{image_path} is not found")

logging.info(f"Image shape: {img.shape}")
validate_channels(img, args.channels)

if mask_path.is_file():
    roi = tifffile.imread(mask_path)
else:
    raise FileNotFoundError(f"{mask_path} is not found")

stats = regionprops_table(roi, img, properties=("label", "intensity_mean"))
labels = stats["label"].astype(np.intp, copy=False)
means = intensity_means_for_channels(stats, img, args.channels)
thresholds = np.asarray(args.thres, dtype=means.dtype)
to_discard = np.any(means < thresholds, axis=1)
discard_labels = labels[to_discard]

if discard_labels.size:
    remove_label = np.zeros(int(roi.max()) + 1, dtype=bool)
    remove_label[discard_labels] = True
    roi[remove_label[roi]] = 0

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

    if fn_ori != fn:
        warnings.warn(f"File {fn_ori} already exists, saving to {fn}")
    tifffile.imwrite(fn, roi)
