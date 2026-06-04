#!/usr/bin/env python
import tifffile
import pandas as pd
import numpy as np
from skimage.measure import regionprops
import warnings
import argparse 
from pathlib import Path
from sklearn.decomposition import PCA
import os
from collections import OrderedDict
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def rotate(mat, angle_deg, center=None, dims=(0, 1)):
    """
    Rotate points in the plane defined by dims by angle_deg degrees.

    Parameters
    ----------
    mat : (n, d) array
        Input coordinates, e.g. PCA coordinates.
    angle_deg : float
        Rotation angle in degrees. Positive = counterclockwise.
    center : array-like of shape (2,), optional
        Center of rotation in the selected dims.
        If None, uses the centroid in those dims.
    dims : tuple of int
        Which two columns to rotate. Default: (0, 1).

    Returns
    -------
    out : (n, d) array
        Rotated copy of mat.
    """
    mat = np.asarray(mat, dtype=float)
    out = mat.copy()

    Y = out[:, dims]

    if center is None:
        center = Y.mean(axis=0)
    else:
        center = np.asarray(center, dtype=float)

    theta = np.deg2rad(angle_deg)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s],
                  [s,  c]])

    out[:, dims] = (Y - center) @ R.T + center
    return out


def standardize(mat, dims=(0, 1)):
    """
    Scale each selected axis independently to span [-1, 1].

    Parameters
    ----------
    mat : (n, d) array
    dims : tuple of int
        Which columns to rescale.

    Returns
    -------
    out : (n, d) array
    """
    mat = np.asarray(mat, dtype=float)
    out = mat.copy()

    for d in dims:
        x = out[:, d]
        xmin, xmax = x.min(), x.max()
        if np.isclose(xmax, xmin):
            out[:, d] = 0.0
        else:
            out[:, d] = 2 * (x - xmin) / (xmax - xmin) - 1

    return out

threshold = 800
control_condition = 'LacZ OE'
sample_metadata = OrderedDict([
    ('LacZ_1_b1', {
       'rotation': 95,
       'condition': 'LacZ OE',
       'roi_path': 'batch1/LacZ_1_D_Ey_Kn_cp_masks.tif',
       'signal_path': 'batch1/LacZ_1_D_Ey_Kn.tif'
    }),
    ('LacZ_2_b1', {
       'rotation': -115,
       'condition': 'LacZ OE',
       'roi_path': 'batch1/LacZ_2_D_Ey_Kn_cp_masks.tif',
       'signal_path': 'batch1/LacZ_2_D_Ey_Kn.tif'
    }),
    ('Vsx1_1_b1', {
       'rotation': -95,
       'condition': 'Vsx1 GOF',
       'roi_path': 'batch1/Vsx1_1_D_Ey_Kn_cp_masks.tif',
       'signal_path': 'batch1/Vsx1_1_D_Ey_Kn.tif'
    }),
    ('Vsx1_2_b1', {
       'rotation': -90,
       'condition': 'Vsx1 GOF',
       'roi_path': 'batch1/Vsx1_2_D_Ey_Kn_cp_masks.tif',
       'signal_path': 'batch1/Vsx1_2_D_Ey_Kn.tif'
    }),
    ('Vsx1_3_b1', {
       'rotation': -90,
       'condition': 'Vsx1 GOF',
       'roi_path': 'batch1/Vsx1_3_D_Ey_Kn_cp_masks.tif',
       'signal_path': 'batch1/Vsx1_3_D_Ey_Kn.tif'
    }),
    ('LacZ_1_b2', {
       'rotation': 90,
       'condition': 'LacZ OE',
       'roi_path': 'batch2/LacZ_1_D_Ey_Kn_cp_masks.tif',
       'signal_path': 'batch2/LacZ_1_D_Ey_Kn.tif'
    }),
    ('LacZ_2_b2', {
       'rotation': -95,
       'condition': 'LacZ OE',
       'roi_path': 'batch2/LacZ_2_D_Ey_Kn_cp_masks.tif',
       'signal_path': 'batch2/LacZ_2_D_Ey_Kn.tif'
    }),
    ('LacZ_3_b2', {
       'rotation': -90,
       'condition': 'LacZ OE',
       'roi_path': 'batch2/LacZ_3_D_Ey_Kn_cp_masks.tif',
       'signal_path': 'batch2/LacZ_3_D_Ey_Kn.tif'
    }),
    ('LacZ_4_b2', {
        'rotation': 90,
        'condition': 'LacZ OE',
        'roi_path': 'batch2/LacZ_4_D_Ey_Kn_cp_masks.tif',
        'signal_path': 'batch2/LacZ_4_D_Ey_Kn.tif'
    }),
    ('Vsx1_1_b2', {
        'rotation': -105,
        'condition': 'Vsx1 GOF',
        'roi_path': 'batch2/Vsx1_1_D_Ey_Kn_cp_masks.tif',
        'signal_path': 'batch2/Vsx1_1_D_Ey_Kn.tif'
    }),
    ('Vsx1_2_b2', {
        'rotation': 90,
        'condition': 'Vsx1 GOF',
        'roi_path': 'batch2/Vsx1_2_D_Ey_Kn_cp_masks.tif',
        'signal_path': 'batch2/Vsx1_2_D_Ey_Kn.tif'
    }),
    ('Vsx1_3_b2', {
        'rotation': -100,
        'condition': 'Vsx1 GOF',
        'roi_path': 'batch2/Vsx1_3_D_Ey_Kn_cp_masks.tif',
        'signal_path': 'batch2/Vsx1_3_D_Ey_Kn.tif'
    })
])

sample_results = OrderedDict()
control_depth_ranges = []

for m in sample_metadata:
  logger.info("Processing sample %s", m)
  roi = tifffile.imread(sample_metadata[m]['roi_path'])
  signal = tifffile.imread(sample_metadata[m]['signal_path'])

  if signal.ndim == 4:
    # Fiji saves channel as the second dimension, but scikit-image expects it
    # as the last dimension. Reorder axes accordingly.
    signal = np.moveaxis(signal, 1, -1)

  stats = regionprops(roi, signal)

  out = {
      'label': [], 'centroid_x': [], 'centroid_y': [], 'centroid_z': [],
      'area': [], 'aspherity': []
  }

  if signal.ndim == 4:
    for i in range(signal.shape[3]):
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
        
    if signal.ndim == 3:
        out['int'].append(i.intensity_mean)
    else:
        for j in range(signal.shape[3]):
            out[f'int_{j}'].append(i.intensity_mean[j])

    out['aspherity'].append(i.axis_major_length/i.equivalent_diameter_area)

  outdf = pd.DataFrame(out)

  pca = PCA(n_components=3, svd_solver='auto')
  pcaobj = pca.fit_transform(
      outdf[['centroid_x', 'centroid_y', 'centroid_z']].to_numpy()
  )

  # D appears to be always on in shallow/nascent cells, so we slice on PC3
  # to only consider cells that have turned D off in the Vsx domain in
  # control samples.
  depth_ptp = np.ptp(pcaobj[:, 2])
  if np.isclose(depth_ptp, 0):
    raise ValueError(f"Sample {m} has no PC3 depth spread")

  if sample_metadata[m]['condition'] == control_condition:
    with_Tm5e = pcaobj[outdf['int_0'] < threshold, 2]
    if len(with_Tm5e) == 0:
      raise ValueError(
        f"Control sample {m} has no int_0 values below threshold {threshold}"
      )
    with_Tm5e_min = np.percentile(with_Tm5e, 10)
    with_Tm5e_max = np.percentile(with_Tm5e, 90)
    control_depth_ranges.append((
         (with_Tm5e_min - pcaobj[:, 2].min()) / depth_ptp,
         (with_Tm5e_max - pcaobj[:, 2].min()) / depth_ptp
    ))

  # Rotate axes to better align samples
  pcamat2 = standardize(rotate(pcaobj, sample_metadata[m]['rotation']))
  sample_results[m] = {
    'outdf': outdf,
    'pcaobj': pcaobj,
    'pcamat2': pcamat2,
    'depth_ptp': depth_ptp,
  }

if not control_depth_ranges:
  raise ValueError(f"No control samples found for condition {control_condition!r}")

mean_min, mean_max = np.mean(control_depth_ranges, axis=0)

for m, result in sample_results.items():
  outdf = result['outdf']
  pcaobj = result['pcaobj']
  pcamat2 = result['pcamat2']
  depth_ptp = result['depth_ptp']

  with_Tm5e_min = mean_min * depth_ptp + pcaobj[:, 2].min()
  with_Tm5e_max = mean_max * depth_ptp + pcaobj[:, 2].min()
  depth_to_keep = (pcaobj[:, 2] > with_Tm5e_min) & (pcaobj[:, 2] < with_Tm5e_max)

  outdf['in_plane'] = depth_to_keep
  outdf['sample'] = m
  outdf['type'] = sample_metadata[m]['condition']
  outdf['x_std'] = pcamat2[:, 0]
  outdf['y_std'] = pcamat2[:, 1]
  outdf['z_std'] = pcamat2[:, 2]
  outdf['rotation'] = sample_metadata[m]['rotation']
  outdf['int_threshold'] = threshold

  in_plane_rows = int(np.sum(outdf['in_plane']))
  logger.info(
    "Sample %s: output_rows=%d, with_Tm5e_min=%.6f, with_Tm5e_max=%.6f, in_plane_rows=%d",
    m,
    len(outdf),
    with_Tm5e_min,
    with_Tm5e_max,
    in_plane_rows,
  )
  
  if not os.path.isdir('result'):
    os.mkdir('result')
    
  outdf.to_csv(f'result/{m}_preprocessed.csv', index=False)
