#!/usr/bin/env python
import tifffile
import pandas as pd
import numpy as np
from skimage.measure import regionprops
import warnings
import logging
import sys
import matplotlib.pyplot as plt
import scipy.stats as scstats
from sklearn.decomposition import PCA
import os
from collections import OrderedDict

logging.basicConfig(
  level=logging.INFO,
  format='%(message)s',
  handlers=[logging.StreamHandler(sys.stdout)]
)
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

sample_metadata = OrderedDict([
  ('mChKD_1', {
    'rotation': -95,
    'thres': 25000,
    'roi_path': 'data/mChi_1_run_0.224_cp_masks.tif',
    'signal_path': 'data/mChi_1_ct_run.tif',
    'condition': 'mCherry RNAi'
  }),
  ('mChKD_2', {
    'rotation': 95,
    'thres': 25000,
    'roi_path': 'data/mChi_2_run_0.241_cp_masks.tif',
    'signal_path': 'data/mChi_2_ct_run.tif',
    'condition': 'mCherry RNAi'
  }),
  ('mChKD_3', {
    'rotation': -90,
    'thres': 25000,
    'roi_path': 'data/mChi_3_run_0.241_cp_masks.tif',
    'signal_path': 'data/mChi_3_ct_run.tif',
    'condition': 'mCherry RNAi'
  }),
  ('VsxDKD_1', {
    'rotation': -105,
    'thres': 25000,
    'roi_path': 'data/VKD_1_run_0.241_cp_masks.tif',
    'signal_path': 'data/VKD_1_ct_run.tif',
    'condition': 'Vsx1/2 RNAi'
  }),
  ('VsxDKD_2', {
    'rotation': -100,
    'thres': 25000,
    'roi_path': 'data/VKD_2_run_0.241_cp_masks.tif',
    'signal_path': 'data/VKD_2_ct_run.tif',
    'condition': 'Vsx1/2 RNAi'
  }),
  ('VsxDKD_3', {
    'rotation': -85,
    'thres': 25000,
    'roi_path': 'data/VKD_3_run_0.241_cp_masks.tif',
    'signal_path': 'data/VKD_3_ct_run.tif',
    'condition': 'Vsx1/2 RNAi'
  }),
])

for m in sample_metadata:
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

  # Rotate axes to better align samples
  pca = PCA(n_components=3, svd_solver='auto')
  pcaobj = pca.fit_transform(
      outdf[['centroid_x', 'centroid_y', 'centroid_z']].to_numpy()
  )

  with_ct = pcaobj[outdf['int_0'] > sample_metadata[m]['thres'], 2]
  with_ct_min = np.percentile(with_ct, 10)
  with_ct_max = np.percentile(with_ct, 90)
  depth_to_keep = [pcaobj[i, 2] < with_ct_max and pcaobj[i, 2] > with_ct_min for i in range(pcaobj.shape[0])]

  pcamat2 = standardize(rotate(pcaobj, sample_metadata[m]['rotation']))
  outdf['in_plane'] = depth_to_keep
  outdf['sample'] = m
  outdf['type'] = sample_metadata[m]['condition']
  outdf['x_std'] = pcamat2[:, 0]
  outdf['y_std'] = pcamat2[:, 1]
  outdf['z_std'] = pcamat2[:, 2]
  outdf['rotation'] = sample_metadata[m]['rotation']
  outdf['int_threshold'] = sample_metadata[m]['thres']

  in_plane_rows = int(np.sum(outdf['in_plane']))
  logger.info(
    "Sample %s: output_rows=%d, with_ct_min=%.6f, with_ct_max=%.6f, in_plane_rows=%d",
    m,
    len(outdf),
    with_ct_min,
    with_ct_max,
    in_plane_rows,
  )
  
  if not os.path.isdir('result'):
    os.mkdir('result')
    
  outdf.to_csv(f'result/{m}_preprocessed.csv', index=False)