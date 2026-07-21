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

sample_metadata = OrderedDict([
   ('mCherryi_1', {
      'condition': 'mCherry RNAi',
      'roi_path': './batch1-pxb-RNAi/mCherryi_1_0.24_cp_masks.tif',
      'bsh_roi_path': './batch1-pxb-RNAi-Bsh/mCherryi_1_0.24_cp_masks.tif',
      'signal_path': './batch1-pxb-RNAi/mCherryi_1_0.24.tif',
      'rotation': 110
   }),
   ('mCherryi_2', {
      'condition': 'mCherry RNAi',
      'roi_path': './batch1-pxb-RNAi/mCherryi_2_0.24_cp_masks.tif',
      'bsh_roi_path': './batch1-pxb-RNAi-Bsh/mCherryi_2_0.24_cp_masks.tif',
      'signal_path': './batch1-pxb-RNAi/mCherryi_2_0.24.tif',
      'rotation': -90
   }),
   ('mCherryi_3', {
     'condition': 'mCherry RNAi',
     'roi_path': './batch1-pxb-RNAi/mCherryi_3_0.24_cp_masks.tif',
     'bsh_roi_path': './batch1-pxb-RNAi-Bsh/mCherryi_3_0.24_cp_masks.tif',
     'signal_path': './batch1-pxb-RNAi/mCherryi_3_0.24.tif',
     'rotation': -75
   }),
   ('Vsx1i_1', {
      'condition': 'Vsx1 NE RNAi',
      'roi_path': './batch1-pxb-RNAi/Vsx1i_1_0.24_cp_masks.tif',
      'bsh_roi_path': './batch1-pxb-RNAi-Bsh/Vsx1i_1_0.24_cp_masks.tif',
      'signal_path': './batch1-pxb-RNAi/Vsx1i_1_0.24.tif',
      'rotation': 75
   }),
   ('Vsx1i_2', {
      'condition': 'Vsx1 NE RNAi',
      'roi_path': './batch1-pxb-RNAi/Vsx1i_2_0.24_cp_masks.tif',
      'bsh_roi_path': './batch1-pxb-RNAi-Bsh/Vsx1i_2_0.24_cp_masks.tif',
      'signal_path': './batch1-pxb-RNAi/Vsx1i_2_0.24.tif',
      'rotation': 90
   }),
   ('Vsx1i_3', {
      'condition': 'Vsx1 NE RNAi',
      'roi_path': './batch1-pxb-RNAi/Vsx1i_3_0.24_cp_masks.tif',
      'bsh_roi_path': './batch1-pxb-RNAi-Bsh/Vsx1i_3_0.24_cp_masks.tif',
      'signal_path': './batch1-pxb-RNAi/Vsx1i_3_0.24.tif',
      'rotation': -75
   }),
])

for m in sample_metadata:
  logger.info("Processing sample %s", m)
  roi = tifffile.imread(sample_metadata[m]['roi_path'])
  bsh_roi = tifffile.imread(sample_metadata[m]['bsh_roi_path'])
  signal = tifffile.imread(sample_metadata[m]['signal_path'])

  bsh_num = np.count_nonzero(np.unique(bsh_roi))

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


  pcamat2 = standardize(rotate(pcaobj, sample_metadata[m]['rotation']))
  outdf['sample'] = m
  outdf['type'] = sample_metadata[m]['condition']
  outdf['bsh_num'] = bsh_num
  outdf['x_std'] = pcamat2[:, 0]
  outdf['y_std'] = pcamat2[:, 1]
  outdf['z_std'] = pcamat2[:, 2]
  outdf['r'] = np.hypot(outdf['x_std'], outdf['y_std'])
  outdf['theta'] = np.arctan2(outdf['y_std'], outdf['x_std'])
  outdf['rotation'] = sample_metadata[m]['rotation']

  logger.info(
    "Sample %s: output_rows=%d, bsh_num=%d",
    m,
    len(outdf),
    bsh_num,
  )
  
  if not os.path.isdir('result'):
    os.mkdir('result')
    
  outdf.to_csv(f'result/{m}_preprocessed.csv', index=False)
