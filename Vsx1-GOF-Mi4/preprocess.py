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

threshold = 22500
sample_metadata = OrderedDict([
  ('LacZ_1', {
    'rotation': 85,
    'roi_path': 'data/LacZ_1_run_0.241_cp_masks.tif',
    'signal_path': 'data/LacZ_1_run_ct.tif',
    'condition': 'LacZ OE'
  }),
  ('LacZ_2', {
    'rotation': -90,
    'roi_path': 'data/LacZ_2_run_0.24_cp_masks.tif',
    'signal_path': 'data/LacZ_2_run_ct.tif',
    'condition': 'LacZ OE'
  }),
  ('LacZ_3', {
    'rotation': 85,
    'roi_path': 'data/LacZ_3_run_0.241_cp_masks.tif',
    'signal_path': 'data/LacZ_3_run_ct.tif',
    'condition': 'LacZ OE'
  }),
  ('LacZ_4', {
    'rotation': -95,
    'roi_path': 'data/LacZ_4_run_0.241_cp_masks.tif',
    'signal_path': 'data/LacZ_4_run_ct.tif',
    'condition': 'LacZ OE'
  }),
  ('Vsx1_1', {
    'rotation': 68,
    'roi_path': 'data/Vsx1_1_run_0.241_cp_masks.tif',
    'signal_path': 'data/Vsx1_1_run_ct.tif',
    'condition': 'Vsx1 OE'
  }),
  ('Vsx1_2', {
    'rotation': 120,
    'roi_path': 'data/Vsx1_2_run_0.213_cp_masks.tif',
    'signal_path': 'data/Vsx1_2_run_ct.tif',
    'condition': 'Vsx1 OE'
  }),
  ('Vsx1_3', {
    'rotation': 28,
    'roi_path': 'data/Vsx1_3_run_0.24_cp_masks.tif',
    'signal_path': 'data/Vsx1_3_run_ct.tif',
    'condition': 'Vsx1 OE'
  }),
  ('Vsx1_4', {
    'rotation': 95,
    'roi_path': 'data/Vsx1_4_run_0.24_cp_masks.tif',
    'signal_path': 'data/Vsx1_4_run_ct.tif',
    'condition': 'Vsx1 OE'
  }),
])

depth_of_interest = {'min': [], 'max': []}

for m in sample_metadata:
    sname = sample_metadata[m]['signal_path']
    roi = tifffile.imread(sample_metadata[m]['roi_path'])
    signal = tifffile.imread(sname)

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

    with_ct = pcaobj[outdf['int_1'] > threshold, 2]

    depth_ptp = np.ptp(pcaobj[:, 2])

    if sample_metadata[m]['condition'] == 'LacZ OE':
      # Keep the trimmed center of Ct range to quantify only the cells that are
      # more mature
      with_ct_min = np.percentile(with_ct, 10)
      with_ct_max = np.percentile(with_ct, 90)
      inplane_range = (
         (with_ct_min - pcaobj[:, 2].min()) / depth_ptp,
         (with_ct_max - pcaobj[:, 2].min()) / depth_ptp
      )
      depth_of_interest['min'].append(inplane_range[0])
      depth_of_interest['max'].append(inplane_range[1])
    else:
      # For Vsx1 OE, Ct is almost completely gone, so we use the mean of LacZ
      # OE range to subsample the quantification to make it somewhat comparable.
      # This is not ideal, but we don't have a better option given the data.
      mean_min = np.mean(depth_of_interest['min'])
      mean_max = np.mean(depth_of_interest['max'])
      with_ct_max = mean_max * depth_ptp + pcaobj[:, 2].min()
      with_ct_min = mean_min * depth_ptp + pcaobj[:, 2].min()
    
    depth_to_keep = [pcaobj[i, 2] < with_ct_max and pcaobj[i, 2] > with_ct_min for i in range(pcaobj.shape[0])]
    outdf['in_plane'] = depth_to_keep


    pcamat2 = standardize(rotate(pcaobj, sample_metadata[m]['rotation']))
    outdf['sample'] = sname[5:11]
    outdf['type'] = sample_metadata[m]['condition']
    outdf['x_std'] = pcamat2[:, 0]
    outdf['y_std'] = pcamat2[:, 1]
    outdf['z_std'] = pcamat2[:, 2]
    outdf['rotation'] = sample_metadata[m]['rotation']
    outdf['int_threshold'] = threshold

    in_plane_rows = int(np.sum(outdf['in_plane']))
    logger.info(
      "Sample %s: output_rows=%d, with_ct_min=%.6f, with_ct_max=%.6f, in_plane_rows=%d",
      sname[5:11],
      len(outdf),
      with_ct_min,
      with_ct_max,
      in_plane_rows,
    )
    
    if not os.path.isdir('result'):
      os.mkdir('result')
      
    outdf.to_csv(f'result/{sname[5:11]}_preprocessed.csv', index=False)