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
from scipy.stats import gaussian_kde

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
    ('CantonS_1', {
        'rotation': -85,
        'yflip': False,
        'roi_path': 'data/CantonS_1_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/CantonS_1_RFP_Bsh_Dimm.tif',
        'condition': 'CantonS'
    }),
    ('CantonS_2', {
        'rotation': -80,
        'yflip': False,
        'roi_path': 'data/CantonS_2_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/CantonS_2_RFP_Bsh_Dimm.tif',
        'condition': 'CantonS'
    }),
    ('CantonS_3', {
        'rotation': 95,
        'yflip': True,
        'roi_path': 'data/CantonS_3_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/CantonS_3_RFP_Bsh_Dimm.tif',
        'condition': 'CantonS'
    }),
    ('CantonS_4', {
        'rotation': 100,
        'yflip': True,
        'roi_path': 'data/CantonS_4_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/CantonS_4_RFP_Bsh_Dimm.tif',
        'condition': 'CantonS'
    }),
    ('CantonS_5', {
        'rotation': 95,
        'yflip': True,
        'roi_path': 'data/CantonS_5_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/CantonS_5_RFP_Bsh_Dimm.tif',
        'condition': 'CantonS'
    }),
    ('CantonS_6', {
        'rotation': -90,
        'yflip': False,
        'roi_path': 'data/CantonS_6_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/CantonS_6_RFP_Bsh_Dimm.tif',
        'condition': 'CantonS'
    }
    ),
    ('CantonS_7', {
        'rotation': -90,
        'yflip': False,
        'roi_path': 'data/CantonS_7_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/CantonS_7_RFP_Bsh_Dimm.tif',
        'condition': 'CantonS'
    }),
    ('sgBi_1', {
        'rotation': -100,
        'yflip': True,
        'roi_path': 'data/Bi_1_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/Bi_1_RFP_Bsh_Dimm.tif',
        'condition': 'sgBi'
    }),
    ('sgBi_2', {
        'rotation': 85,
        'yflip': True,
        'roi_path': 'data/Bi_2_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/Bi_2_RFP_Bsh_Dimm.tif',
        'condition': 'sgBi'
    }),
    ('sgBi_3', {
        'rotation': -100,
        'yflip': True,
        'roi_path': 'data/Bi_3_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/Bi_3_RFP_Bsh_Dimm.tif',
        'condition': 'sgBi'
    }),
    ('sgBi_4', {
        'rotation': -90,
        'yflip': False,
        'roi_path': 'data/Bi_4_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/Bi_4_RFP_Bsh_Dimm.tif',
        'condition': 'sgBi'
    }),
    ('sgBi_5', {
        'rotation': 80,
        'yflip': True,
        'roi_path': 'data/Bi_5_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/Bi_5_RFP_Bsh_Dimm.tif',
        'condition': 'sgBi'
    }),
    ('sgBi_6', {
        'rotation': 90,
        'yflip': True,
        'roi_path': 'data/Bi_6_RFP_Bsh_Dimm_cp_masks.tif',
        'signal_path': 'data/Bi_6_RFP_Bsh_Dimm.tif',
        'condition': 'sgBi'
    }),
])

dimm_thres = 500
rfp_thres = 1000

# In the case of interactive reruns, this prevents appending to an existing
# outdf variable. In a clean run, this will do nothing.
try: del outdf
except NameError: pass
try: del depth_of_interest
except NameError: pass

# Process each sample and concatenate results into outdf
depth_of_interest = {'min': [], 'max': []}
for sample_name, meta in sample_metadata.items():
    logger.info(f'Processing {sample_name}...')
    roi = tifffile.imread(meta['roi_path'])
    signal = tifffile.imread(meta['signal_path'])

    if signal.ndim > 3:
      signal = np.moveaxis(signal, 1, -1)
    
    stats = regionprops(roi, signal)
    temp = {'label': [], 'centroid_x': [], 'centroid_y': [], 'centroid_z': [], 'area': [], 'aspherity': []}

    if signal.ndim == 4:
        for i in range(signal.shape[3]):
            temp[f'int_{i}'] = []
    else:
        temp['int'] = []


    for i in stats:
        temp['label'].append(i.label)
        temp['centroid_x'].append(i.centroid[0])
        temp['centroid_y'].append(i.centroid[1])
        temp['centroid_z'].append(i.centroid[2])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            temp['area'].append(i.area_convex)
            
        if signal.ndim == 3:
            temp['int'].append(i.intensity_mean)
        else:
            for j in range(signal.shape[3]):
                temp[f'int_{j}'].append(i.intensity_mean[j])

        temp['aspherity'].append(i.axis_major_length/i.equivalent_diameter_area)

    tempdf = pd.DataFrame(temp)
    pca = PCA(n_components=3, svd_solver='auto')
    pcaobj = pca.fit_transform(tempdf[['centroid_x', 'centroid_y', 'centroid_z']].to_numpy())

    pcamat2 = standardize(rotate(pcaobj, meta['rotation']))
    if meta['yflip']:
        pcamat2[:, 1] *= -1

    tempdf[['pca_x', 'pca_y', 'pca_z']] = pcamat2
    # tempdf = tempdf[tempdf['int_2'] > dimm_thres]
    tempdf['sample'] = sample_name
    tempdf['condition'] = meta['condition']

    with_dimm = pcaobj[tempdf['int_2'] > dimm_thres, 2]
    depth_ptp = np.ptp(pcaobj[:, 2])

    # Define the depth of interest as where Dimm+ cells are found in the
    # control samples. This is used to exclude cells that are too young and
    # haven't turned on Dimm yet.
    if meta['condition'] == 'CantonS':
      # Keep the trimmed center of Dimm range to quantify only the cells that are
      # more mature
      with_dimm_min = np.percentile(with_dimm, 10)
      with_dimm_max = np.percentile(with_dimm, 90)
      inplane_range = (
         (with_dimm_min - pcaobj[:, 2].min()) / depth_ptp,
         (with_dimm_max - pcaobj[:, 2].min()) / depth_ptp
      )
      depth_of_interest['min'].append(inplane_range[0])
      depth_of_interest['max'].append(inplane_range[1])
    else:
      # For Bi sKO, Dimm is mostly gone, so we use the mean of CantonS
      # range to subsample the quantification to make it somewhat comparable.
      # This is not ideal, but we don't have a better option given the data.
      mean_min = np.mean(depth_of_interest['min'])
      mean_max = np.mean(depth_of_interest['max'])
      with_dimm_max = mean_max * depth_ptp + pcaobj[:, 2].min()
      with_dimm_min = mean_min * depth_ptp + pcaobj[:, 2].min()

    depth_to_keep = [pcaobj[i, 2] < with_dimm_max and pcaobj[i, 2] > with_dimm_min for i in range(pcaobj.shape[0])]
    tempdf['in_plane'] = depth_to_keep
    
    try:
        outdf = pd.concat([outdf, tempdf], ignore_index=True)
    except NameError:
        outdf = tempdf.copy()

    del tempdf

# Take only CantonS (control) TE cells to define an ROI where TE cells are
# expected.
TE_df = outdf[outdf['condition'] == 'CantonS']
TE_df = TE_df[TE_df['int_2'] > dimm_thres]
TE_df = TE_df[TE_df['in_plane'] == True]

coords = np.vstack([TE_df['pca_x'], TE_df['pca_y']])
kde = gaussian_kde(coords, bw_method='silverman')

# Evaluate KDE on a grid to decide on a threshold for TE classification
X, Y = np.mgrid[-1:1:1000j, -1:1:1000j]
grid_coords = np.vstack([X.ravel(), Y.ravel()])
density = kde(grid_coords)
Z = np.reshape(density, X.shape)
threshold = Z.max() * 0.50

# Classify all cells as inROI or not based on KDE density
outdf['inROI'] = kde(outdf[['pca_x', 'pca_y']].to_numpy().T) > threshold

roidf = outdf[outdf['inROI'] == True][outdf['in_plane'] == True].copy()
roidf['TE'] = roidf['int_2'] > dimm_thres
roidf['RFP'] = roidf['int_0'] > rfp_thres
  
if not os.path.isdir('result'):
    os.mkdir('result')
      
# Compute TE ratio among RFP+ cells in the ROI for each sample
rtab = pd.crosstab(
  roidf[roidf['RFP'] == True]['TE'], roidf[roidf['RFP'] == True]['sample']
).T
rtab['ratio'] = rtab[True] / (rtab[True] + rtab[False])
rtab['condition'] = [i.split('_')[0] for i in rtab.index.tolist()]

outdf.to_csv(f'result/all_bsh_cells.csv', index=False)
roidf.to_csv(f'result/in_roi_bsh_cells.csv', index=False)
rtab.to_csv(f'result/te_ratios.csv', index=True)
