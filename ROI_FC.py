

import os
import nibabel as nib
import numpy as np
import pandas as pd

def load_nifti_file(file_path):
    return nib.load(file_path).get_fdata()

def compute_functional_connectivity(time_series):
    corr_matrix = np.corrcoef(time_series)
    return corr_matrix

def fisher_r_to_z(corr_matrix):
    
    z_matrix = np.arctanh(corr_matrix)
    return z_matrix

def compute_overlap(roi, atlas_region):
    intersection = np.logical_and(roi, atlas_region)
    return np.sum(intersection) / np.sum(atlas_region)

def main(fmri_dir, roi_dir, atlas_path, additional_roi_path, output_dir):
    atlas = load_nifti_file(atlas_path)
    additional_roi = load_nifti_file(additional_roi_path)
    atlas_regions = [atlas == i for i in range(1, 1001)]
    atlas_regions.append(additional_roi)  

    
    fmri_files = sorted([f for f in os.listdir(fmri_dir) if f.endswith('.nii')])
    roi_files = sorted([f for f in os.listdir(roi_dir) if f.endswith('.nii')])

    for fmri_file, roi_file in zip(fmri_files, roi_files):
        fmri_path = os.path.join(fmri_dir, fmri_file)
        roi_path = os.path.join(roi_dir, roi_file)

        fmri_data = load_nifti_file(fmri_path)  
        roi_data = load_nifti_file(roi_path)

        time_series = []
        valid_regions = []

        for i, region in enumerate(atlas_regions):
            region = np.array(region, dtype=bool)
            region_voxels = fmri_data[region]  # region_voxels shape: (num_voxels, T)
            if region_voxels.shape[0] == 0:
                continue
            region_time_series = np.mean(region_voxels, axis=0)  # region_time_series shape: (T,)
            time_series.append(region_time_series)
            valid_regions.append(i)

        time_series = np.array(time_series)  # time_series shape: (num_regions, T)
        connectivity_matrix = compute_functional_connectivity(time_series)
        z_connectivity_matrix = fisher_r_to_z(connectivity_matrix)

        for idx in valid_regions:
            if idx == len(atlas_regions) - 1:
                # Skip the additional ROI
                continue

            overlap = compute_overlap(roi_data, atlas_regions[idx])
            if overlap > 0.95:
                z_connectivity_matrix[idx, :] = 0
                z_connectivity_matrix[:, idx] = 0

        np.fill_diagonal(z_connectivity_matrix, 0) 

        output_file = os.path.join(output_dir, f'{os.path.splitext(fmri_file)[0]}_connectivity.csv')
        df = pd.DataFrame(z_connectivity_matrix)
        df.to_csv(output_file, index=False, header=False)

if __name__ == "__main__":
    fmri_dir = r'G:\depression_repair\sanzu\zhu\rsdata\1'
    roi_dir = r'G:\depression_repair\sanzu\zhu\gs_Resliced'
    atlas_path = r'G:\depression_repair\Resliced_Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_1mm.nii'
    additional_roi_path = r'G:\depression_repair\Resliced_321236.nii' 
    output_dir = r'G:\depression_repair\sanzu\zhu\1.connqugsi'
    os.makedirs(output_dir, exist_ok=True)
    main(fmri_dir, roi_dir, atlas_path, additional_roi_path, output_dir)

