import numpy as np
import os
import h5py
from scipy.stats import zscore
from scipy.signal import hilbert
import fbpca

def load_mat_4d(fp):
    with h5py.File(fp, "r") as f:
        key = list(f.keys())[0]
        arr = np.array(f[key])
    return arr

def prepare_data(img4d):
    shape = img4d.shape
    if len(shape) != 4:
        raise ValueError(f"Expected 4D fMRI, got shape {shape}")
    T = shape[0]
    spatial_shape = shape[1:]
    V = np.prod(spatial_shape)
    data = img4d.reshape(T, V)
    return data, spatial_shape

def get_valid_mask(data):
    valid_mask = np.isfinite(data).all(axis=0)
    std = data.std(axis=0)
    var_mask = std > 1e-8
    return valid_mask & var_mask

def do_complex_pca_first_component(data):
    data_z = zscore(data, axis=0)
    data_h = hilbert(data_z, axis=0).conj()
    U, s, Va = fbpca.pca(data_h, k=1)
    loadings = Va.T @ np.diag(s)
    loadings /= np.sqrt(data_h.shape[0] - 1)
    phase = np.angle(loadings.T)[0]
    return phase

def main():
    folder = r"E:\coding\DATA"
    files = [os.path.join(folder, x) for x in os.listdir(folder) if x.endswith(".mat")]
    total = len(files)
    
    if total == 0:
        print("No .mat files found in", folder)
        return

    print(f"Found {total} subjects, loading data...")

    all_data = []
    common_mask = None
    spatial_shape = None

    for idx, fp in enumerate(files):
        img4d = load_mat_4d(fp)
        data, sp_shape = prepare_data(img4d)
        
        if spatial_shape is None:
            spatial_shape = sp_shape
        
        mask = get_valid_mask(data)
        if common_mask is None:
            common_mask = mask
        else:
            common_mask = common_mask & mask
        
        all_data.append(data)
        print(f"Loaded {idx+1}/{total}: {os.path.basename(fp)}, T={data.shape[0]}")

    print(f"Common valid voxels: {common_mask.sum()}")

    print("Concatenating time series...")
    concat_data = []
    for data in all_data:
        concat_data.append(data[:, common_mask])
    concat_data = np.vstack(concat_data)
    print(f"Concatenated data shape: {concat_data.shape} (T_total x V)")

    print("Running complex PCA on concatenated data...")
    phase = do_complex_pca_first_component(concat_data)

    full_phase = np.zeros(common_mask.shape, dtype=np.float32)
    full_phase[common_mask] = phase
    vol_phase = full_phase.reshape(spatial_shape)

    np.save("latency_concatenated.npy", vol_phase)
    print("Done.")
    print("Latency map saved as latency_concatenated.npy")
    print("Shape:", vol_phase.shape)

if __name__ == "__main__":
    main()
