import pickle
import numpy as np
import nibabel as nib
import os

DATA_FOLDER = r"E:\coding\DATA"
OUTPUT_FOLDER = r"E:\coding\DATA\nifti_output"

def process_one_subject(pkl_path, out_folder):
    subj_name = os.path.splitext(os.path.basename(pkl_path))[0]
    
    with open(pkl_path, "rb") as f:
        res = pickle.load(f)

    loadings = res["loadings"]
    mask = res["mask"]
    vol_shape = res["vol_shape"]

    n_components = loadings.shape[0]

    for i in range(n_components):
        comp = loadings[i]
        full_volume = np.zeros(mask.shape)
        full_volume[mask] = comp
        vol_3d = full_volume.reshape(vol_shape)
        img = nib.Nifti1Image(vol_3d, affine=np.eye(4))
        out_path = os.path.join(out_folder, f"{subj_name}_component_{i+1}.nii.gz")
        nib.save(img, out_path)

    return n_components

def main():
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    files = [f for f in os.listdir(DATA_FOLDER) if f.endswith(".pkl")]
    total = len(files)

    if total == 0:
        print("No .pkl files found in", DATA_FOLDER)
        return

    for idx, fname in enumerate(files):
        pkl_path = os.path.join(DATA_FOLDER, fname)
        n_comp = process_one_subject(pkl_path, OUTPUT_FOLDER)
        print(f"Processed {idx+1}/{total}: {fname} ({n_comp} components)")

    print("All subjects done.")

if __name__ == "__main__":
    main()
