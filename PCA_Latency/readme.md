# fMRI Latency PCA Analysis

A toolkit for fMRI latency analysis based on Complex PCA. Extracts temporal delay information between brain regions using Hilbert transform and PCA.

## Files

| File | Language | Description |
|------|----------|-------------|
| `cpca_latency_pipeline.py` | Python | Concatenates multi-subject time series, computes unified latency map |
| `run_pca_fmri.m` | MATLAB | Per-subject GM/WM latency analysis |
| `make_pca_nifti.py` | Python | Converts PCA results to NIfTI format |

## Method

1. Z-score normalization
2. Hilbert transform (analytic signal)
3. Complex PCA to extract first principal component
4. Compute phase angle of loadings as latency values

## Dependencies

### Python
```
numpy
scipy
h5py
nibabel
fbpca
```

### MATLAB
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

## Usage

### cpca_latency_pipeline.py

Concatenates time series from all subjects and performs a single PCA analysis.

```bash
python cpca_latency_pipeline.py
```

**Input:** `E:\coding\DATA\*.mat` (4D fMRI data, dimension `[X, Y, Z, T]`)

**Output:** `latency_concatenated.npy` (whole-brain latency map)

### run_pca_fmri.m

Performs GM/WM/whole-brain PCA latency analysis for each subject separately.

```matlab
run_pca_fmri
```

**Input:**
- `E:\coding\DATA\*.mat` (4D fMRI data)
- `E:\coding\DATA\gm_mask_mean.mat` (gray matter probability map)
- `E:\coding\DATA\wm_mask_mean.mat` (white matter probability map)

**Output:** `E:\coding\DATA\pca_output\`
- `{subj}_latency.mat` - latency data
- `{subj}_GM_latency.png` - gray matter latency map
- `{subj}_WM_latency.png` - white matter latency map
- `{subj}_GM_WB_latency.png` - GM portion from whole-brain PCA
- `{subj}_WM_WB_latency.png` - WM portion from whole-brain PCA

### make_pca_nifti.py

Converts pickle format PCA results to NIfTI format.

```bash
python make_pca_nifti.py
```

**Input:** `E:\coding\DATA\*.pkl`

**Output:** `E:\coding\DATA\nifti_output\{subj}_component_{n}.nii.gz`

## Data Format

### Input Requirements

- 4D fMRI data: `.mat` format, dimension `[X, Y, Z, T]`

### Output Latency Map

- Value range: `[-π, π]` radians
- Positive values indicate relative delay, negative values indicate relative advance

## Two Analysis Strategies

| Strategy | Script | Use Case |
|----------|--------|----------|
| Concatenated | `cpca_latency_pipeline.py` | Group-level unified latency pattern |
| Per-subject | `run_pca_fmri.m` | Preserves individual differences for statistical analysis |

## Directory Structure

```
E:\coding\DATA\
├── subject1.mat
├── subject2.mat
├── ...
├── gm_mask_mean.mat
├── wm_mask_mean.mat
├── nifti_output\
│   └── {subj}_component_{n}.nii.gz
└── pca_output\
    ├── {subj}_latency.mat
    └── {subj}_*_latency.png
```

