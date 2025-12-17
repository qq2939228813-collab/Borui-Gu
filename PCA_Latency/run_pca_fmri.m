%% Batch Latency PCA with GM / WM mask
clear; clc;

data_folder = 'E:\coding\DATA';
output_folder = fullfile(data_folder, 'pca_output');
mask_folder = data_folder;

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% Load masks
wm = load(fullfile(mask_folder, 'wm_mask_mean.mat'));
gm = load(fullfile(mask_folder, 'gm_mask_mean.mat'));

wm = squeeze(struct2array(wm));
gm = squeeze(struct2array(gm));

wm_mask = (wm > 0.8);
gm_mask = (gm > 0.5);

wm_mask = wm_mask(:);
gm_mask = gm_mask(:);

fprintf('WM mask vox = %d, GM mask vox = %d\n', sum(wm_mask), sum(gm_mask));

%% Get all mat files
file_list = dir(fullfile(data_folder, '*.mat'));
file_list = file_list(~contains({file_list.name}, 'mask'));

total_files = length(file_list);
fprintf('Found %d fMRI files\n', total_files);

%% Process each subject
for subj_idx = 1:total_files
    fname = file_list(subj_idx).name;
    fpath = fullfile(data_folder, fname);
    [~, subj_name, ~] = fileparts(fname);
    
    fprintf('\n[%d/%d] Processing: %s\n', subj_idx, total_files, fname);
    
    fmri_data = load(fpath);
    fn = fieldnames(fmri_data);
    fmri = fmri_data.(fn{1});
    
    [X, Y, Z, T] = size(fmri);
    fprintf('  Size = [%d %d %d %d]\n', X, Y, Z, T);
    
    fmri2 = reshape(fmri, [], T)';
    
    % GM-only PCA
    GM_vol = run_mask_pca(fmri2, gm_mask, X, Y, Z);
    gm_outname = fullfile(output_folder, [subj_name '_GM_latency.png']);
    plot_latency(GM_vol, ['GM Latency - ' subj_name], gm_outname);
    
    % WM-only PCA
    WM_vol = run_mask_pca(fmri2, wm_mask, X, Y, Z);
    wm_outname = fullfile(output_folder, [subj_name '_WM_latency.png']);
    plot_latency(WM_vol, ['WM Latency - ' subj_name], wm_outname);
    
    % Whole brain PCA
    wb_mask = gm_mask | wm_mask;
    WB_vol = run_mask_pca(fmri2, wb_mask, X, Y, Z);
    
    GM_WB = WB_vol .* reshape(gm_mask, [X,Y,Z]);
    gm_wb_outname = fullfile(output_folder, [subj_name '_GM_WB_latency.png']);
    plot_latency(GM_WB, ['GM from WB - ' subj_name], gm_wb_outname);
    
    WM_WB = WB_vol .* reshape(wm_mask, [X,Y,Z]);
    wm_wb_outname = fullfile(output_folder, [subj_name '_WM_WB_latency.png']);
    plot_latency(WM_WB, ['WM from WB - ' subj_name], wm_wb_outname);
    
    % Save latency data as .mat
    latency_data = struct();
    latency_data.GM_latency = GM_vol;
    latency_data.WM_latency = WM_vol;
    latency_data.WB_latency = WB_vol;
    latency_data.GM_from_WB = GM_WB;
    latency_data.WM_from_WB = WM_WB;
    
    mat_outname = fullfile(output_folder, [subj_name '_latency.mat']);
    save(mat_outname, '-struct', 'latency_data');
    fprintf('  Saved: %s\n', mat_outname);
    
    fprintf('  Done: %s\n', subj_name);
end

fprintf('\nAll %d subjects processed.\n', total_files);

%% PCA function
function latency3D = run_mask_pca(fmri2, mask, X, Y, Z)
    idx = find(mask);
    data_mask = fmri2(:, idx);
    
    good = std(data_mask) > 0 & all(isfinite(data_mask));
    data_mask = data_mask(:, good);
    idx = idx(good);
    
    data_mask = zscore(data_mask);
    data_hil = hilbert(data_mask);
    
    [coeff, ~, ~] = pca(data_hil, 'NumComponents', 1);
    loading = coeff(:,1);
    latency = angle(loading);
    
    vol = zeros(X*Y*Z, 1);
    vol(idx) = latency;
    latency3D = reshape(vol, [X, Y, Z]);
end

%% Plot function
function plot_latency(vol, title_str, pngname)
    slice_ids = 27 : 2 : (27 + 2*23);
    
    fig = figure('Position', [100 100 1200 800], 'Visible', 'off');
    colormap(jet);
    
    for i = 1:length(slice_ids)
        subplot(4, 6, i);
        img = rot90(vol(:,:,slice_ids(i)), 1);
        imagesc(img);
        axis image off;
        title(sprintf('Z=%d', slice_ids(i)));
    end
    
    h = colorbar('Position',[0.92 0.1 0.02 0.8]);
    ylabel(h, 'Latency (rad)');
    sgtitle(title_str, 'FontSize', 14);
    
    saveas(fig, pngname);
    close(fig);
end
