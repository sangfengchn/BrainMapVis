%% setup major directories
clear; clc
tmpinfo = what('./');
proj_dir = tmpinfo.path(1:43);
hcp_dir = [proj_dir '/data'];
surf_dir = [proj_dir '/matlab/surfmodels'];

%% setup major toolboxes
fs_home = [proj_dir '/matlab/freesurfer'];
%connectome computation system
ccs_matlab = [proj_dir '/matlab/core'];
ccs_vistool = [proj_dir '/matlab/vistlbx'];
%hcp workbench
hcpwkbc_dir = [proj_dir '/connectome_workbench/macos'];
wbcmd = [hcpwkbc_dir '/bin_macosx64/wb_command'];
ConteAtlas_dir = [surf_dir '/32k_ConteAtlas_v2'];
%cifti toolbox
cifti_matlab = [proj_dir '/matlab/cifti-matlab-master'];

%% add the paths to matlab
addpath(genpath(ccs_matlab)) %ccs matlab scripts
addpath(genpath(ccs_vistool)) %ccs matlab scripts
addpath(genpath(cifti_matlab)) %cifti paths
addpath(genpath(fs_home)) %freesurfer matlab scripts

%% load the geometry of the 32k_ConteAtlas
Conte32k_lh = gifti([ConteAtlas_dir ...
    '/Conte69.L.inflated.32k_fs_LR.surf.gii']);
nVertices_lh = size(Conte32k_lh.vertices,1);
surfConte69_lh.tri = Conte32k_lh.faces;
surfConte69_lh.coord = Conte32k_lh.vertices'; 
% right hemisphere
Conte32k_rh = gifti([ConteAtlas_dir ...
    '/Conte69.R.inflated.32k_fs_LR.surf.gii']);
nVertices_rh = size(Conte32k_rh.vertices,1);
surfConte69_rh.tri = Conte32k_rh.faces;
surfConte69_rh.coord = Conte32k_rh.vertices';

%% Load Colormaps
cmap_sdmean = jet(256); cmap_sdmean(1,:) = 0.5;
cmap_tdmean = jet(256); cmap_tdmean(127:129,:) = 0.5;
%matplot default
fcolor = [ccs_vistool '/colormaps_reference_00.hires.png'];
cmap_pus = ccs_extractcolormaps(fcolor);
cmap_var = cmap_pus{1}; cmap_var(1,:) = 0.5;
%caret single direction
fcolor = [ccs_vistool '/fimages/Purple-Red_Seg10_Caret.png'];
cmap_caret = ccs_mkcolormap(fcolor);
cmap_icc = cmap_caret; cmap_icc(1,:) = 0.5; 

%% loop subjects from HCP1200 release
data_dir = [hcp_dir '/1200'];
rest_labels = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', ...
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};
bands_label = {'slow6','slow5','slow4','slow3','slow2','slow1'};
num_bands = numel(bands_label);
metrics = {'alff','reho','vmhc'};
fig_dir = [proj_dir '/figures/100610'];
%cifti toolbox
surfstat_matlab = [proj_dir '/matlab/surfstat'];
addpath(genpath(surfstat_matlab))
subject = '100610';
sub_dir = [data_dir '/' subject '/MNINonLinear/Results'];
for idx_rest=1:4
    disp(['subject ' subject ': ' rest_labels{idx_rest}])
    func_dir = [sub_dir '/' rest_labels{idx_rest}];
    if ~exist([fig_dir '/' rest_labels{idx_rest}],'dir')
        mkdir(fig_dir, rest_labels{idx_rest});
    end
    %loop bands
    for idx_bd=3 %num_bands
            %reho
            fmetric = [func_dir '/boldreho_cort.' bands_label{idx_bd} '.mat'];
            load(fmetric); tmpstats = [reho_cc_lh; reho_cc_rh];
            cmax = max(tmpstats); cmin = min(tmpstats);
            if cmin<0 
                cmap_range = [-cmax cmax]*0.8;
                cmap_tmp = cmap_tdmean;
            else 
                cmap_range = [0 cmax]*0.8;
                cmap_tmp = cmap_sdmean;
            end
            figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
            SurfStatView(reho_cc_lh, surfConte69_lh, ' ', 'white', 'true'); 
            colormap(cmap_tmp); SurfStatColLim(cmap_range);
            set(gcf, 'PaperPositionMode', 'auto');            
            figout = [fig_dir '/' rest_labels{idx_rest} ...
                '/reho_cc.' bands_label{idx_bd} '.lh.png'];
            print('-dpng', '-r300', figout); close
            figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
            SurfStatView(reho_cc_rh, surfConte69_rh, ' ', 'white', 'true'); 
            colormap(cmap_tmp); SurfStatColLim(cmap_range);
            set(gcf, 'PaperPositionMode', 'auto');
            figout = [fig_dir '/' rest_labels{idx_rest} ...
                '/reho_cc.' bands_label{idx_bd} '.rh.png'];
            print('-dpng', '-r300', figout); close
    end
end
