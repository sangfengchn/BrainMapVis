%% setup major directories
clear; clc
disk_dir = '/Volumes/LaCie_1';
hcp_dir = [disk_dir '/HCP'];
parc_dir = [disk_dir '/Schaefer2018_LocalGlobal'];
proj_dir = [hcp_dir '/zprojects/reliability'];
surf_dir = [hcp_dir '/zprojects/surfmodels'];

%% setup major toolboxes
fs_home = '/Applications/freesurfer/7.2';
%connectome computation system
ccs_dir = '/Users/xinian.zuo/Projects/CCS';
ccs_matlab = [ccs_dir '/matlab'];
ccs_vistool = [ccs_dir '/vistool'];
%hcp workbench
hcpwkbc_dir = [hcp_dir '/workbench'];
wbcmd = [hcpwkbc_dir '/macosx64_apps/' ...
    'wb_command.app/Contents/MacOS/wb_command'];
ConteAtlas_dir = [surf_dir '/32k_ConteAtlas_v2'];
%cifti toolbox
cifti_matlab = [proj_dir '/matlab/cifti-matlab-master'];

%% add the paths to matlab
addpath(genpath(ccs_matlab)) %ccs matlab scripts
addpath(genpath(ccs_vistool)) %ccs matlab scripts
addpath(genpath(cifti_matlab)) %cifti paths
addpath(genpath([fs_home '/matlab'])) %freesurfer matlab scripts

%% load the geometry of the 32k_ConteAtlas
Conte32k_lh = gifti([ConteAtlas_dir ...
    '/Conte69.L.inflated.32k_fs_LR.surf.gii']);
nVertices_lh = size(Conte32k_lh.vertices,1);
surfConte69_lh.tri = Conte32k_lh.faces;
surfConte69_lh.coord = Conte32k_lh.vertices'; 

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

%% read the subject list
subjects_list = [proj_dir '/info/subjects_3t7t.xlsx'];
[status,sheets,xlFormat] = xlsfinfo(subjects_list);
[num,txt,raw] = xlsread(subjects_list,sheets{1});
subjects_id = num(:,1); nsubs = numel(subjects_id);

%% loop subjects from HCP1200 release
data_dir = [hcp_dir '/7T'];
rest_labels = {'rfMRI_REST1_7T_PA', 'rfMRI_REST2_7T_AP', ...
    'rfMRI_REST3_7T_PA', 'rfMRI_REST4_7T_AP'};
bands_label = {'slow6','slow5','slow4','slow3','slow2','slow1'};
num_bands = numel(bands_label);
metrics = {'alff','reho','vmhc'};
fig_dir = [proj_dir '/figures/100610'];
for idx_sub=1
    subject = num2str(subjects_id(idx_sub));
    sub_dir = [data_dir '/' subject '/MNINonLinear/Results'];
    for idx_rest=1:4
        disp(['subject ' subject ': ' rest_labels{idx_rest}])
        func_dir = [sub_dir '/' rest_labels{idx_rest}];
        if ~exist([fig_dir '/' rest_labels{idx_rest}],'dir')
            mkdir(fig_dir, rest_labels{idx_rest});
        end
        %loop bands
        for idx_bd=1:num_bands
            %alff
            fmetric = [func_dir '/boldalff_cort.' bands_label{idx_bd} '.mat'];
            if exist(fmetric,'file')
                load(fmetric); tmpstats = [alff_cc_lh; alff_cc_rh];
                cmax = max(tmpstats); cmin = min(tmpstats);
                if cmin<0 
                    cmap_range = [-cmax cmax]*0.8;
                    cmap_tmp = cmap_tdmean;
                else 
                    cmap_range = [0 cmax]*0.8;
                    cmap_tmp = cmap_sdmean;
                end
                figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
                SurfStatView(alff_cc_lh, surfConte69_lh, ' ', 'white', 'true'); 
                colormap(cmap_tmp); SurfStatColLim(cmap_range);
                set(gcf, 'PaperPositionMode', 'auto');            
                figout = [fig_dir '/' rest_labels{idx_rest} ...
                    '/alff_cc.' bands_label{idx_bd} '.lh.png'];
                print('-dpng', '-r300', figout); close
                figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
                SurfStatView(alff_cc_rh, surfConte69_rh, ' ', 'white', 'true'); 
                colormap(cmap_tmp); SurfStatColLim(cmap_range);
                set(gcf, 'PaperPositionMode', 'auto');
                figout = [fig_dir '/' rest_labels{idx_rest} ...
                    '/alff_cc.' bands_label{idx_bd} '.rh.png'];
                print('-dpng', '-r300', figout); close
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
                %vmhc
                fmetric = [func_dir '/boldvmhc_cort.' bands_label{idx_bd} '.mat'];
                load(fmetric); tmpstats = [vmhc_cc_lh; vmhc_cc_rh];
                cmax = max(tmpstats); cmin = min(tmpstats);
                if cmin<0 
                    cmap_range = [-cmax cmax]*0.8;
                    cmap_tmp = cmap_tdmean;
                else 
                    cmap_range = [0 cmax]*0.8;
                    cmap_tmp = cmap_sdmean;
                end
                figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
                SurfStatView(vmhc_cc_lh, surfConte69_lh, ' ', 'white', 'true'); 
                colormap(cmap_tmp); SurfStatColLim(cmap_range);
                set(gcf, 'PaperPositionMode', 'auto');            
                figout = [fig_dir '/' rest_labels{idx_rest} ...
                    '/vmhc_cc.' bands_label{idx_bd} '.lh.png'];
                print('-dpng', '-r300', figout); close
                figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
                SurfStatView(vmhc_cc_rh, surfConte69_rh, ' ', 'white', 'true'); 
                colormap(cmap_tmp); SurfStatColLim(cmap_range);
                set(gcf, 'PaperPositionMode', 'auto');
                figout = [fig_dir '/' rest_labels{idx_rest} ...
                    '/vmhc_cc.' bands_label{idx_bd} '.rh.png'];
                print('-dpng', '-r300', figout); close
            end
        end
    end
end
