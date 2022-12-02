%% setup major directories
clear; clc
disk_dir = '/Volumes/LaCie_1';
hcp_dir = [disk_dir '/HCP'];
parc_dir = [disk_dir '/Schaefer2018_LocalGlobal'];
proj_dir = [hcp_dir '/zprojects/reliability'];
surf_dir = [hcp_dir '/zprojects/surfmodels'];
fig_dir = [proj_dir '/group/figures'];

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
    '/Conte69.L.midthickness.32k_fs_LR.surf.gii']);
nVertices_lh = size(Conte32k_lh.vertices,1);
Conte32k_rh = gifti([ConteAtlas_dir ...
    '/Conte69.R.midthickness.32k_fs_LR.surf.gii']);
nVertices_rh = size(Conte32k_rh.vertices,1);

%% read the subject list
subjects_list = [proj_dir '/info/subjects_3t7t.xlsx'];
[status,sheets,xlFormat] = xlsfinfo(subjects_list);
[num,txt,raw] = xlsread(subjects_list,sheets{1});
load([proj_dir '/info/subidx_3t7t_rest4sess.mat'])
subjects_id = num(idx_subs_rest_3t7t,1); 
nsubs = numel(subjects_id);

%% Subjects from HCP1200 release
data_dir = [hcp_dir '/7T'];
rest_labels = {'rfMRI_REST1_7T_PA', 'rfMRI_REST2_7T_AP', ...
    'rfMRI_REST3_7T_PA', 'rfMRI_REST4_7T_AP'};
num_sess = numel(rest_labels);
bands_label = {'slow6','slow5','slow4','slow3','slow2','slow1'};
num_bands = numel(bands_label);
DR_cort_lh = zeros(nVertices_lh,7,nsubs,num_sess);
DR_cort_rh = zeros(nVertices_rh,7,nsubs,num_sess);
DR_timeseries = cell(nsubs,num_sess);
%loop subjects
for idx_sub=1:nsubs
    subject = num2str(subjects_id(idx_sub));
    sub_dir = [data_dir '/' subject '/MNINonLinear/Results'];
    for idx_rest=1:num_sess
        func_dir = [sub_dir '/' rest_labels{idx_rest}];
        for idx_bd=3 %slow-4
            disp(['subject ' subject ': ' rest_labels{idx_rest} ...
                ' maps from ' bands_label{idx_bd} ' ...'])
            %dual regression
            fmap = [func_dir '/bolddualreg_cort.' bands_label{idx_bd} '.mat'];
            tmpdualreg = load(fmap);
            DR_cort_lh(:,:,idx_sub,idx_rest) = tmpdualreg.dr_cc_lh;
            DR_cort_rh(:,:,idx_sub,idx_rest) = tmpdualreg.dr_cc_rh;
            DR_timeseries{idx_sub,idx_rest} = tmpdualreg.dr_timeseries;
        end
    end
end

%% save data for subsequent stats
fdualreg = [proj_dir '/metrics/dualreg_cc_slow4_7t.mat'];
save(fdualreg,'DR_cort_lh','DR_cort_rh','DR_timeseries');
