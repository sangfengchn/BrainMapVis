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
data_dir = [hcp_dir '/1200'];
rest_labels = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', ...
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};
num_sess = numel(rest_labels);
idx_subcort = [3:9 12:21]; 
num_subcort = numel(idx_subcort);
bands_label = {'slow6','slow5','slow4','slow3','slow2','slow1'};
num_bands = numel(bands_label);
%ALFF
ALFF_cort_lh = zeros(nVertices_lh,nsubs,num_sess);
ALFF_cort_rh = zeros(nVertices_rh,nsubs,num_sess);
%ReHo
ReHo_cort_lh = zeros(nVertices_lh,nsubs,num_sess);
ReHo_cort_rh = zeros(nVertices_rh,nsubs,num_sess);
%VMHC
VMHC_cort_lh = zeros(nVertices_lh,nsubs,num_sess);
VMHC_cort_rh = zeros(nVertices_rh,nsubs,num_sess);
%loop subjects
for idx_sub=1:nsubs
    subject = num2str(subjects_id(idx_sub));
    sub_dir = [data_dir '/' subject '/MNINonLinear/Results'];
    for idx_rest=1:num_sess
        func_dir = [sub_dir '/' rest_labels{idx_rest}];
        for idx_bd=3 %slow-4
            disp(['subject ' subject ': ' rest_labels{idx_rest} ...
                ' maps from ' bands_label{idx_bd} ' ...'])
            %alff
            fmap = [func_dir '/boldalff_cort.' bands_label{idx_bd} '.mat'];
            tmpalff = load(fmap);
            ALFF_cort_lh(:,idx_sub,idx_rest) = tmpalff.alff_cc_lh(:,1);
            ALFF_cort_rh(:,idx_sub,idx_rest) = tmpalff.alff_cc_rh(:,1);
            %reho
            fmap = [func_dir '/boldreho_cort.' bands_label{idx_bd} '.mat'];
            tmpreho = load(fmap);
            ReHo_cort_lh(:,idx_sub,idx_rest) = tmpreho.reho_cc_lh(:,1);
            ReHo_cort_rh(:,idx_sub,idx_rest) = tmpreho.reho_cc_rh(:,1);
            %vmhc
            fmap = [func_dir '/boldvmhc_cort.' bands_label{idx_bd} '.mat'];
            tmpvmhc = load(fmap);
            VMHC_cort_lh(:,idx_sub,idx_rest) = tmpvmhc.vmhc_cc_lh(:,1);
            VMHC_cort_rh(:,idx_sub,idx_rest) = tmpvmhc.vmhc_cc_rh(:,1);
        end
    end
end

%% save data for subsequent stats
falff = [proj_dir '/metrics/alff_cc_slow4_3t.mat'];
save(falff,'ALFF_cort_lh','ALFF_cort_rh');
freho = [proj_dir '/metrics/reho_cc_slow4_3t.mat'];
save(freho,'ReHo_cort_lh','ReHo_cort_rh');
fvmhc = [proj_dir '/metrics/vmhc_cc_slow4_3t.mat'];
save(fvmhc,'VMHC_cort_lh','VMHC_cort_rh');
