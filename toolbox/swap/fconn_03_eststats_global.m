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

%% global descriptive stats - 3T
metrics_dir = [proj_dir '/metrics'];
%ALFF
fmap = [metrics_dir '/alff_cc_slow4_3t.mat'];
tmpalff = load(fmap);
alff_cc = [tmpalff.ALFF_cort_lh; tmpalff.ALFF_cort_rh];
alff_cc(alff_cc==0) = nan;
gALFF_3t_cort_mean = squeeze(mean(alff_cc,1,'omitnan'));
gALFF_3t_cort_std = squeeze(std(alff_cc,0,1,'omitnan'));            
%ReHo
fmap = [metrics_dir '/reho_cc_slow4_3t.mat'];
tmpreho = load(fmap);
reho_cc = [tmpreho.ReHo_cort_lh; tmpreho.ReHo_cort_rh];
reho_cc(reho_cc==0) = nan;
gReHo_3t_cort_mean = squeeze(mean(reho_cc,1,'omitnan'));
gReHo_3t_cort_std = squeeze(std(reho_cc,0,1,'omitnan'));
%VMHC
fmap = [metrics_dir '/vmhc_cc_slow4_3t.mat'];
tmpvmhc = load(fmap);
vmhc_cc = [tmpvmhc.VMHC_cort_lh; tmpvmhc.VMHC_cort_rh];
vmhc_cc(vmhc_cc==0) = nan;
gVMHC_3t_cort_mean = squeeze(mean(vmhc_cc,1,'omitnan'));
gVMHC_3t_cort_std = squeeze(std(vmhc_cc,0,1,'omitnan'));

%% global descriptive stats - 7T
metrics_dir = [proj_dir '/metrics'];
%ALFF
fmap = [metrics_dir '/alff_cc_slow4_7t.mat'];
tmpalff = load(fmap);
alff_cc = [tmpalff.ALFF_cort_lh; tmpalff.ALFF_cort_rh];
alff_cc(alff_cc==0) = nan;
gALFF_7t_cort_mean = squeeze(mean(alff_cc,1,'omitnan'));
gALFF_7t_cort_std = squeeze(std(alff_cc,0,1,'omitnan'));            
%ReHo
fmap = [metrics_dir '/reho_cc_slow4_7t.mat'];
tmpreho = load(fmap);
reho_cc = [tmpreho.ReHo_cort_lh; tmpreho.ReHo_cort_rh];
reho_cc(reho_cc==0) = nan;
gReHo_7t_cort_mean = squeeze(mean(reho_cc,1,'omitnan'));
gReHo_7t_cort_std = squeeze(std(reho_cc,0,1,'omitnan'));
%VMHC
fmap = [metrics_dir '/vmhc_cc_slow4_7t.mat'];
tmpvmhc = load(fmap);
vmhc_cc = [tmpvmhc.VMHC_cort_lh; tmpvmhc.VMHC_cort_rh];
vmhc_cc(vmhc_cc==0) = nan;
gVMHC_7t_cort_mean = squeeze(mean(vmhc_cc,1,'omitnan'));
gVMHC_7t_cort_std = squeeze(std(vmhc_cc,0,1,'omitnan'));

%% save data for prism graph: gALFF
%REST1 - day1
gmALFF_REST1_3t = squeeze(mean(gALFF_3t_cort_mean(:,1:2),2));
gmALFF_REST1_7t = squeeze(mean(gALFF_7t_cort_mean(:,1:2),2));
gsALFF_REST1_3t = squeeze(mean(gALFF_3t_cort_std(:,1:2),2));
gsALFF_REST1_7t = squeeze(mean(gALFF_7t_cort_std(:,1:2),2));
%REST2 - day2
gmALFF_REST2_3t = squeeze(mean(gALFF_3t_cort_mean(:,3:4),2));
gmALFF_REST2_7t = squeeze(mean(gALFF_7t_cort_mean(:,3:4),2));
gsALFF_REST2_3t = squeeze(mean(gALFF_3t_cort_std(:,3:4),2));
gsALFF_REST2_7t = squeeze(mean(gALFF_7t_cort_std(:,3:4),2));

%% save data for prism graph: gReHo
%REST1 - day1
gmReHo_REST1_3t = squeeze(mean(gReHo_3t_cort_mean(:,1:2),2));
gmReHo_REST1_7t = squeeze(mean(gReHo_7t_cort_mean(:,1:2),2));
gsReHo_REST1_3t = squeeze(mean(gReHo_3t_cort_std(:,1:2),2));
gsReHo_REST1_7t = squeeze(mean(gReHo_7t_cort_std(:,1:2),2));
%REST2 - day2
gmReHo_REST2_3t = squeeze(mean(gReHo_3t_cort_mean(:,3:4),2));
gmReHo_REST2_7t = squeeze(mean(gReHo_7t_cort_mean(:,3:4),2));
gsReHo_REST2_3t = squeeze(mean(gReHo_3t_cort_std(:,3:4),2));
gsReHo_REST2_7t = squeeze(mean(gReHo_7t_cort_std(:,3:4),2));

%% save data for prism graph: gVMHC
%REST1 - day1
gmVMHC_REST1_3t = squeeze(mean(gVMHC_3t_cort_mean(:,1:2),2));
gmVMHC_REST1_7t = squeeze(mean(gVMHC_7t_cort_mean(:,1:2),2));
gsVMHC_REST1_3t = squeeze(mean(gVMHC_3t_cort_std(:,1:2),2));
gsVMHC_REST1_7t = squeeze(mean(gVMHC_7t_cort_std(:,1:2),2));
%REST2 - day2
gmVMHC_REST2_3t = squeeze(mean(gVMHC_3t_cort_mean(:,3:4),2));
gmVMHC_REST2_7t = squeeze(mean(gVMHC_7t_cort_mean(:,3:4),2));
gsVMHC_REST2_3t = squeeze(mean(gVMHC_3t_cort_std(:,3:4),2));
gsVMHC_REST2_7t = squeeze(mean(gVMHC_7t_cort_std(:,3:4),2));
