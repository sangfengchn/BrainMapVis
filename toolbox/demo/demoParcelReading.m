%% major directory settings
clear all; clc
disk_dir = '/Volumes/LaCie';
hcp_dir = [disk_dir '/HCP'];
parc_dir = [disk_dir '/Schaefer2018_LocalGlobal'];
proj_dir = [hcp_dir '/zprojects/reliability'];
fig_dir = [proj_dir '/group/figures'];

%% major toolbox settings
fs_home = '/Applications/freesurfer/freesurfer53';
%connectome computation system
ccs_dir = '/Users/mac/Projects/CCS';
ccs_matlab = [ccs_dir '/matlab'];
ccs_vistool = [ccs_dir '/vistool'];
%hcp workbench
hcpwkbc_dir = [ccs_dir '/extool/hcpworkbench'];
wbcmd = [hcpwkbc_dir ...
    '/macosx64_apps/wb_command.app/Contents/MacOS/wb_command'];
atlas_dir = [hcpwkbc_dir '/resources/32k_ConteAtlas_v2'];
%cifti toolbox
cifti_matlab = [proj_dir '/matlab/cifti-matlab-master'];

%% add the paths to matlab
addpath(genpath(ccs_matlab)) %ccs matlab scripts
addpath(genpath(ccs_vistool)) %ccs matlab scripts
addpath(genpath(cifti_matlab)) %cifti paths
addpath(genpath([fs_home '/matlab'])) %freesurfer matlab scripts

%% read brain parcellation
num_parc = 400;
fparc = [parc_dir '/Parcellations/HCP/fslr32k/cifti' ...
    '/Schaefer2018_' num2str(num_parc) 'Parcels_7Networks_order.dlabel.nii'];
xlgp400 = ft_read_cifti(fparc,'mapname','array');
