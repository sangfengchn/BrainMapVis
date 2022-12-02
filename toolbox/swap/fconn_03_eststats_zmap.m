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
age = num(idx_subs_rest_3t7t,2);
sex = num(idx_subs_rest_3t7t,3);

%% global descriptive stats - 3T
metrics_dir = [proj_dir '/metrics'];
%ALFF
fmap = [metrics_dir '/alff_cc_slow4_3t.mat'];
alff_cc = load(fmap);
alff_mask_lh = ones(nVertices_lh,1);
alff_mask_rh = ones(nVertices_lh,1);
for subid=1:nsubs
    for sessid=1:4
        tmpalff_lh = alff_cc.ALFF_cort_lh(:,subid,sessid);
        tmpalff_rh = alff_cc.ALFF_cort_rh(:,subid,sessid);
        alff_mask_lh(tmpalff_lh==0) = 0;
        alff_mask_rh(tmpalff_rh==0) = 0;
    end
end
alff_mask = [alff_mask_lh; alff_mask_rh];
fmask = [metrics_dir '/alff_cc_slow4_mask_3t.mat'];
save(fmask,'alff_mask')
tmpALFF = [alff_cc.ALFF_cort_lh; alff_cc.ALFF_cort_rh];
tmpALFF_masked = tmpALFF(alff_mask==1,:,:);
[tmpZ,tmpmu,tmpsigma] = zscore(tmpALFF_masked);
alff_cc_zmap = zeros(size(tmpALFF));
alff_cc_zmap(alff_mask==1,:,:) = tmpZ;
fzmap = [metrics_dir '/alff_cc_slow4_zmap_3t.mat'];
save(fzmap,'alff_cc_zmap')
gALFF_3t_cort_mean = squeeze(tmpmu);
gALFF_3t_cort_std = squeeze(tmpsigma);
fgmap = [metrics_dir '/alff_cc_slow4_global_3t.mat'];
save(fgmap,'gALFF_3t_cort_mean','gALFF_3t_cort_std')
%ReHo
fmap = [metrics_dir '/reho_cc_slow4_3t.mat'];
reho_cc = load(fmap);
reho_mask_lh = ones(nVertices_lh,1);
reho_mask_rh = ones(nVertices_lh,1);
for subid=1:nsubs
    for sessid=1:4
        tmpreho_lh = reho_cc.ReHo_cort_lh(:,subid,sessid);
        tmpreho_rh = reho_cc.ReHo_cort_rh(:,subid,sessid);
        reho_mask_lh(tmpreho_lh==0) = 0;
        reho_mask_rh(tmpreho_rh==0) = 0;
    end
end
reho_mask = [reho_mask_lh; reho_mask_rh];
fmask = [metrics_dir '/reho_cc_slow4_mask_3t.mat'];
save(fmask,'reho_mask')
tmpReHo = [reho_cc.ReHo_cort_lh; reho_cc.ReHo_cort_rh];
tmpReHo_masked = tmpReHo(reho_mask==1,:,:);
[tmpZ,tmpmu,tmpsigma] = zscore(tmpReHo_masked);
reho_cc_zmap = zeros(size(tmpReHo));
reho_cc_zmap(reho_mask==1,:,:) = tmpZ;
fzmap = [metrics_dir '/reho_cc_slow4_zmap_3t.mat'];
save(fzmap,'reho_cc_zmap')
gReHo_3t_cort_mean = squeeze(tmpmu);
gReHo_3t_cort_std = squeeze(tmpsigma);
fgmap = [metrics_dir '/reho_cc_slow4_global_3t.mat'];
save(fgmap,'gReHo_3t_cort_mean','gReHo_3t_cort_std')
%VMHC
fmap = [metrics_dir '/vmhc_cc_slow4_3t.mat'];
vmhc_cc = load(fmap);
vmhc_mask_lh = ones(nVertices_lh,1);
vmhc_mask_rh = ones(nVertices_lh,1);
for subid=1:nsubs
    for sessid=1:4
        tmpvmhc_lh = vmhc_cc.VMHC_cort_lh(:,subid,sessid);
        tmpvmhc_rh = vmhc_cc.VMHC_cort_rh(:,subid,sessid);
        vmhc_mask_lh(tmpvmhc_lh==0) = 0;
        vmhc_mask_rh(tmpvmhc_rh==0) = 0;
    end
end
vmhc_mask = [vmhc_mask_lh; vmhc_mask_rh];
fmask = [metrics_dir '/vmhc_cc_slow4_mask_3t.mat'];
save(fmask,'vmhc_mask')
tmpVMHC = [vmhc_cc.VMHC_cort_lh; vmhc_cc.VMHC_cort_rh];
tmpVMHC_masked = tmpVMHC(vmhc_mask==1,:,:);
[tmpZ,tmpmu,tmpsigma] = zscore(tmpVMHC_masked);
vmhc_cc_zmap = zeros(size(tmpVMHC));
vmhc_cc_zmap(vmhc_mask==1,:,:) = tmpZ;
fzmap = [metrics_dir '/vmhc_cc_slow4_zmap_3t.mat'];
save(fzmap,'vmhc_cc_zmap')
gVMHC_3t_cort_mean = squeeze(tmpmu);
gVMHC_3t_cort_std = squeeze(tmpsigma);
fgmap = [metrics_dir '/vmhc_cc_slow4_global_3t.mat'];
save(fgmap,'gVMHC_3t_cort_mean','gVMHC_3t_cort_std')

%% global descriptive stats - 7T
metrics_dir = [proj_dir '/metrics'];
%ALFF
fmap = [metrics_dir '/alff_cc_slow4_7t.mat'];
alff_cc = load(fmap);
alff_mask_lh = ones(nVertices_lh,1);
alff_mask_rh = ones(nVertices_lh,1);
for subid=1:nsubs
    for sessid=1:4
        tmpalff_lh = alff_cc.ALFF_cort_lh(:,subid,sessid);
        tmpalff_rh = alff_cc.ALFF_cort_rh(:,subid,sessid);
        alff_mask_lh(tmpalff_lh==0) = 0;
        alff_mask_rh(tmpalff_rh==0) = 0;
    end
end
alff_mask = [alff_mask_lh; alff_mask_rh];
fmask = [metrics_dir '/alff_cc_slow4_mask_7t.mat'];
save(fmask,'alff_mask')
tmpALFF = [alff_cc.ALFF_cort_lh; alff_cc.ALFF_cort_rh];
tmpALFF_masked = tmpALFF(alff_mask==1,:,:);
[tmpZ,tmpmu,tmpsigma] = zscore(tmpALFF_masked);
alff_cc_zmap = zeros(size(tmpALFF));
alff_cc_zmap(alff_mask==1,:,:) = tmpZ;
fzmap = [metrics_dir '/alff_cc_slow4_zmap_7t.mat'];
save(fzmap,'alff_cc_zmap')
gALFF_7t_cort_mean = squeeze(tmpmu);
gALFF_7t_cort_std = squeeze(tmpsigma);
fgmap = [metrics_dir '/alff_cc_slow4_global_7t.mat'];
save(fgmap,'gALFF_7t_cort_mean','gALFF_7t_cort_std')
%ReHo
fmap = [metrics_dir '/reho_cc_slow4_7t.mat'];
reho_cc = load(fmap);
reho_mask_lh = ones(nVertices_lh,1);
reho_mask_rh = ones(nVertices_lh,1);
for subid=1:nsubs
    for sessid=1:4
        tmpreho_lh = reho_cc.ReHo_cort_lh(:,subid,sessid);
        tmpreho_rh = reho_cc.ReHo_cort_rh(:,subid,sessid);
        reho_mask_lh(tmpreho_lh==0) = 0;
        reho_mask_rh(tmpreho_rh==0) = 0;
    end
end
reho_mask = [reho_mask_lh; reho_mask_rh];
fmask = [metrics_dir '/reho_cc_slow4_mask_7t.mat'];
save(fmask,'reho_mask')
tmpReHo = [reho_cc.ReHo_cort_lh; reho_cc.ReHo_cort_rh];
tmpReHo_masked = tmpReHo(reho_mask==1,:,:);
[tmpZ,tmpmu,tmpsigma] = zscore(tmpReHo_masked);
reho_cc_zmap = zeros(size(tmpReHo));
reho_cc_zmap(reho_mask==1,:,:) = tmpZ;
fzmap = [metrics_dir '/reho_cc_slow4_zmap_7t.mat'];
save(fzmap,'reho_cc_zmap')
gReHo_7t_cort_mean = squeeze(tmpmu);
gReHo_7t_cort_std = squeeze(tmpsigma);
fgmap = [metrics_dir '/reho_cc_slow4_global_7t.mat'];
save(fgmap,'gReHo_7t_cort_mean','gReHo_7t_cort_std')
%VMHC
fmap = [metrics_dir '/vmhc_cc_slow4_7t.mat'];
vmhc_cc = load(fmap);
vmhc_mask_lh = ones(nVertices_lh,1);
vmhc_mask_rh = ones(nVertices_lh,1);
for subid=1:nsubs
    for sessid=1:4
        tmpvmhc_lh = vmhc_cc.VMHC_cort_lh(:,subid,sessid);
        tmpvmhc_rh = vmhc_cc.VMHC_cort_rh(:,subid,sessid);
        vmhc_mask_lh(tmpvmhc_lh==0) = 0;
        vmhc_mask_rh(tmpvmhc_rh==0) = 0;
    end
end
vmhc_mask = [vmhc_mask_lh; vmhc_mask_rh];
fmask = [metrics_dir '/vmhc_cc_slow4_mask_7t.mat'];
save(fmask,'vmhc_mask')
tmpVMHC = [vmhc_cc.VMHC_cort_lh; vmhc_cc.VMHC_cort_rh];
tmpVMHC_masked = tmpVMHC(vmhc_mask==1,:,:);
[tmpZ,tmpmu,tmpsigma] = zscore(tmpVMHC_masked);
vmhc_cc_zmap = zeros(size(tmpVMHC));
vmhc_cc_zmap(vmhc_mask==1,:,:) = tmpZ;
fzmap = [metrics_dir '/vmhc_cc_slow4_zmap_7t.mat'];
save(fzmap,'vmhc_cc_zmap')
gVMHC_7t_cort_mean = squeeze(tmpmu);
gVMHC_7t_cort_std = squeeze(tmpsigma);
fgmap = [metrics_dir '/vmhc_cc_slow4_global_7t.mat'];
save(fgmap,'gVMHC_7t_cort_mean','gVMHC_7t_cort_std')

%% save data for prism graph: gALFF
%REST1 - day1
gmALFF_REST1_3t = mean(gALFF_3t_cort_mean(:,1:2),2);
gmALFF_REST1_7t = mean(gALFF_7t_cort_mean(:,1:2),2);
gsALFF_REST1_3t = mean(gALFF_3t_cort_std(:,1:2),2);
gsALFF_REST1_7t = mean(gALFF_7t_cort_std(:,1:2),2);
%REST2 - day2
gmALFF_REST2_3t = mean(gALFF_3t_cort_mean(:,3:4),2);
gmALFF_REST2_7t = mean(gALFF_7t_cort_mean(:,3:4),2);
gsALFF_REST2_3t = mean(gALFF_3t_cort_std(:,3:4),2);
gsALFF_REST2_7t = mean(gALFF_7t_cort_std(:,3:4),2);

%% save data for prism graph: gReHo
%REST1 - day1
gmReHo_REST1_3t = mean(gReHo_3t_cort_mean(:,1:2),2);
gmReHo_REST1_7t = mean(gReHo_7t_cort_mean(:,1:2),2);
gsReHo_REST1_3t = mean(gReHo_3t_cort_std(:,1:2),2);
gsReHo_REST1_7t = mean(gReHo_7t_cort_std(:,1:2),2);
%REST2 - day2
gmReHo_REST2_3t = mean(gReHo_3t_cort_mean(:,3:4),2);
gmReHo_REST2_7t = mean(gReHo_7t_cort_mean(:,3:4),2);
gsReHo_REST2_3t = mean(gReHo_3t_cort_std(:,3:4),2);
gsReHo_REST2_7t = mean(gReHo_7t_cort_std(:,3:4),2);

%% save data for prism graph: gVMHC
%REST1 - day1
gmVMHC_REST1_3t = mean(gVMHC_3t_cort_mean(:,1:2),2);
gmVMHC_REST1_7t = mean(gVMHC_7t_cort_mean(:,1:2),2);
gsVMHC_REST1_3t = mean(gVMHC_3t_cort_std(:,1:2),2);
gsVMHC_REST1_7t = mean(gVMHC_7t_cort_std(:,1:2),2);
%REST2 - day2
gmVMHC_REST2_3t = mean(gVMHC_3t_cort_mean(:,3:4),2);
gmVMHC_REST2_7t = mean(gVMHC_7t_cort_mean(:,3:4),2);
gsVMHC_REST2_3t = mean(gVMHC_3t_cort_std(:,3:4),2);
gsVMHC_REST2_7t = mean(gVMHC_7t_cort_std(:,3:4),2);
