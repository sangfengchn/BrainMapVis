%% setup major directories
clear; clc
disk_dir = '/Volumes/LaCie_1';
hcp_dir = [disk_dir '/HCP'];
parc_dir = [disk_dir '/Schaefer2018_LocalGlobal'];
proj_dir = [hcp_dir '/zprojects/reliability'];
surf_dir = [hcp_dir '/zprojects/surfmodels'];
fig_dir = [proj_dir '/figures'];

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
numRSNs = 7; numSess = 4;

%% global descriptive stats - 3T
metrics_dir = [proj_dir '/metrics'];
fmap = [metrics_dir '/dualreg_cc_slow4_3t.mat'];
dualreg_cc = load(fmap);
gDR_3t_cort_mean = zeros(numRSNs,nsubs,numSess);
gDR_3t_cort_std = zeros(numRSNs,nsubs,numSess);
%building mask
dualreg_mask_lh = ones(nVertices_lh,1);
dualreg_mask_rh = ones(nVertices_lh,1);
for idxRSN=1:numRSNs
    for subid=1:nsubs
        for sessid=1:4
            tmpdr_lh = squeeze(dualreg_cc.DR_cort_lh(:,idxRSN,subid,sessid));
            tmpdr_rh = squeeze(dualreg_cc.DR_cort_rh(:,idxRSN,subid,sessid));
            dualreg_mask_lh(tmpdr_lh==0) = 0;
            dualreg_mask_rh(tmpdr_rh==0) = 0;
        end
    end
end
dualreg_mask = [dualreg_mask_lh; dualreg_mask_rh];
fmask = [metrics_dir '/dualreg_cc_slow4_mask_3t.mat'];
save(fmask,'dualreg_mask')
%building zmaps
nVertices = nVertices_lh + nVertices_rh;
dualreg_cc_zmap = zeros(nVertices,numRSNs,nsubs,numSess);    
for idxRSN=1:numRSNs
    tmpDR = [squeeze(dualreg_cc.DR_cort_lh(:,idxRSN,:,:)); ...
        squeeze(dualreg_cc.DR_cort_rh(:,idxRSN,:,:))];
    tmpDR_masked = tmpDR(dualreg_mask==1,:,:);
    [tmpZ,~,~] = zscore(tmpDR_masked);
    dualreg_cc_zmap(dualreg_mask==1,idxRSN,:,:) = tmpZ;
end
fzmap = [metrics_dir '/dualreg_cc_slow4_zmap_3t.mat'];
save(fzmap,'dualreg_cc_zmap')

%% global descriptive stats - 7T
metrics_dir = [proj_dir '/metrics'];
fmap = [metrics_dir '/dualreg_cc_slow4_7t.mat'];
dualreg_cc = load(fmap);
gDR_7t_cort_mean = zeros(numRSNs,nsubs,numSess);
gDR_7t_cort_std = zeros(numRSNs,nsubs,numSess);
%building mask
dualreg_mask_lh = ones(nVertices_lh,1);
dualreg_mask_rh = ones(nVertices_lh,1);
for idxRSN=1:numRSNs
    for subid=1:nsubs
        for sessid=1:4
            tmpdr_lh = squeeze(dualreg_cc.DR_cort_lh(:,idxRSN,subid,sessid));
            tmpdr_rh = squeeze(dualreg_cc.DR_cort_rh(:,idxRSN,subid,sessid));
            dualreg_mask_lh(tmpdr_lh==0) = 0;
            dualreg_mask_rh(tmpdr_rh==0) = 0;
        end
    end
end
dualreg_mask = [dualreg_mask_lh; dualreg_mask_rh];
fmask = [metrics_dir '/dualreg_cc_slow4_mask_7t.mat'];
save(fmask,'dualreg_mask')
%building zmaps
nVertices = nVertices_lh + nVertices_rh;
dualreg_cc_zmap = zeros(nVertices,numRSNs,nsubs,numSess);
for idxRSN=1:numRSNs
    tmpDR = [squeeze(dualreg_cc.DR_cort_lh(:,idxRSN,:,:)); ...
        squeeze(dualreg_cc.DR_cort_rh(:,idxRSN,:,:))];
    tmpDR_masked = tmpDR(dualreg_mask==1,:,:);
    [tmpZ,~,~] = zscore(tmpDR_masked);
    dualreg_cc_zmap(dualreg_mask==1,idxRSN,:,:) = tmpZ;
end
fzmap = [metrics_dir '/dualreg_cc_slow4_zmap_7t.mat'];
save(fzmap,'dualreg_cc_zmap')

%% save data for prism graph: gDR the spatial maps
% %REST1 - day1
% gmDR_REST1_3t = squeeze(mean(gDR_3t_cort_mean(:,:,1:2),3));
% gmDR_REST1_7t = squeeze(mean(gDR_7t_cort_mean(:,:,1:2),3));
% gsDR_REST1_3t = squeeze(mean(gDR_3t_cort_std(:,:,1:2),3));
% gsDR_REST1_7t = squeeze(mean(gDR_7t_cort_std(:,:,1:2),3));
% for idxRSN=1:numRSNs
%     tmpmat = [gmDR_REST1_3t(idxRSN,:)' gmDR_REST1_7t(idxRSN,:)'];
%     fprism = [fig_dir '/gmDR_REST1_RSN' num2str(idxRSN) '_check.csv'];
%     writematrix(tmpmat,fprism)
%     tmpmat = [gsDR_REST1_3t(idxRSN,:)' gsDR_REST1_7t(idxRSN,:)'];
%     fprism = [fig_dir '/gsDR_REST1_RSN' num2str(idxRSN) '_check.csv'];
%     writematrix(tmpmat,fprism)
% end
% 
% %REST2 - day2
% gmDR_REST2_3t = squeeze(mean(gDR_3t_cort_mean(:,:,3:4),3));
% gmDR_REST2_7t = squeeze(mean(gDR_7t_cort_mean(:,:,3:4),3));
% gsDR_REST2_3t = squeeze(mean(gDR_3t_cort_std(:,:,3:4),3));
% gsDR_REST2_7t = squeeze(mean(gDR_7t_cort_std(:,:,3:4),3));
% for idxRSN=1:numRSNs
%     tmpmat = [gmDR_REST2_3t(idxRSN,:)' gmDR_REST2_7t(idxRSN,:)'];
%     fprism = [fig_dir '/gmDR_REST2_RSN' num2str(idxRSN) '_check.csv'];
%     writematrix(tmpmat,fprism)
%     tmpmat = [gsDR_REST2_3t(idxRSN,:)' gsDR_REST2_7t(idxRSN,:)'];
%     fprism = [fig_dir '/gsDR_REST2_RSN' num2str(idxRSN) '_check.csv'];
%     writematrix(tmpmat,fprism)
% end
