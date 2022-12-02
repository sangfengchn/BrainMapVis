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

%% global descriptive stats - 3T
metrics_dir = [proj_dir '/metrics'];
fmap = [metrics_dir '/dualreg_cc_slow4_3t.mat'];
tmpdr = load(fmap);
%space
dr_cc = [tmpdr.DR_cort_lh; tmpdr.DR_cort_rh];
dr_cc(dr_cc==0) = nan;
gDR_3t_cort_mean = squeeze(mean(dr_cc,1,'omitnan'));
gDR_3t_cort_std = squeeze(std(dr_cc,0,1,'omitnan'));            
%time
numRSNs = 7; numSess = 4;
DR_timeseries = tmpdr.DR_timeseries;
gDR_3t_ts_mean = zeros(numRSNs,nsubs,numSess);
gDR_3t_ts_std = zeros(numRSNs,nsubs,numSess);
DR_3t_tsCoRR = zeros(numRSNs,numRSNs,nsubs,numSess);
gDR_3t_tscorr_mean = zeros(nsubs,numSess,2);%pos vs neg
gDR_3t_tscorr_std = zeros(nsubs,numSess,2);%pos vs neg
for idxSub=1:nsubs
    for idxSess=1:numSess
        tmpts = DR_timeseries{idxSub,idxSess};
        tmpR = tril(ccs_core_fastcorr(tmpts,tmpts),-1);
        DR_3t_tsCoRR(:,:,idxSub,idxSess) = tmpR;
        gDR_3t_tscorr_mean(idxSub,idxSess,1) = mean(tmpR(tmpR>0));
        gDR_3t_tscorr_mean(idxSub,idxSess,2) = mean(tmpR(tmpR<0));
        gDR_3t_tscorr_std(idxSub,idxSess,1) = std(tmpR(tmpR>0));
        gDR_3t_tscorr_std(idxSub,idxSess,2) = std(tmpR(tmpR<0));
        gDR_3t_ts_mean(:,idxSub,idxSess) = mean(tmpts);
        gDR_3t_ts_std(:,idxSub,idxSess) = std(tmpts);
    end
end
fmap = [metrics_dir '/dualreg_tscorr_slow4_3t.mat'];
save(fmap,'DR_3t_tsCoRR')

%% global descriptive stats - 7T
metrics_dir = [proj_dir '/metrics'];
fmap = [metrics_dir '/dualreg_cc_slow4_7t.mat'];
tmpdr = load(fmap);
%space
dr_cc = [tmpdr.DR_cort_lh; tmpdr.DR_cort_rh];
dr_cc(dr_cc==0) = nan;
gDR_7t_cort_mean = squeeze(mean(dr_cc,1,'omitnan'));
gDR_7t_cort_std = squeeze(std(dr_cc,0,1,'omitnan'));            
%time
numRSNs = 7; numSess = 4;
DR_timeseries = tmpdr.DR_timeseries;
gDR_7t_ts_mean = zeros(numRSNs,nsubs,numSess);
gDR_7t_ts_std = zeros(numRSNs,nsubs,numSess);
DR_7t_tsCoRR = zeros(numRSNs,numRSNs,nsubs,numSess);
gDR_7t_tscorr_mean = zeros(nsubs,numSess,2);%pos vs neg
gDR_7t_tscorr_std = zeros(nsubs,numSess,2);%pos vs neg
for idxSub=1:nsubs
    for idxSess=1:numSess
        tmpts = DR_timeseries{idxSub,idxSess};
        tmpR = tril(ccs_core_fastcorr(tmpts,tmpts),-1);
        DR_7t_tsCoRR(:,:,idxSub,idxSess) = tmpR;
        gDR_7t_tscorr_mean(idxSub,idxSess,1) = mean(tmpR(tmpR>0));
        gDR_7t_tscorr_mean(idxSub,idxSess,2) = mean(tmpR(tmpR<0));
        gDR_7t_tscorr_std(idxSub,idxSess,1) = std(tmpR(tmpR>0));
        gDR_7t_tscorr_std(idxSub,idxSess,2) = std(tmpR(tmpR<0));
        gDR_7t_ts_mean(:,idxSub,idxSess) = mean(tmpts);
        gDR_7t_ts_std(:,idxSub,idxSess) = std(tmpts);
    end
end
fmap = [metrics_dir '/dualreg_tscorr_slow4_7t.mat'];
save(fmap,'DR_7t_tsCoRR')

%% save data for prism graph: gDR the spatial maps
%REST1 - day1
gmDR_REST1_3t = squeeze(mean(gDR_3t_cort_mean(:,:,1:2),3));
gmDR_REST1_7t = squeeze(mean(gDR_7t_cort_mean(:,:,1:2),3));
gsDR_REST1_3t = squeeze(mean(gDR_3t_cort_std(:,:,1:2),3));
gsDR_REST1_7t = squeeze(mean(gDR_7t_cort_std(:,:,1:2),3));
for idxRSN=1:numRSNs
    tmpmat = [gmDR_REST1_3t(idxRSN,:)' gmDR_REST1_7t(idxRSN,:)'];
    fprism = [fig_dir '/gmDR_REST1_RSN' num2str(idxRSN) '.csv'];
    writematrix(tmpmat,fprism)
    tmpmat = [gsDR_REST1_3t(idxRSN,:)' gsDR_REST1_7t(idxRSN,:)'];
    fprism = [fig_dir '/gsDR_REST1_RSN' num2str(idxRSN) '.csv'];
    writematrix(tmpmat,fprism)
end

%REST2 - day2
gmDR_REST2_3t = squeeze(mean(gDR_3t_cort_mean(:,:,3:4),3));
gmDR_REST2_7t = squeeze(mean(gDR_7t_cort_mean(:,:,3:4),3));
gsDR_REST2_3t = squeeze(mean(gDR_3t_cort_std(:,:,3:4),3));
gsDR_REST2_7t = squeeze(mean(gDR_7t_cort_std(:,:,3:4),3));
for idxRSN=1:numRSNs
    tmpmat = [gmDR_REST2_3t(idxRSN,:)' gmDR_REST2_7t(idxRSN,:)'];
    fprism = [fig_dir '/gmDR_REST2_RSN' num2str(idxRSN) '.csv'];
    writematrix(tmpmat,fprism)
    tmpmat = [gsDR_REST2_3t(idxRSN,:)' gsDR_REST2_7t(idxRSN,:)'];
    fprism = [fig_dir '/gsDR_REST2_RSN' num2str(idxRSN) '.csv'];
    writematrix(tmpmat,fprism)
end

%% save data for prism graph: gDR the time series
%REST1 - day1
gmDR_REST1_3t = squeeze(mean(gDR_3t_ts_mean(:,:,1:2),3));
gmDR_REST1_7t = squeeze(mean(gDR_7t_ts_mean(:,:,1:2),3));
gsDR_REST1_3t = squeeze(mean(gDR_3t_ts_std(:,:,1:2),3));
gsDR_REST1_7t = squeeze(mean(gDR_7t_ts_std(:,:,1:2),3));
for idxRSN=1:numRSNs
    tmpmat = [gmDR_REST1_3t(idxRSN,:)' gmDR_REST1_7t(idxRSN,:)'];
    fprism = [fig_dir '/gmDR_REST1_ts_RSN' num2str(idxRSN) '.csv'];
    writematrix(tmpmat,fprism)
    tmpmat = [gsDR_REST1_3t(idxRSN,:)' gsDR_REST1_7t(idxRSN,:)'];
    fprism = [fig_dir '/gsDR_REST1_ts_RSN' num2str(idxRSN) '.csv'];
    writematrix(tmpmat,fprism)
end

%REST2 - day2
gmDR_REST2_3t = squeeze(mean(gDR_3t_ts_mean(:,:,3:4),3));
gmDR_REST2_7t = squeeze(mean(gDR_7t_ts_mean(:,:,3:4),3));
gsDR_REST2_3t = squeeze(mean(gDR_3t_ts_std(:,:,3:4),3));
gsDR_REST2_7t = squeeze(mean(gDR_7t_ts_std(:,:,3:4),3));
for idxRSN=1:numRSNs
    tmpmat = [gmDR_REST2_3t(idxRSN,:)' gmDR_REST2_7t(idxRSN,:)'];
    fprism = [fig_dir '/gmDR_REST2_ts_RSN' num2str(idxRSN) '.csv'];
    writematrix(tmpmat,fprism)
    tmpmat = [gsDR_REST2_3t(idxRSN,:)' gsDR_REST2_7t(idxRSN,:)'];
    fprism = [fig_dir '/gsDR_REST2_ts_RSN' num2str(idxRSN) '.csv'];
    writematrix(tmpmat,fprism)
end

%% save data for prism graph: gDR the temporal correlation
%REST1 - day1 - pos corr
gmprDR_REST1_3t = squeeze(mean(gDR_3t_tscorr_mean(:,1:2,1),2));
gmprDR_REST1_7t = squeeze(mean(gDR_7t_tscorr_mean(:,1:2,1),2));
gsprDR_REST1_3t = squeeze(mean(gDR_3t_tscorr_std(:,1:2,1),2));
gsprDR_REST1_7t = squeeze(mean(gDR_7t_tscorr_std(:,1:2,1),2));
tmpmat = [gmprDR_REST1_3t gmprDR_REST1_7t];
fprism = [fig_dir '/gmpDR_REST1_tscorr.csv'];
writematrix(tmpmat,fprism)
tmpmat = [gsprDR_REST1_3t gsprDR_REST1_7t];
fprism = [fig_dir '/gspDR_REST1_tscorr.csv'];
writematrix(tmpmat,fprism)

%REST1 - day1 - neg corr
gmnrDR_REST1_3t = squeeze(mean(gDR_3t_tscorr_mean(:,1:2,2),2));
gmnrDR_REST1_7t = squeeze(mean(gDR_7t_tscorr_mean(:,1:2,2),2));
gsnrDR_REST1_3t = squeeze(mean(gDR_3t_tscorr_std(:,1:2,2),2));
gsnrDR_REST1_7t = squeeze(mean(gDR_7t_tscorr_std(:,1:2,2),2));
tmpmat = [gmnrDR_REST1_3t gmnrDR_REST1_7t];
fprism = [fig_dir '/gmnDR_REST1_tscorr.csv'];
writematrix(tmpmat,fprism)
tmpmat = [gsnrDR_REST1_3t gsnrDR_REST1_7t];
fprism = [fig_dir '/gsnDR_REST1_tscorr.csv'];
writematrix(tmpmat,fprism)

%REST1 - day2 - pos corr
gmprDR_REST2_3t = squeeze(mean(gDR_3t_tscorr_mean(:,3:4,1),2));
gmprDR_REST2_7t = squeeze(mean(gDR_7t_tscorr_mean(:,3:4,1),2));
gsprDR_REST2_3t = squeeze(mean(gDR_3t_tscorr_std(:,3:4,1),2));
gsprDR_REST2_7t = squeeze(mean(gDR_7t_tscorr_std(:,3:4,1),2));
tmpmat = [gmprDR_REST2_3t gmprDR_REST2_7t];
fprism = [fig_dir '/gmpDR_REST2_tscorr.csv'];
writematrix(tmpmat,fprism)
tmpmat = [gsprDR_REST2_3t gsprDR_REST2_7t];
fprism = [fig_dir '/gspDR_REST2_tscorr.csv'];
writematrix(tmpmat,fprism)

%REST1 - day2 - neg corr
gmnrDR_REST2_3t = squeeze(mean(gDR_3t_tscorr_mean(:,3:4,2),2));
gmnrDR_REST2_7t = squeeze(mean(gDR_7t_tscorr_mean(:,3:4,2),2));
gsnrDR_REST2_3t = squeeze(mean(gDR_3t_tscorr_std(:,3:4,2),2));
gsnrDR_REST2_7t = squeeze(mean(gDR_7t_tscorr_std(:,3:4,2),2));
tmpmat = [gmnrDR_REST2_3t gmnrDR_REST2_7t];
fprism = [fig_dir '/gmnDR_REST2_tscorr.csv'];
writematrix(tmpmat,fprism)
tmpmat = [gsnrDR_REST2_3t gsnrDR_REST2_7t];
fprism = [fig_dir '/gsnDR_REST2_tscorr.csv'];
writematrix(tmpmat,fprism)
