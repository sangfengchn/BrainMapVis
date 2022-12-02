%% setup major directories
clear all; clc
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
    '/Conte69.L.inflated.32k_fs_LR.surf.gii']);
nVertices_lh = size(Conte32k_lh.vertices,1);
surfConte69_lh.tri = Conte32k_lh.faces;
surfConte69_lh.coord = Conte32k_lh.vertices'; 
%----------------------------------------------
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
%div
fcolor = [ccs_vistool '/colormaps_reference_03.hires.png'];
cmap_div = ccs_extractcolormaps(fcolor);
%blue-white-red
cmap_bwr = cmap_div{3}; %cmap_bwr(1,:) = 0;
%cool-warm
cmap_cw = cmap_div{4};
%caret single direction
fcolor = [ccs_vistool '/fimages/Purple-Red_Seg10_Caret.png'];
cmap_caret = ccs_mkcolormap(fcolor);
cmap_icc = cmap_caret; cmap_icc(1,:) = 0.5;
%from Ting's NI paper
JEcmaps = load('colormaps/colormap.mat');
cmap_spt = JEcmaps.cmap.spectram2(end:-1:1,:);
cmap_spt(127:129,:) = 0.5;

%% read the subject list
subjects_list = [proj_dir '/info/subjects_3t7t.xlsx'];
[status,sheets,xlFormat] = xlsfinfo(subjects_list);
[num,txt,raw] = xlsread(subjects_list,sheets{1});
load([proj_dir '/info/subidx_3t7t_rest4sess.mat'])
subjects_id = num(idx_subs_rest_3t7t,1); 
nsubs = numel(subjects_id);

%% paired t-tests on ALFF: 7T versus 3T
metrics_dir = [proj_dir '/metrics'];
%load masks
fmap = [metrics_dir '/alff_cc_slow4_mask_3t.mat'];
alff_mask_3t = load(fmap);
fmap = [metrics_dir '/alff_cc_slow4_mask_7t.mat'];
alff_mask_7t = load(fmap);
alff_mask = alff_mask_3t.alff_mask.*alff_mask_7t.alff_mask;
%load zmaps
fmap = [metrics_dir '/alff_cc_slow4_zmap_3t.mat'];
tmpalff_3t = load(fmap);
alff_3t_day1 = squeeze(mean(tmpalff_3t.alff_cc_zmap(alff_mask==1,:,1:2),3));
alff_3t_day2 = squeeze(mean(tmpalff_3t.alff_cc_zmap(alff_mask==1,:,3:4),3));
fmap = [metrics_dir '/alff_cc_slow4_zmap_7t.mat'];
tmpalff_7t = load(fmap);
alff_7t_day1 = squeeze(mean(tmpalff_7t.alff_cc_zmap(alff_mask==1,:,1:2),3));
alff_7t_day2 = squeeze(mean(tmpalff_7t.alff_cc_zmap(alff_mask==1,:,3:4),3));
%statistic tests
numtests = nnz(alff_mask_3t.alff_mask);
p = ones(numtests,2);
t = zeros(numtests,2);
for vtxid=1:numtests
    %day1
    tmp3t = alff_3t_day1(vtxid,:);
    tmp7t = alff_7t_day1(vtxid,:);
    [~,p(vtxid,1),~,tmpstats] = ttest(tmp7t,tmp3t);
    t(vtxid,1) = tmpstats.tstat;
    %day2
    tmp3t = alff_3t_day2(vtxid,:);
    tmp7t = alff_7t_day2(vtxid,:);
    [~,p(vtxid,2),~,tmpstats] = ttest(tmp7t,tmp3t);
    t(vtxid,2) = tmpstats.tstat;
end
pBFR_thr = 0.05/numtests;
for dayID=1:2
    pBFR = p(:,dayID)*numtests;
    sigBFR = -sign(t(:,dayID)).*log10(pBFR);
    sigBFR(p(:,dayID)>pBFR_thr) = 0;
    %ready for render
    sigmap = zeros(nVertices_lh+nVertices_rh,1);
    sigmap(alff_mask==1) = sigBFR;
    cmax = max(abs(sigmap));
    cmap_range = [-cmax cmax];
    cmap_tmp = cmap_cw;
    %render lh
    tmap_lh = sigmap(1:nVertices_lh);
    fprefix = ['alff_cc_slow4_7Tvs3T.day' num2str(dayID)];
    figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
    SurfStatView(tmap_lh, surfConte69_lh, ' ', 'white', 'true'); 
    colormap(cmap_tmp); SurfStatColLim(cmap_range);
    set(gcf, 'PaperPositionMode', 'auto');            
    figout = [fig_dir '/' fprefix '.lh.png'];
    print('-dpng', '-r300', figout); close
    %render rh
    tmap_rh = sigmap((1+nVertices_lh):end);
    figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
    SurfStatView(tmap_rh, surfConte69_rh, ' ', 'white', 'true'); 
    colormap(cmap_tmp); SurfStatColLim(cmap_range);
    set(gcf, 'PaperPositionMode', 'auto');            
    figout = [fig_dir '/' fprefix '.rh.png'];
    print('-dpng', '-r300', figout); close
    %save sig maps
    fmat = [proj_dir '/stats/' fprefix '.mat'];
    save(fmat,'tmap_lh','tmap_rh')
end

%% paired t-tests on ReHo: 7T versus 3T
metrics_dir = [proj_dir '/metrics'];
%load masks
fmap = [metrics_dir '/reho_cc_slow4_mask_3t.mat'];
reho_mask_3t = load(fmap);
fmap = [metrics_dir '/reho_cc_slow4_mask_7t.mat'];
reho_mask_7t = load(fmap);
reho_mask = reho_mask_3t.reho_mask.*reho_mask_7t.reho_mask;
%load zmaps
fmap = [metrics_dir '/reho_cc_slow4_zmap_3t.mat'];
tmpreho_3t = load(fmap);
reho_3t_day1 = squeeze(mean(tmpreho_3t.reho_cc_zmap(reho_mask==1,:,1:2),3));
reho_3t_day2 = squeeze(mean(tmpreho_3t.reho_cc_zmap(reho_mask==1,:,3:4),3));
fmap = [metrics_dir '/reho_cc_slow4_zmap_7t.mat'];
tmpreho_7t = load(fmap);
reho_7t_day1 = squeeze(mean(tmpreho_7t.reho_cc_zmap(reho_mask==1,:,1:2),3));
reho_7t_day2 = squeeze(mean(tmpreho_7t.reho_cc_zmap(reho_mask==1,:,3:4),3));
%statistic tests
numtests = nnz(reho_mask_3t.reho_mask);
p = ones(numtests,2);
t = zeros(numtests,2);
for vtxid=1:numtests
    %day1
    tmp3t = reho_3t_day1(vtxid,:);
    tmp7t = reho_7t_day1(vtxid,:);
    [~,p(vtxid,1),~,tmpstats] = ttest(tmp7t,tmp3t);
    t(vtxid,1) = tmpstats.tstat;
    %day2
    tmp3t = reho_3t_day2(vtxid,:);
    tmp7t = reho_7t_day2(vtxid,:);
    [~,p(vtxid,2),~,tmpstats] = ttest(tmp7t,tmp3t);
    t(vtxid,2) = tmpstats.tstat;
end
pBFR_thr = 0.05/numtests;
for dayID=1:2
    pBFR = p(:,dayID)*numtests;
    sigBFR = -sign(t(:,dayID)).*log10(pBFR);
    sigBFR(p(:,dayID)>pBFR_thr) = 0;
    %ready for render
    sigmap = zeros(nVertices_lh+nVertices_rh,1);
    sigmap(reho_mask==1) = sigBFR;
    cmax = max(abs(sigmap));
    cmap_range = [-cmax cmax];
    cmap_tmp = cmap_cw;
    fprefix = ['reho_cc_slow4_7Tvs3T.day' num2str(dayID)];
    %render lh
    tmap_lh = sigmap(1:nVertices_lh);
    figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
    SurfStatView(tmap_lh, surfConte69_lh, ' ', 'white', 'true'); 
    colormap(cmap_tmp); SurfStatColLim(cmap_range);
    set(gcf, 'PaperPositionMode', 'auto');            
    figout = [fig_dir '/' fprefix '.lh.png'];
    print('-dpng', '-r300', figout); close
    %render rh
    tmap_rh = sigmap((1+nVertices_lh):end);
    figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
    SurfStatView(tmap_rh, surfConte69_rh, ' ', 'white', 'true'); 
    colormap(cmap_tmp); SurfStatColLim(cmap_range);
    set(gcf, 'PaperPositionMode', 'auto');            
    figout = [fig_dir '/' fprefix '.rh.png'];
    print('-dpng', '-r300', figout); close
    %save sig maps
    fmat = [proj_dir '/stats/' fprefix '.mat'];
    save(fmat,'tmap_lh','tmap_rh')
end

%% paired t-tests on VMHC: 7T versus 3T
metrics_dir = [proj_dir '/metrics'];
%load masks
fmap = [metrics_dir '/vmhc_cc_slow4_mask_3t.mat'];
vmhc_mask_3t = load(fmap);
fmap = [metrics_dir '/vmhc_cc_slow4_mask_7t.mat'];
vmhc_mask_7t = load(fmap);
vmhc_mask = vmhc_mask_3t.vmhc_mask.*vmhc_mask_7t.vmhc_mask;
%load zmaps
fmap = [metrics_dir '/vmhc_cc_slow4_zmap_3t.mat'];
tmpvmhc_3t = load(fmap);
vmhc_3t_day1 = squeeze(mean(tmpvmhc_3t.vmhc_cc_zmap(vmhc_mask==1,:,1:2),3));
vmhc_3t_day2 = squeeze(mean(tmpvmhc_3t.vmhc_cc_zmap(vmhc_mask==1,:,3:4),3));
fmap = [metrics_dir '/vmhc_cc_slow4_zmap_7t.mat'];
tmpvmhc_7t = load(fmap);
vmhc_7t_day1 = squeeze(mean(tmpvmhc_7t.vmhc_cc_zmap(vmhc_mask==1,:,1:2),3));
vmhc_7t_day2 = squeeze(mean(tmpvmhc_7t.vmhc_cc_zmap(vmhc_mask==1,:,3:4),3));
%statistic tests
numtests = nnz(vmhc_mask_3t.vmhc_mask);
p = ones(numtests,2);
t = zeros(numtests,2);
for vtxid=1:numtests
    %day1
    tmp3t = vmhc_3t_day1(vtxid,:);
    tmp7t = vmhc_7t_day1(vtxid,:);
    [~,p(vtxid,1),~,tmpstats] = ttest(tmp7t,tmp3t);
    t(vtxid,1) = tmpstats.tstat;
    %day2
    tmp3t = vmhc_3t_day2(vtxid,:);
    tmp7t = vmhc_7t_day2(vtxid,:);
    [~,p(vtxid,2),~,tmpstats] = ttest(tmp7t,tmp3t);
    t(vtxid,2) = tmpstats.tstat;
end
pBFR_thr = 0.05/numtests;
for dayID=1:2
    pBFR = p(:,dayID)*numtests;
    sigBFR = -sign(t(:,dayID)).*log10(pBFR);
    sigBFR(p(:,dayID)>pBFR_thr) = 0;
    fprefix = ['vmhc_cc_slow4_7Tvs3T.day' num2str(dayID)];
    %ready for render
    sigmap = zeros(nVertices_lh+nVertices_rh,1);
    sigmap(vmhc_mask==1) = sigBFR;
    sigmap(abs(sigmap)>20) = 0;
    cmax = max(abs(sigmap));
    cmap_range = [-cmax cmax];
    cmap_tmp = cmap_cw;
    %render lh
    tmap_lh = sigmap(1:nVertices_lh);
    figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
    SurfStatView(tmap_lh, surfConte69_lh, ' ', 'white', 'true'); 
    colormap(cmap_tmp); SurfStatColLim(cmap_range);
    set(gcf, 'PaperPositionMode', 'auto');            
    figout = [fig_dir '/' fprefix '.lh.png'];
    print('-dpng', '-r300', figout); close
    %render rh
    tmap_rh = sigmap((1+nVertices_lh):end);
    figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
    SurfStatView(tmap_rh, surfConte69_rh, ' ', 'white', 'true'); 
    colormap(cmap_tmp); SurfStatColLim(cmap_range);
    set(gcf, 'PaperPositionMode', 'auto');            
    figout = [fig_dir '/' fprefix '.rh.png'];
    print('-dpng', '-r300', figout); close
    %save sig maps
    fmat = [proj_dir '/stats/' fprefix '.mat'];
    save(fmat,'tmap_lh','tmap_rh')
end
