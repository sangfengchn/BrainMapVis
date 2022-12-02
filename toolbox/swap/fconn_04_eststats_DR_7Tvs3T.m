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
JEcmaps = load(['colormaps/colormap.mat']);
cmap_spt = JEcmaps.cmap.spectram2(end:-1:1,:);
cmap_spt(127:129,:) = 0.5;

%% read the subject list
subjects_list = [proj_dir '/info/subjects_3t7t.xlsx'];
[status,sheets,xlFormat] = xlsfinfo(subjects_list);
[num,txt,raw] = xlsread(subjects_list,sheets{1});
load([proj_dir '/info/subidx_3t7t_rest4sess.mat'])
subjects_id = num(idx_subs_rest_3t7t,1); 
nsubs = numel(subjects_id);
numRSNs = 7; numSess = 4;

%% paired t-tests on DR maps: 7T versus 3T
metrics_dir = [proj_dir '/metrics'];
%load masks
fmap = [metrics_dir '/dualreg_cc_slow4_mask_3t.mat'];
dualreg_mask_3t = load(fmap);
fmap = [metrics_dir '/dualreg_cc_slow4_mask_7t.mat'];
dualreg_mask_7t = load(fmap);
dualreg_mask = dualreg_mask_3t.dualreg_mask.*dualreg_mask_7t.dualreg_mask;
%load zmaps
fmap = [metrics_dir '/dualreg_cc_slow4_zmap_3t.mat'];
tmpDR_3t = load(fmap);
dualreg_3t = squeeze(mean(tmpDR_3t.dualreg_cc_zmap(dualreg_mask==1,:,:,:),4));
fmap = [metrics_dir '/dualreg_cc_slow4_zmap_7t.mat'];
tmpDR_7t = load(fmap);
dualreg_7t = squeeze(mean(tmpDR_7t.dualreg_cc_zmap(dualreg_mask==1,:,:,:),4));
%statistic tests
numtests = nnz(dualreg_mask);
p = ones(numtests,numRSNs);
t = zeros(numtests,numRSNs);
for idxRSN=1:numRSNs
    dr3T = squeeze(dualreg_3t(:,idxRSN,:));
    dr7T = squeeze(dualreg_7t(:,idxRSN,:));
    for vtxid=1:numtests
        if ~mod(vtxid,500) 
            disp(['Network' num2str(idxRSN) ': ', num2str(vtxid/numtests*100) ...
            ' percent vertices are processed ...'])
        end
        tmp3t = dr3T(vtxid,:);
        tmp7t = dr7T(vtxid,:);
        [~,p(vtxid,idxRSN),~,tmpstats] = ttest(tmp7t,tmp3t);
        t(vtxid,idxRSN) = tmpstats.tstat;
    end
end
pBFR_thr = 0.05/numtests;
pBFR = p*numtests;
sigBFR = -sign(t).*log10(pBFR);
sigBFR(p>pBFR_thr) = 0;
sigmap = zeros(nVertices_lh+nVertices_rh,7);
sigmap(dualreg_mask==1,:) = sigBFR;
sigmap_lh = sigmap(1:nVertices_lh,:);
sigmap_rh = sigmap((1+nVertices_lh):end,:);
%save stats maps
fprefix = 'dualreg_cc_slow4_7Tvs3T';
fmat = [proj_dir '/stats/' fprefix '.mat'];
save(fmat,'sigmap_lh','sigmap_rh','dualreg_mask','p','t')

%% Rendering the statistical maps
for idxRSN=1:numRSNs
    cmax = max(abs(sigmap(:,idxRSN)));
    cmap_range = [-cmax cmax];
    cmap_tmp = cmap_cw;
    %render lh
    tmap_lh = sigmap_lh(:,idxRSN);
    figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
    SurfStatView(tmap_lh, surfConte69_lh, ' ', 'white', 'true'); 
    colormap(cmap_tmp); SurfStatColLim(cmap_range);
    set(gcf, 'PaperPositionMode', 'auto');            
    figout = [fig_dir '/dualreg_network' num2str(idxRSN) '_slow4_7Tvs3T.lh.png'];
    print('-dpng', '-r300', figout); close
    %render rh
    tmap_rh = sigmap_rh(:,idxRSN);
    figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
    SurfStatView(tmap_rh, surfConte69_rh, ' ', 'white', 'true'); 
    colormap(cmap_tmp); SurfStatColLim(cmap_range);
    set(gcf, 'PaperPositionMode', 'auto');            
    figout = [fig_dir '/dualreg_network' num2str(idxRSN) '_slow4_7Tvs3T.rh.png'];
    print('-dpng', '-r300', figout); close
end