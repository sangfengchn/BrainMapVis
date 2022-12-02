clear all; clc
hcp_dir = '/Volumes/LaCie/HCP';
ccs_dir = '/Users/mac/Projects/CCS';
cbig_dir = '/Volumes/DJICopilot/CBIG-master';
spm_dir = [ccs_dir '/extool/spm12'];
ccs_matlab = [ccs_dir '/matlab'];
ccs_vistool = [ccs_dir '/vistool'];
fs_home = '/Applications/freesurfer/freesurfer53';
cifti_matlab = [ccs_dir '/extool/cifti-matlab'];
hcpwkbc_dir = [ccs_dir '/extool/hcpworkbench'];
wbcmd = [hcpwkbc_dir ...
    '/macosx64_apps/wb_command.app/Contents/MacOS/wb_command'];
atlas_dir = [hcpwkbc_dir '/resources/32k_ConteAtlas_v2'];
work_dir = '/Volumes/LaCie/HCP';
%add the paths to matlab function
addpath(genpath(ccs_matlab)) %ccs matlab scripts
addpath(genpath(ccs_vistool)) %ccs matlab scripts
addpath(genpath(cifti_matlab)) %cifti paths
addpath(genpath([fs_home '/matlab'])) %freesurfer matlab scripts

%% resting satte FMRI names
rest_name = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', ...
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL'}; 
nscans = numel(rest_name);

%% load symmetric surfaces - make Conte69 surface structure
fSURF = [atlas_dir '/Conte69.L.inflated.32k_fs_LR.surf.gii'];
lh_inflated = gifti(fSURF); 
numVertices_lh = size(lh_inflated.vertices,1); %lh
surfConte69_lh.tri = lh_inflated.faces;
surfConte69_lh.coord = lh_inflated.vertices'; 

fSURF = [atlas_dir '/Conte69.R.inflated.32k_fs_LR.surf.gii'];
rh_inflated = gifti(fSURF); 
numVertices_rh = size(rh_inflated.vertices,1); %rh
surfConte69_rh.tri = rh_inflated.faces;
surfConte69_rh.coord = rh_inflated.vertices'; 

%% load aparc parcellation
subject = '103818';
subj_dir = [work_dir '/TRT/' subject];
faparc = [subj_dir '/MNINonLinear/fsaverage_LR32k/' ...
    subject '.aparc.32k_fs_LR.dlabel.nii'];
aparc = ft_read_cifti(faparc);
aparc_lh = aparc.x103818_aparc(aparc.brainstructure==1);
aparc_lh(isnan(aparc_lh)) = 0;
num_parcs_lh = max(aparc_lh);
aparc_rh = aparc.x103818_aparc(aparc.brainstructure==2);
aparc_rh(isnan(aparc_rh)) = 0;
num_parcs_rh = max(aparc_rh) - max(aparc_lh);
%read annotation from FS5
fannot_aparc_lh = [fs_home '/subjects/fsaverage/label/lh.aparc.annot'];
[vertices_lh,labels_lh,colortable_lh] = read_annotation(fannot_aparc_lh);
fannot_aparc_rh = [fs_home '/subjects/fsaverage/label/rh.aparc.annot'];
[vertices_rh,labels_rh,colortable_rh] = read_annotation(fannot_aparc_rh);
%parcel labels and colors - lh
colortable_parcs_lh = colortable_lh.table(2:end,1:3);
names_parcs_lh = cell(num_parcs_lh,2);
for nameID=1:num_parcs_lh
    tmpname = colortable_lh.struct_names{nameID+1};
    names_parcs_lh{nameID,1} = tmpname;
    names_parcs_lh{nameID,2} = ['L.' upper(tmpname(1:3:end))];
end
%parcel labels and colors - rh
colortable_parcs_rh = colortable_rh.table(2:end,1:3);
names_parcs_rh = cell(num_parcs_rh,2);
for nameID=1:num_parcs_rh
    tmpname = colortable_rh.struct_names{nameID+1};
    names_parcs_rh{nameID,1} = tmpname;
    names_parcs_rh{nameID,2} = ['R.' upper(tmpname(1:3:end))];
end
%colormap
cmap_fsLR = [[128 128 128]; colortable_parcs_lh; colortable_parcs_rh]/256;
%double check with render surfaces - lh
figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
SurfStatView(aparc_lh, surfConte69_lh, ' ', 'white', 'true'); 
colormap(cmap_fsLR); SurfStatColLim([0 70]);
set(gcf, 'PaperPositionMode', 'auto');
%double check with render surfaces - rh
figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
SurfStatView(aparc_rh, surfConte69_rh, ' ', 'white', 'true'); 
colormap(cmap_fsLR); SurfStatColLim([0 70]);
set(gcf, 'PaperPositionMode', 'auto');
%save the updates of fsLR parcellation
fdk2006 = [ccs_dir '/parcellation/DK2006_32k_fsLR.mat'];
save(fdk2006, 'aparc_lh', 'colortable_parcs_lh', 'names_parcs_lh', ...
    'aparc_rh', 'colortable_parcs_rh', 'names_parcs_rh');

%% Build individual brain graph
load([ccs_dir '/parcellation/DK2006_32k_fsLR.mat']);
num_parcs = num_parcs_lh + num_parcs_rh;
corr_mat = zeros(num_parcs,num_parcs);
for sid=1:nscans
    %load bold time series
    func_dir_name = ['MNINonLinear/Results/' rest_name{sid}];
    disp(['Loading rfMRI for subject ' subject ': ' rest_name{sid} ' ...'])
    func_dir = [work_dir '/TRT/' subject '/' func_dir_name];
    fbold = [func_dir '/' rest_name{sid} ...
        '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
    boldts = ft_read_cifti(fbold);
    boldts_lh = boldts.dtseries(boldts.brainstructure==1,:);
    boldts_rh = boldts.dtseries(boldts.brainstructure==2,:);
    tslength = size(boldts_lh,2);
    clear boldts
    %compute the representative ts - lh
    tmp_parcts_lh = zeros(tslength,num_parcs_lh);
    for parcID=1:num_parcs_lh
        tmpidx = find(aparc_lh==parcID);
        tmpts = boldts_lh(tmpidx,:);
        tmp_parcts_lh(:,parcID) = mean(tmpts);
    end
    %compute the representative ts - rh
    tmp_parcts_rh = zeros(tslength,num_parcs_rh);
    for parcID=1:num_parcs_rh
        tmpidx = find(aparc_rh==(parcID+35));
        tmpts = boldts_rh(tmpidx,:);
        tmp_parcts_rh(:,parcID) = mean(tmpts);
    end
    %compute the correlation matrix
    tmp_parcts = [tmp_parcts_lh tmp_parcts_rh];
    tmpcorr = ccshcp_core_fastcorr(tmp_parcts,tmp_parcts);
    corr_mat = corr_mat + atanh(tmpcorr);
end
corr_mat = tanh(corr_mat/4);
%cut off the CC parcel (not part of the cortex)

%% Estimate adjacency matrix - binary
edge_density = 0.1; centers = -1:0.001:1;
number_edges = zeros(size(centers));
tmpTRI = tril(corr_mat,-1);
tmpNNZ = tmpTRI(tmpTRI~=0);
number_edges = number_edges + hist(tmpNNZ(:), centers);
cdf_edges = cumsum(number_edges)/sum(number_edges);
idx_corr_thr = find(cdf_edges >= (1-edge_density));
corr_thr = centers(idx_corr_thr(1));
corr_mat_thr = corr_mat;
corr_mat_thr(corr_mat<corr_thr) = 0;
corr_mat_thr(corr_mat>=corr_thr) = 1;
adj_mat = sparse(corr_mat_thr) - speye(num_parcs);
spy(adj_mat); nnz(full(sum(adj_mat)))


%% Render a biograph
node_names = cell(num_parcs,1);
node_desc = cell(num_parcs,1);
node_rgbs = zeros(num_parcs,3);
for nodeID=1:num_parcs_lh
    node_names{nodeID} = names_parcs_lh{nodeID,2};
    node_desc{nodeID} = names_parcs_lh{nodeID,1};
    node_rgbs(nodeID,:) = colortable_parcs_lh(nodeID,:)/256;
    node_names{nodeID+num_parcs_lh} = names_parcs_rh{nodeID,2};
    node_desc{nodeID+num_parcs_lh} = names_parcs_rh{nodeID,1};
    node_rgbs(nodeID+num_parcs_lh,:) = colortable_parcs_rh(nodeID,:)/256;
end
%left hemisphere
adj_mat_lh = adj_mat(1:num_parcs_lh,1:num_parcs_lh);
node_names_lh = node_names(1:num_parcs_lh);
node_desc_lh = node_desc(1:num_parcs_lh);
gObj_lh = biograph(adj_mat_lh);
for nodeID=1:num_parcs_lh
    gObj_lh.nodes(nodeID).Size = [80 60];
    gObj_lh.nodes(nodeID).TextColor = [1 1 1];
    gObj_lh.nodes(nodeID).Color = node_rgbs(nodeID,:);
    gObj_lh.nodes(nodeID).Label = node_names_lh{nodeID};
    gObj_lh.nodes(nodeID).Description = node_desc_lh{nodeID};
end
lh_gObjFig = view(gObj_lh);
%right hemisphere
adj_mat_rh = adj_mat((1+num_parcs_lh):end,(1+num_parcs_lh):end);
node_names_rh = node_names((1+num_parcs_lh):end);
node_desc_rh = node_desc((1+num_parcs_lh):end);
gObj_rh = biograph(adj_mat_rh);
for nodeID=1:num_parcs_rh
    gObj_rh.nodes(nodeID).Size = [80 60];
    gObj_rh.nodes(nodeID).TextColor = [1 1 1];
    gObj_rh.nodes(nodeID).Color = node_rgbs(nodeID+num_parcs_lh,:);
    gObj_rh.nodes(nodeID).Label = node_names_rh{nodeID};
    gObj_rh.nodes(nodeID).Description = node_desc_rh{nodeID};
end
rh_gObjFig = view(gObj_rh);
%full cortex
gObj = biograph(adj_mat);
for nodeID=1:num_parcs
    gObj.nodes(nodeID).Size = [80 60];
    gObj.nodes(nodeID).TextColor = [1 1 1];
    gObj.nodes(nodeID).Color = node_rgbs(nodeID,:);
    gObj.nodes(nodeID).Label = node_names{nodeID};
    gObj.nodes(nodeID).Description = node_desc{nodeID};
end
gObjFig = view(gObj);