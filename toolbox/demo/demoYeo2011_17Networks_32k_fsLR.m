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

%% load Yeo17-networks parcellation
fyeo2011 = [ccs_dir '/parcellation/Yeo114_32k_fsLR.mat'];
yeo114 = load(fyeo2011);
labels_fsLR_lh = yeo114.lh_labels_fslr;
nVertices_lh = numel(labels_fsLR_lh);
labels_fsLR_rh = yeo114.rh_labels_fslr;
nVertices_rh = numel(labels_fsLR_rh);
%read annotation from FS5
fannot_yeo17_lh = [ccs_dir '/parcellation/ParcelsYeo2011/fsaverage' ...
    '/label/lh.Yeo2011_17Networks_N1000.split_components.annot'];
[vertices_lh,labels_lh,colortable_lh] = read_annotation(fannot_yeo17_lh);
fannot_yeo17_rh = [ccs_dir '/parcellation/ParcelsYeo2011/fsaverage' ...
    '/label/rh.Yeo2011_17Networks_N1000.split_components.annot'];
[vertices_rh,labels_rh,colortable_rh] = read_annotation(fannot_yeo17_rh);
%parcel labels and colors - lh
colortable_fsLR_lh = colortable_lh.table(1:58,1:3);
names_fsLR_lh = cell(58,2);
%colormap
cmap_fsLR = [colortable_fsLR_lh; colortable_fsLR_rh(2:end,:)]/255;
%cmap_fsLR(1,:) = 0.5;
for nameID=1:58
    tmpname = colortable_lh.struct_names{nameID};
    idxUDL = strfind(tmpname,'_');
    if numel(idxUDL) > 2
        idxNetwork = (idxUDL(2)+1):(idxUDL(3)-1);
    else
        idxNetwork = (idxUDL(2)+1):length(tmpname);
    end
    names_fsLR_lh{nameID,1} = tmpname(idxNetwork);
    if numel(idxUDL) > 2
        idxRegion = (idxUDL(3)+1):length(tmpname);
        names_fsLR_lh{nameID,2} = tmpname(idxRegion);
    else
        names_fsLR_lh{nameID,2} = names_fsLR_lh{nameID,1};
    end
end
%parcel labels and colors - rh
colortable_fsLR_rh = colortable_rh.table([1 59:end],1:3);
names_fsLR_rh = cell(58,2);
for nameID=1:58
    if nameID > 1
        tmpname = colortable_rh.struct_names{nameID+57};
    else
        tmpname = colortable_rh.struct_names{nameID};
    end
    idxUDL = strfind(tmpname,'_');
    if numel(idxUDL) > 2
        idxNetwork = (idxUDL(2)+1):(idxUDL(3)-1);
    else
        idxNetwork = (idxUDL(2)+1):length(tmpname);
    end
    names_fsLR_rh{nameID,1} = tmpname(idxNetwork);
    if numel(idxUDL) > 2
        idxRegion = (idxUDL(3)+1):length(tmpname);
        names_fsLR_rh{nameID,2} = tmpname(idxRegion);
    else
        names_fsLR_rh{nameID,2} = names_fsLR_rh{nameID,1};
    end
end
%double check with render surfaces - lh
figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
SurfStatView(labels_fsLR_lh, surfConte69_lh, ' ', 'white', 'true'); 
colormap(cmap_fsLR); SurfStatColLim([0 114]);
set(gcf, 'PaperPositionMode', 'auto');
%double check with render surfaces - rh
figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
SurfStatView(labels_fsLR_rh, surfConte69_rh, ' ', 'white', 'true'); 
colormap(cmap_fsLR); SurfStatColLim([0 114]);
set(gcf, 'PaperPositionMode', 'auto');
%save the updates of fsLR parcellation
fyeo2011 = [ccs_dir '/parcellation/Yeo114_17Networks_32k_fsLR.mat'];
save(fyeo2011, 'labels_fsLR_lh', 'colortable_fsLR_lh', 'names_fsLR_lh', ...
    'labels_fsLR_rh', 'colortable_fsLR_rh', 'names_fsLR_rh');
