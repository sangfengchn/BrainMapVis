%% setup major directories
clear; clc
disk_dir = '/Volumes/LaCie';
hcp_dir = [disk_dir '/HCP'];
parc_dir = [disk_dir '/Schaefer2018_LocalGlobal'];
proj_dir = [hcp_dir '/zprojects/reliability'];
fig_dir = [proj_dir '/group/figures'];

%% setup major toolboxes
fs_home = '/Applications/freesurfer/freesurfer60';
%connectome computation system
ccs_dir = '/Users/mac/Projects/CCS';
ccs_matlab = [ccs_dir '/matlab'];
ccs_vistool = [ccs_dir '/vistool'];
%hcp workbench
hcpwkbc_dir = [hcp_dir '/workbench'];
wbcmd = [hcpwkbc_dir '/macosx64_apps/' ...
    'wb_command.app/Contents/MacOS/wb_command'];
ConteAtlas_dir = [hcp_dir '/zprojects/surfmodels/32k_ConteAtlas_v2'];
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
Conte32k_rh = gifti([ConteAtlas_dir ...
    '/Conte69.R.midthickness.32k_fs_LR.surf.gii']);

%% build lh graph: up to 2*8 = 16mm distance
FV.vertices = Conte32k_lh.vertices ; FV.faces = Conte32k_lh.faces; 
nVertices_lh = size(Conte32k_lh.vertices,1); edge_lh = mesh_adjacency(FV); 
%searching neighbours (nrbs)
lh_nbrs = cell(nVertices_lh,8);
for nstep=1:8
    if nstep > 1
        tmpedge = tmpedge*edge_lh;
    else
        tmpedge = edge_lh;
    end
    for k=1:nVertices_lh
        if mod(k,1000)==0
            disp(['Counting neighbours for ' num2str(nstep) '-step: ' ...
                num2str(round(100*k/nVertices_lh)) ...
                ' percentage vertices completed...'])
        end
        lh_nbrs{k,nstep} = find(tmpedge(k,:)>0);
    end
end

%% build rh graph: up to 2*8 = 16mm distance
FV.vertices = Conte32k_rh.vertices ; FV.faces = Conte32k_rh.faces; 
nVertices_rh = size(Conte32k_rh.vertices,1); edge_rh = mesh_adjacency(FV); 
%searching neighbours (nrbs)
rh_nbrs = cell(nVertices_rh,8);
for nstep=1:8
    if nstep > 1
        tmpedge = tmpedge*edge_rh;
    else
        tmpedge = edge_rh;
    end
    for k=1:nVertices_rh
        if mod(k,1000)==0
            disp(['Counting neighbours for ' num2str(nstep) '-step: ' ...
                num2str(round(100*k/nVertices_rh)) ...
                ' percentage vertices completed...'])
        end
        rh_nbrs{k,nstep} = find(tmpedge(k,:)>0);
    end
end

%% render connection matrix
subplot(121), spy(edge_lh); subplot(122), spy(edge_rh)
edge = [edge_lh sparse(nVertices_lh, nVertices_rh); ...
    sparse(nVertices_lh, nVertices_rh) edge_rh];
figure, spy(edge)

%% save neighbour lists
save('Conte69_32k_surfgraph.mat', 'lh_nbrs', 'rh_nbrs')
