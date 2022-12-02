function [idx_surfroi,coord_seed] = ccs_core_seedsurfroi(seed_coord, ...
    ConteAtlas_dir, SurfGraph, sizeROI)
%CCS_CORE_SEEDSURFROI Locate the vertices of the region of interests (ROI)
%   seeded by the given coordination (x,y,z) in MNI space
% Input:   
%   seed_coord - Coordinates for the seed vertex
%   ConteAtlas_dir - Directory for the Conte Surface Atlas
%   SurfGraph - Graphs contain neighbour information
%   sizeROI - steps of neighbours of the ROI (1-8)
% Output:
%   idx_surfroi - vertices' indices within the ROI (two-step neighbours)
%   coord_seed - seed coordinates of the seed ROI
% Credit:
%   Xi-Nian Zuo (Beijing Normal University)

if nargin<4
    sizeROI = 1;
end
% left hemisphere
Conte32k_lh = gifti([ConteAtlas_dir ...
    '/Conte69.L.midthickness.32k_fs_LR.surf.gii']);
nVertices_lh = size(Conte32k_lh.vertices,1);
Coordinates_lh = Conte32k_lh.vertices'; 
% right hemisphere
Conte32k_rh = gifti([ConteAtlas_dir ...
    '/Conte69.R.midthickness.32k_fs_LR.surf.gii']);
nVertices_rh = size(Conte32k_rh.vertices,1);
Coordinates_rh = Conte32k_rh.vertices';
%brain surface graphs
surfgraph = load(SurfGraph);

%locate the vertex nearest to the seed
if size(seed_coord,1)~=3
    seed_coord = seed_coord';
end
seed_dist_lh = vecnorm(Coordinates_lh - repmat(seed_coord,1,nVertices_lh));
[dist_min_lh, idx_min_lh] = min(seed_dist_lh);
seed_dist_rh = vecnorm(Coordinates_rh - repmat(seed_coord,1,nVertices_rh));
[dist_min_rh, idx_min_rh] = min(seed_dist_rh);
%find the vertices in the ROI
if dist_min_lh < dist_min_rh
    coord_seed = Coordinates_lh(:,idx_min_lh);
    idx_surfroi = surfgraph.lh_nbrs{idx_min_lh,sizeROI};
else
    coord_seed = Coordinates_rh(:,idx_min_rh);
    idx_surfroi = surfgraph.rh_nbrs{idx_min_rh,sizeROI};
end
