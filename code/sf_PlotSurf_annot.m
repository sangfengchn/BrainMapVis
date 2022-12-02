clc, clear; close all;
% add path
addpath('toolbox/gifti-master');

% load surface and annot file
lh_Surf = gifti('fsaverage_lh_inflatedwhite.surf.gii');
[vertics, label, colortable] = read_annotation('lh.aparc_HCP_MMP1.freesurfer.annot');
% colortable.table = uint64(colortable.table);

lh_Atlas.cdata = label;
lh_Atlas.labels.key = colortable.table(:, 5)'; % red, green, blue, alpha, index                        
lh_Atlas.labels.rgba = colortable.table(:, 1:4);
lh_Atlas.labels.rgba(:, 1:3) = lh_Atlas.labels.rgba(:, 1:3) / 255;
lh_Atlas.labels.rgba(:, 4) = 1 - lh_Atlas.labels.rgba(:, 4);

lh_SurfView.tri = lh_Surf.faces;
lh_SurfView.coord = lh_Surf.vertices';

% load color from .label.gii and transform the fit format
colors = zeros(numel(lh_Atlas.cdata), 4); % rgba
for i = lh_Atlas.labels.key
    tmpIdx = find(lh_Atlas.labels.key == i);
    tmpIdx_Color = lh_Atlas.labels.rgba(tmpIdx, :);
    tmpIdxMask_cdata = lh_Atlas.cdata == i;
    colors(tmpIdxMask_cdata, :) = repmat(tmpIdx_Color, [numel(find(tmpIdxMask_cdata)), 1]);
end

% plotting and save figure
figure('Units', 'pixel', 'Position', [100 100 800 800]); 
patch_Surf = trisurf(lh_SurfView.tri, ...
    lh_SurfView.coord(1,:), ...
    lh_SurfView.coord(2,:), ...
    lh_SurfView.coord(3,:));
patch_Surf.EdgeColor = 'interp';
patch_Surf.FaceColor = 'interp';
patch_Surf.FaceVertexCData = repmat([207, 210, 207] / 255, [numel(lh_SurfView.coord(1, :)), 1]);

view(-90, 0);
daspect([1 1 1]);
axis vis3d off;
lighting flat;
material dull;
shading flat;
camlight('headlight', 'infinite');
axis tight;
hold on; % overlap tow patches by using hold on/off

patch_Atlas = trisurf(lh_SurfView.tri, ...
    lh_SurfView.coord(1,:), ...
    lh_SurfView.coord(2,:), ...
    lh_SurfView.coord(3,:));
% set vertex color
patch_Atlas.EdgeColor = 'interp';
patch_Atlas.FaceColor = 'interp';
patch_Atlas.FaceVertexCData = colors(:, 1:3);
% set vertex alpha
patch_Atlas.EdgeAlpha = 'interp';
patch_Atlas.FaceAlpha = 'interp';
patch_Atlas.AlphaDataMapping = 'none';
patch_Atlas.FaceVertexAlphaData = colors(:, 4);
material dull;
shading flat;
hold off;
% alpha(patch_Atlas, 0.5);
% exportgraphics(gcf, 'left_area.png', 'Resolution', 1000, 'BackgroundColor', 'none');
% close;