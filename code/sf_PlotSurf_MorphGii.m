clc; clear; close all;

% add toolbox path
addpath('toolbox/gifti-master');
addpath('toolbox/core');

% load surface file and measure file
lh_Surf = gifti('fsaverage_lh_inflatedwhite.surf.gii');
lh_Thickness = gifti('surf/s12.lh.area.resampled.sub-12500_T1w.gii');
% lh_Surf = gifti('surf/lh.central.sub-12500_T1w.gii');
% lh_Thickness = gifti('surf/lh.thickness.sub-12500_T1w.gii');


% generate the new data for visualization
lh_SurfView.tri = lh_Surf.faces;
lh_SurfView.coord = lh_Surf.vertices';

% load color map figure, set the color of NaN
cmap_caret2 = ccs_mkcolormap('toolbox/colormaps/Spectral.png');
colorNan = [207, 210, 207] / 255;

% generate colormap by cdata linearly, and set color of NaN
idx_CDataNan = isnan(lh_Thickness.cdata);
color_Vertex = sf_GetColors(cmap_caret2, lh_Thickness.cdata(idx_CDataNan == false));
colors = zeros(numel(lh_Thickness.cdata), 3);
colors(idx_CDataNan == false, :) = color_Vertex;
colors(idx_CDataNan, :) = repmat(colorNan, [numel(find(idx_CDataNan)), 1]);

% plotting and save figure
figure('Units', 'pixel', 'Position', [100 100 800 800]); 
f1 = trisurf(lh_SurfView.tri, ...
    lh_SurfView.coord(1,:), ...
    lh_SurfView.coord(2,:), ...
    lh_SurfView.coord(3,:));
f1.EdgeColor = 'interp';
f1.FaceColor = 'interp';
f1.FaceVertexCData = colors;
view(-90, 0);
daspect([1 1 1]);
axis vis3d off;
lighting flat;
material dull;
shading flat;
camlight('headlight', 'infinite');
axis tight;
exportgraphics(gcf, 'left_area.png', 'Resolution', 1000, 'BackgroundColor', 'none');
close;

% function: generate colormap
function colors = sf_GetColors(colormap, x)
xMin = min(x);
xMax = max(x);
xIdx = uint8((x - xMin) / (xMax - xMin) * 255) + 1;
colors = colormap(xIdx, :);
end