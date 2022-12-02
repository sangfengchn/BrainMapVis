function [mask,idx_mask] = ccs_core_boldmask(boldts,tdim)
%CCS_CORE_BOLDMASK generate spatial mask for 2D timeseries data
%   boldts - bold time series matrix
%   tdim - the dimension of time
if numel(size(boldts))==2
    tmpstd = std(boldts,0,tdim);
    mask = zeros(numel(tmpstd),1);
    idx_mask = find(tmpstd > 0);
    mask(idx_mask) = 1;
else
    disp('Please ensure the input data is 2D ...')
end

end

