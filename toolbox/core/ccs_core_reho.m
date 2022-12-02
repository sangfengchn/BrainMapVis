function [ reho ] = ccs_core_reho(tsmat, nbrs)
%CCS_CORE_REHO MATLAB version of ReHo Computation
%   Detailed explanation goes here
%   INPUTS
%       tsmat - a TxN matrix of time series
%       nbrs - neighbour cells
if isempty(nbrs)
    reho = ccs_core_kcc(tsmat);
else
    [~, nsp] = size(tsmat); %num of time points, num of spatial points
    reho = zeros(nsp,1);
    for ii=1:nsp        
        tmp_ts = squeeze(tsmat(:,ii));
        if std(tmp_ts) > 0
            tmp_nbrs = nbrs{ii};
            nbrs_ts = squeeze(tsmat(:,tmp_nbrs));% total neighbor vertices
            ts = [tmp_ts nbrs_ts(:,std(nbrs_ts,1,1)>0)]; 
            reho(ii) = ccs_core_kcc(ts);
        end
   end
end

