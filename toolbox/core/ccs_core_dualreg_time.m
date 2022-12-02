function [network_maps] = ccs_core_dualreg_time(...
    rfmrivol_dr, time_regressors, type)
%CCS_CORE_DUALREG_TIME Perform time dual regression of rfMRI images on 
% the predefined regressors.
%Input:
%   rfmrivol_dr -- the preprocessed resting state FMRI data (masked)
%   sp_regressors -- spatial regressiors for the first regression (cell)
%   type -- multiple regression or sigle map regression ('multiple' or
%   'single')

if nargin < 3
    type = 'multiple';
end

%get basic information
[numVertices, numTRs] = size(rfmrivol_dr);
[~, numTIs] = size(time_regressors);
time_series = zeros(numTRs,numTIs);
network_maps = zeros(numVertices, numTIs);

%dr2: temporal regression
tmpmDesignDR = [ones(numTRs,1) zscore(time_regressors)];
for vtxID=1:numVertices
    tmpY = rfmrivol_dr(vtxID,:);
    if strcmp(type, 'multiple')
        tmpB = regress(tmpY', tmpmDesignDR);
        network_maps(vtxID,:) = tmpB(2:end);
    else
        for cogID=1:numTIs
            tmpsDesignDR = [ones(numTRs,1) zscore(time_series(:,cogID))];
            tmpB = regress(tmpY', tmpsDesignDR);
            network_maps(vtxID,cogID) = tmpB(2);
        end
    end
end

end