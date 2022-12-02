function [time_series] = ccs_core_dualreg_space(...
    rfmrivol_dr, sp_regressors, type)
%CCS_CORE_DUALREG_SPACE Perform the space dual regression of rfMRI images 
% on the predefined regressors.
%Input:
%   rfmrivol_dr -- the preprocessed resting state FMRI data (masked)
%   sp_regressors -- spatial regressiors for the first regression (cell)
%   type -- multiple regression or sigle map regression ('multiple' or
%   'single')
if nargin < 3
    type = 'single';
end

%get basic information
[~, numTRs] = size(rfmrivol_dr);
[numVertices_tmpt, numSPs] = size(sp_regressors);
time_series = zeros(numTRs,numSPs);

%dr1: spatial regression
tmpmDesignDR = [ones(numVertices_tmpt,1) sp_regressors];
for trID=1:numTRs 
    tmpY = rfmrivol_dr(:,trID);
    if strcmp(type, 'multiple')
        tmpB = regress(tmpY, tmpmDesignDR);
        time_series(trID,:) = tmpB(2:end);
    else %one-by-one
        for cogID=1:numSPs
            tmpsDesignDR = [ones(numVertices_tmpt,1) sp_regressors(:,cogID)];
            tmpB = regress(tmpY, tmpsDesignDR);
            time_series(trID,cogID) = tmpB(2);
        end
    end
end

end