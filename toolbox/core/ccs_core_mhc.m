function [ mhc ] = ccs_core_mhc(ts_lh, ts_rh)
%CCS_CORE_MHC MATLAB version of Mirrored Hemispheric Connectivity
%   Detailed explanation goes here
%   INPUTS
%       ts_lh - a TxN matrix of time series for left hemi
%       ts_rh - a TxN matrix of time series for right hemi
[ntp_lh, nsp_lh] = size(ts_lh);
[ntp_rh, nsp_rh] = size(ts_rh);
if ntp_lh ~= ntp_rh
    disp('The time points must be same between the two hemi.')
else
    if nsp_lh ~= nsp_rh
        disp('The space points must be same between the two hemi.')
    else
        z_ts_lh = zscore(ts_lh,1);
        z_ts_rh = zscore(ts_rh,1);
        mhc = sum(z_ts_lh.*z_ts_rh)/ntp_lh;
    end
end

end

