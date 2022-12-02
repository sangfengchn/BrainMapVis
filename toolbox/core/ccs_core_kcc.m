function [ kcc ] = ccs_core_kcc(tsmat)
%CCS_CORE_KCC Kendal's W Computation
%
%   Detailed explanation:
%    INPUT:
%       tsmat -- original time series matrix
%   Credits:
%      Xi-Nian Zuo, PhD of Applied Mathematics
%      Beijing Normal University
%      Email: xinian.zuo@bnu.edu.cn
%      https://zuoxinian.github.io

[ntp, nsp] = size(tsmat); %num of time points, num of spatial points
[~,I]=sort(tsmat); [~,R]=sort(I);
S=sum(sum(R,2).^2)-ntp*mean(sum(R,2)).^2;
F=nsp*nsp*(ntp*ntp*ntp-ntp); kcc=12*S/F;

