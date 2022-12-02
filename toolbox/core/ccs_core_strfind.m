function [k] = ccs_core_strfind(str,pat)
%CCS_CORE_STRFIND Summary of this function goes here
%   Detailed explanation goes here
k = [];
numStr = numel(str);
for idxStr=1:numStr
    if strcmp(str{idxStr},pat)
        k = [k idxStr];
    end
end