%% dir settings (may not usable for you and you have to change them...)
clear all; clc
hcp_dir = '/Volumes/LaCie/HCP';
ana_dir = [hcp_dir '/zprojects/reliability'];
fig_dir = [ana_dir '/group/figures'];
ccs_dir = '/Users/mac/Projects/CCS';
ccs_matlab = [ccs_dir '/matlab'];
ccs_vistool = [ccs_dir '/vistool'];
fs_home = '/Applications/freesurfer'; 
fsaverage = 'fsaverage5';

%% Set up the path to matlab function in Freesurfer release
addpath(genpath(ccs_matlab)) %ccs matlab scripts
addpath(genpath(ccs_vistool)) %ccs matlab scripts
addpath(genpath([fs_home '/matlab'])) %freesurfer matlab scripts

%% Read in the list of subjects (TRT)
fdemo = [ana_dir '/info/subjects_trt.xls'];
[~, ~, raw] = xlsread(fdemo,'Sheet1');
subjects = raw(2:end,1); nsubs = numel(subjects);
age = cell2mat(raw(2:end,2));
gender = cell2mat(raw(2:end,3));
trt_interval = cell2mat(raw(2:end,4));
subID_3T = cell2mat(raw(2:end,5));

%% Load FreeSurfer table - TRT release
fTable = [hcp_dir '/TRT/info/FreeSurfer.csv'];
table_trt_retest = readtable(fTable);
%write table to an excel file
writetable(table_trt_retest,...
    [ana_dir '/info/fs_trt.xls'],'Sheet',2)
%get feature names
tmp_names = table_trt_retest.Properties.VariableNames;
clear feature_names
for fID=1:numel(tmp_names)
    tmp_name = tmp_names{fID};
    if fID>1
        feature_names{fID,1} = tmp_name(4:end);
    else
        feature_names{fID,1} = tmp_name;
    end
end

%% Load FreeSurfer table - HCP1200 release
flist_1200 = [hcp_dir '/1200/info/behavioral.csv'];
table_1200 = readtable(flist_1200);
subjects_1200 = table_1200.Subject; 
fTable = [hcp_dir '/1200/info/FreeSurfer.csv'];
table_1200_FS = readtable(fTable);
subjects_FS = table_1200_FS.Subject; 
%get subID_3T_MRI
[subjects_trt, subID_3T_FS, subID_trt] = ...
    intersect(subjects_FS,subjects_1200(subID_3T));
%get test table
table_trt_test = table_1200_FS(subID_3T_FS,:);
%write table to an excel file
writetable(table_trt_test,...
    [ana_dir '/info/fs_trt.xls'],'Sheet', 1)
