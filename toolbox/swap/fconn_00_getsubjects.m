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

%% get subjects list - 1200 release
flist_1200 = [hcp_dir '/1200/info/behavioral.csv'];
table_1200 = readtable(flist_1200);
subjects = table_1200.Subject; 
nsubs = numel(subjects);
gender = table_1200.Gender; 
flist_restricted = [hcp_dir '/1200/info/restricted.csv'];
table_restricted = readtable(flist_restricted);
subjects_restricted = table_restricted.Subject;
%if subjects are identical: YES for 0
subsdiff = sum(subjects - subjects_restricted);
if subsdiff==0
    age = table_restricted.Age_in_Yrs;
end
%numberize gender as sex
sex = zeros(nsubs,1);
% numberize
for sid=1:nsubs
    if strcmp(gender{sid}, 'M')
        sex(sid,1) = 1;
    else
        sex(sid,1) = 0;
    end
end

%% get subjects list -  7T 184 release
flist_7t184 = [hcp_dir '/7T/info/behavioral.csv'];
table_7t184 = readtable(flist_7t184);
subjects_7t184 = table_7t184.Subject;
gender_7t184 = table_7t184.Gender;
x3T = table_7t184.x3T_Full_MR_Compl;
flist_rsct7t = [hcp_dir '/7T/info/restricted.csv'];
table_rsct7t = readtable(flist_rsct7t);
subjects_rsct7t = table_rsct7t.Subject;
%if subjects are identical: YES for 0
subsdiff = sum(subjects_7t184 - subjects_rsct7t);
if subsdiff==0
    age_7t184 = table_rsct7t.Age_in_Yrs;
end
%numberize gender as sex
sex_7t184 = zeros(184,1);
% numberize
for sid=1:184
    if strcmp(gender_7t184{sid}, 'M')
        sex_7t184(sid,1) = 1;
    else
        sex_7t184(sid,1) = 0;
    end
end

%% get subjects list - 100 unrelated
flist_u100 = [hcp_dir '/1200/info/behavioral_unrelated100.csv'];
table_u100 = readtable(flist_u100);
subjects_u100 = table_u100.Subject; 

%% get subjects list - unrelated 3T&7T
familyID_rsct7t = table_rsct7t.Family_ID;
familyID_rsct7t_unique = ccs_core_uniquestrcell(familyID_rsct7t);
numSubs_uniq7t = numel(familyID_rsct7t_unique);
subjetcs_uniq7t = zeros(numSubs_uniq7t,1);
subID_uniq = 0; subID_uniq7t = zeros(numSubs_uniq7t,1);
for famID=1:numSubs_uniq7t
    tmpfamilyID = familyID_rsct7t_unique{famID};
    for subID=1:numel(familyID_rsct7t)
        if strcmp(tmpfamilyID,familyID_rsct7t{subID})
            subID_uniq = subID_uniq + 1;
            subID_uniq7t(subID_uniq) = subID;
            subjetcs_uniq7t(subID_uniq) = subjects_rsct7t(subID);
            break
        end
    end
end
%intersection of 3T and 7T
[subjects_3t7t, subID_1200, subID_uniq7t] = ...
    intersect(subjects,subjetcs_uniq7t);
sex_3t7t = sex(subID_1200);
age_3t7t = age(subID_1200);
%define a table
table_3t7t = table(subjects_3t7t,age_3t7t,sex_3t7t,...
    subID_1200,subID_uniq7t);
table_3t7t.Properties.VariableNames = {'Subject','Age','Gender',...
    'subID_3T','subID_7T'};
%save a table to excell file
%fout = [ana_dir '/info/subjects_3t7t.xls'];
%writetable(table_3t7t,fout,'FileType','spreadsheet')

%% get subjects list - test-retest
flist_trt = [hcp_dir '/TRT/info/behavioral.csv'];
table_trt = readtable(flist_trt);
subjects_trt = table_trt.Subject;
numSubs_trt = numel(subjects_trt);
flist_rsct_trt = [hcp_dir '/TRT/info/restricted.csv'];
table_rsct_trt = readtable(flist_rsct_trt);
trtinterval = table_rsct_trt.TestRetestInterval; %days
%intersection of 3T and TRT
[subjects_trt_final, subID_1200, subID_trt] = ...
    intersect(subjects,subjects_trt);
age_trt = age(subID_1200);
sex_trt = sex(subID_1200);
%define a table
table_trt_final = table(subjects_trt_final,age_trt,sex_trt,trtinterval,...
    subID_1200,subID_trt);
table_trt_final.Properties.VariableNames = {'Subject','Age','Gender',...
    'TRT_interval','subID_3T','subID_TRT'};
%save a table to excel file
%fout = [ana_dir '/info/subjects_trt.xls'];
%writetable(table_trt_final,fout,'FileType','spreadsheet')

%% get subjects list - trt for both 3T and 7T
[subjects_trt_3t7t, subID_3t7t, subID_trt] = ...
    intersect(subjects_3t7t,subjects_trt_final);
age_trt_3t7t = age_trt(subID_trt);
sex_trt_3t7t = sex_trt(subID_trt);
trtinterval_3t7t = trtinterval(subID_trt);
%define a table
table_trt_3t7t_final = table(subjects_trt_3t7t,age_trt_3t7t,...
    sex_trt_3t7t,trtinterval_3t7t,subID_3t7t,subID_trt);
table_trt_final_3t7t.Properties.VariableNames = ...
    {'Subject','Age','Gender','TRT_interval','subID_3T7T','subID_TRT'};
%save a table to excell file
fout = [ana_dir '/info/subjects_trt_3t7t.xls'];
writetable(table_trt_3t7t_final,fout,'FileType','spreadsheet')

%% get subjects list - trt for both 3T and 7T: REST data
proj_dir = ana_dir;
subjects_list = [proj_dir '/info/subjects_3t7t.xlsx'];
[status,sheets,xlFormat] = xlsfinfo(subjects_list);
[num,txt,raw] = xlsread(subjects_list,sheets{1});
subjects_id = num(:,1); nsubs = numel(subjects_id);
%data directories and rest sessions
data_dir_1200 = [hcp_dir '/1200'];
rest1200_labels = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', ...
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};
data_dir_7T = [hcp_dir '/7T'];
rest7T_labels = {'rfMRI_REST1_7T_PA', 'rfMRI_REST2_7T_AP', ...
    'rfMRI_REST3_7T_PA', 'rfMRI_REST4_7T_AP'};
%loop subjects
subs_rest_3t = zeros(nsubs,4);
subs_rest_7t = zeros(nsubs,4);
for idx_sub=1:nsubs
    subject = num2str(subjects_id(idx_sub));
    sub_dir_3t = [data_dir_1200 '/' subject '/MNINonLinear/Results'];
    sub_dir_7t = [data_dir_7T '/' subject '/MNINonLinear/Results'];
    for idx_rest=1:4
        disp(['subject ' subject ': REST' num2str(idx_rest) ' ...'])
        %3T
        func_dir = [sub_dir_3t '/' rest1200_labels{idx_rest}];
        frest = [func_dir '/' rest1200_labels{idx_rest} ...
            '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
        if exist(frest,'file')
            subs_rest_3t(idx_sub,idx_rest) = 1;
        end
        %7T
        func_dir = [sub_dir_7t '/' rest7T_labels{idx_rest}];
        frest = [func_dir '/' rest7T_labels{idx_rest} ...
            '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
        if exist(frest,'file')
            subs_rest_7t(idx_sub,idx_rest) = 1;
        end
    end
end
subs_rest_3t7t = subs_rest_3t.*subs_rest_7t;
tmpsubjects = sum(subs_rest_3t7t,2);
idx_subs_rest_3t7t = find(tmpsubjects==4);
%save subject indices
save([proj_dir '/info/subidx_3t7t_rest4sess_raw.mat'],'idx_subs_rest_3t7t')
