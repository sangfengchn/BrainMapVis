%% dir settings (may not usable for you and you have to change them...)
clear all; clc
hcp_dir = '/Volumes/LaCie/HCP';
proj_dir = [hcp_dir '/zprojects/reliability'];
fig_dir = [proj_dir '/group/figures'];
ccs_dir = '/Users/mac/Projects/CCS';
ccs_matlab = [ccs_dir '/matlab'];
ccs_vistool = [ccs_dir '/vistool'];
fs_home = '/Applications/freesurfer'; 
fsaverage = 'fsaverage5';

%% Set up the path to matlab function in Freesurfer release
addpath(genpath(ccs_matlab)) %ccs matlab scripts
addpath(genpath(ccs_vistool)) %ccs matlab scripts
addpath(genpath([fs_home '/matlab'])) %freesurfer matlab scripts

%% read the subject list
flist_u100 = [hcp_dir, '/1200/info/behavioral_unrelated100.csv'];
table_u100_BH = readtable(flist_u100);
subjects_id = table_u100_BH.Subject;

%% Load FreeSurfer table - HCP1200 release
flist_1200 = [hcp_dir '/1200/info/behavioral.csv'];
table_1200_BH = readtable(flist_1200);
subjects_1200 = table_1200_BH.Subject;
%1200 restricted info
flist_1200R = [hcp_dir '/1200/info/restricted.csv'];
table_1200R_BH = readtable(flist_1200R);
subjects_1200R = table_1200R_BH.Subject; 
%get demographic info
[~, subID_3T_BH, ~] = intersect(subjects_1200,subjects_id);
sex = table_1200_BH.Gender(subID_3T_BH);
[~, subID_3T_RBH, ~] = intersect(subjects_1200R,subjects_id);
age = table_1200R_BH.Age_in_Yrs(subID_3T_RBH);
hand = table_1200R_BH.Handedness(subID_3T_RBH);
%get freesurfer subcortical volumes
fTable = [hcp_dir '/1200/info/FreeSurfer.csv'];
table_1200_FS = readtable(fTable);
subjects_FS = table_1200_FS.Subject; 
[~, subID_3T_FS, ~] = intersect(subjects_FS,subjects_id);
table_u100 = table_1200_FS(subID_3T_FS,:);

%% Estimate Laterality for Subcortical Regions
ICV_u100 = table_u100.FS_IntraCranial_Vol;
idx_subcort = [25:28 32 33 35 36 43:50];
table_u100_subcort = table_u100(:,idx_subcort);
vol_u100_subcort_L = table2array(table_u100_subcort(:,1:8));
vol_u100_subcort_R = table2array(table_u100_subcort(:,9:16));
LvR_subcort = vol_u100_subcort_L - vol_u100_subcort_R;
LaR_subcort = vol_u100_subcort_L + vol_u100_subcort_R;
AI_u100_subcort = 100*(2*LvR_subcort./LaR_subcort);

%% Paired stats tests on laterality
h = zeros(8,1); p = ones(8,1);
ci = zeros(8,2); stats = cell(8,1);
for idx=1:8
    tmpL = vol_u100_subcort_L(:,idx);
    tmpR = vol_u100_subcort_R(:,idx);
    [h(idx),p(idx),ci(idx,:),stats{idx,1}] = ttest(tmpL,tmpR,'Alpha',0.05/8);
end

%% PreStat tests effects of age,sex,hand and icv on AI
[r_age,p_age] = corr(AI_u100_subcort,age);
[r_hand,p_hand] = corr(AI_u100_subcort,hand);
[r_icv,p_icv] = corr(AI_u100_subcort,ICV_u100);

%% Sex-related differences in AI
idx_female = contains(sex,'F');
idx_male = contains(sex,'M');
cov_female_dm = IPN_demean([age(idx_female) ...
    hand(idx_female) ICV_u100(idx_female)]);
cov_male_dm = IPN_demean([age(idx_male) ...
    hand(idx_male) ICV_u100(idx_male)]);
%regress covariates seperately for males and females
hAI = zeros(8,1); pAI = ones(8,1);
ciAI = zeros(8,2); statsAI = cell(8,1);
for idx=1:8
    tmpAI = AI_u100_subcort(:,idx);
    %female
    tmpAI_female = tmpAI(idx_female);
    [~,~,tmpAI_female_adj] = regress(tmpAI_female,cov_female_dm);
    %male
    tmpAI_male = tmpAI(idx_male);
    [~,~,tmpAI_male_adj] = regress(tmpAI_male,cov_male_dm);
    %two-sample t-tests
    [hAI(idx),pAI(idx),ciAI(idx,:),statsAI{idx,1}] = ...
        ttest2(tmpAI_female_adj,tmpAI_male_adj,'Alpha',0.05/8);
end
%regress covariates together
cov_dm = IPN_demean([age hand ICV_u100]);
hAI_all = zeros(8,1); pAI_all = ones(8,1);
ciAI_all = zeros(8,2); statsAI_all = cell(8,1);
for idx=1:8
    tmpAI = AI_u100_subcort(:,idx);
    [~,~,tmpAI_adj] = regress(tmpAI,cov_dm);
    %female
    tmpAI_female_adj = tmpAI_adj(idx_female);
    %male
    tmpAI_male_adj = tmpAI_adj(idx_male);
    %two-sample t-tests
    [hAI_all(idx),pAI_all(idx),ciAI_all(idx,:),statsAI_all{idx,1}] = ...
        ttest2(tmpAI_female_adj,tmpAI_male_adj,'Alpha',0.05/8);
end
