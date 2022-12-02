%% setup major directories
clear all; clc
disk_dir = '/Volumes/LaCie_1';
hcp_dir = [disk_dir '/HCP'];
parc_dir = [disk_dir '/Schaefer2018_LocalGlobal'];
proj_dir = [hcp_dir '/zprojects/reliability'];
surf_dir = [hcp_dir '/zprojects/surfmodels'];
fig_dir = [proj_dir '/group/figures'];

%% setup major toolboxes
fs_home = '/Applications/freesurfer/7.2';
%connectome computation system
ccs_dir = '/Users/xinian.zuo/Projects/CCS';
ccs_matlab = [ccs_dir '/matlab'];
ccs_vistool = [ccs_dir '/vistool'];

%% add the paths to matlab
addpath(genpath(ccs_matlab)) %ccs matlab scripts
addpath(genpath(ccs_vistool)) %ccs matlab scripts
addpath(genpath([fs_home '/matlab'])) %freesurfer matlab scripts

%% read the subject list
subjects_list = [proj_dir '/info/subjects_3t7t.xlsx'];
[status,sheets,xlFormat] = xlsfinfo(subjects_list);
[num,txt,raw] = xlsread(subjects_list,sheets{1});
subjects_id = num(:,1); nsubs = numel(subjects_id);

%% alff/reho/vmhc computation
data_dir = [hcp_dir '/1200'];
rest_labels = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', ...
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};
idx_subcort = [3:9 12:21]; 
num_subcort = numel(idx_subcort);
bands_label = {'slow6','slow5','slow4','slow3','slow2','slow1'};
nsamples = 1200; rep_time = 0.72;
fbands = ccs_core_lfobands(nsamples,rep_time);
num_bands = numel(fbands);
%start timing
tcpu = cputime;
for idx_sub=1:nsubs
    subject = num2str(subjects_id(idx_sub));
    sub_dir = [data_dir '/' subject '/MNINonLinear/Results'];
    for idx_rest=1:4
        disp(['subject ' subject ': ' rest_labels{idx_rest}])
        func_dir = [sub_dir '/' rest_labels{idx_rest}];
        frest = [func_dir '/boldts_subcort.mat'];
        if ~exist(frest,'file')
            disp(['There is no ' rest_labels{idx_rest} ...
                'for subject ' subject])
        else
            %load data
            load(frest)
            boldalff_subcort = zeros(num_subcort,2,num_bands);
            %loop structure
            for idx_sc=1:num_subcort
                tmpboldts = boldts_subcort{idx_sc};
                %loop bands
                for idx_bd=1:num_bands
                    %alff computation
                    disp(['alff and band-pass filtering: ' ...
                        bands_label{idx_bd} ' ...']) 
                    [tmpalff,tmpalff_nor,~] = ccs_core_ampfilt(...
                        tmpboldts',rep_time,fbands{idx_bd}(1),...
                        fbands{idx_bd}(2));
                    boldalff_subcort(idx_sc,1,idx_bd) = mean(tmpalff);
                    boldalff_subcort(idx_sc,2,idx_bd) = mean(tmpalff_nor);
                end
            end
            %update boldmaps
            fboldmaps = [func_dir '/boldmaps_subcort.mat'];
            tmpboldmaps = load(fboldmaps);
            boldvmhc_subcort = tmpboldmaps.boldvmhc_subcort;
            boldreho_subcort = tmpboldmaps.boldreho_subcort;
            %save maps
            save(fboldmaps,'boldvmhc_subcort',...
                'boldalff_subcort','boldreho_subcort');
        end
    end
end
%end time
ecpu = cputime - tcpu

