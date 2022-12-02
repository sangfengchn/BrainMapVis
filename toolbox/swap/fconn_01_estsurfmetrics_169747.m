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
ccs_dir = '/Users/mac/Projects/CCS';
ccs_matlab = [ccs_dir '/matlab'];
ccs_vistool = [ccs_dir '/vistool'];
%hcp workbench
hcpwkbc_dir = [hcp_dir '/workbench'];
wbcmd = [hcpwkbc_dir '/macosx64_apps/' ...
    'wb_command.app/Contents/MacOS/wb_command'];
ConteAtlas_dir = [surf_dir '/32k_ConteAtlas_v2'];
%cifti toolbox
cifti_matlab = [proj_dir '/matlab/cifti-matlab-master'];

%% add the paths to matlab
addpath(genpath(ccs_matlab)) %ccs matlab scripts
addpath(genpath(ccs_vistool)) %ccs matlab scripts
addpath(genpath(cifti_matlab)) %cifti paths
addpath(genpath([fs_home '/matlab'])) %freesurfer matlab scripts

%% load the geometry of the 32k_ConteAtlas
Conte32k_lh = gifti([ConteAtlas_dir ...
    '/Conte69.L.midthickness.32k_fs_LR.surf.gii']);
nVertices_lh = size(Conte32k_lh.vertices,1);
Conte32k_rh = gifti([ConteAtlas_dir ...
    '/Conte69.R.midthickness.32k_fs_LR.surf.gii']);
nVertices_rh = size(Conte32k_rh.vertices,1);
%brain surface graphs
load([surf_dir '/matlab/Conte69_32k_surfgraph.mat'])

%% read the subject list
subjects_list = [proj_dir '/info/subjects_3t7t.xlsx'];
[status,sheets,xlFormat] = xlsfinfo(subjects_list);
[num,txt,raw] = xlsread(subjects_list,sheets{1});
subjects_id = num(:,1); nsubs = numel(subjects_id);

%% loop subjects from HCP1200 release
data_dir = [hcp_dir '/1200'];
rest_labels = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', ...
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};
idx_subcort = [3:9 12:21]; 
num_subcort = numel(idx_subcort);
bands_label = {'slow6','slow5','slow4','slow3','slow2','slow1'};
num_bands = numel(bands_label);
%start timing
tcpu = cputime;
for idx_sub=45 %169747
    subject = num2str(subjects_id(idx_sub));
    sub_dir = [data_dir '/' subject '/MNINonLinear/Results'];
    for idx_rest=4
        disp(['subject ' subject ': ' rest_labels{idx_rest}])
        func_dir = [sub_dir '/' rest_labels{idx_rest}];
        frest = [func_dir '/' rest_labels{idx_rest} ...
            '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
        if ~exist(frest,'file')
            disp(['There is no ' rest_labels{idx_rest} ...
                'for subject ' subject])
        else
            boldts = ft_read_cifti(frest);
            %extract and save ts in subcortex (skip cerebellumn)
            boldts_subcort = cell(num_subcort,1);
            for idx_sc=1:num_subcort
                boldts_subcort{idx_sc} = boldts.dtseries(...
                    boldts.brainstructure==idx_subcort(idx_sc),:);
            end
            fboldts = [func_dir '/boldts_subcort.mat'];
            save(fboldts,'boldts_subcort');
            %time-freq info
            time = boldts.time; nsamples = length(time); %959 not 1200
            rep_time = time(2)-time(1);
            fbands = ccs_core_lfobands(nsamples,rep_time);
            %extract cerebral cortex bold time series and clean up
            boldts_cc_lh = boldts.dtseries(boldts.brainstructure==1,:);
            boldts_cc_rh = boldts.dtseries(boldts.brainstructure==2,:);
            clear boldts boldts_subcort
            %mask cereb cortex
            [mask_cc_lh, idx_mask_cc_lh] = ccs_core_boldmask(boldts_cc_lh,2);
            [mask_cc_rh, idx_mask_cc_rh] = ccs_core_boldmask(boldts_cc_rh,2);
            mask_cc_hemi = mask_cc_lh.*mask_cc_rh;
            idx_mask_hemi = find(mask_cc_hemi==1);
            fboldmask = [func_dir '/boldmask_cort.mat'];
            save(fboldmask,'mask_cc_lh','mask_cc_rh','mask_cc_hemi');
            %loop bands
            for idx_bd=num_bands
                %alff computation
                disp(['alff and band-pass filtering: ' ...
                    bands_label{idx_bd} ' ...']) 
                [tmpalff, tmpts_filt] = ccs_core_ampfilt(...
                    boldts_cc_lh(idx_mask_cc_lh,:)',...
                    rep_time,fbands{idx_bd}(1),fbands{idx_bd}(2));
                alff_cc_lh = mask_cc_lh;
                alff_cc_lh(idx_mask_cc_lh) = tmpalff;
                ts_lh_filt = boldts_cc_lh';
                ts_lh_filt(:,idx_mask_cc_lh) = tmpts_filt;
                clear tmpalff tmpts_filt
                [tmpalff, tmpts_filt] = ccs_core_ampfilt(...
                    boldts_cc_rh(idx_mask_cc_rh,:)',...
                    rep_time,fbands{idx_bd}(1),fbands{idx_bd}(2));
                alff_cc_rh = mask_cc_rh;
                alff_cc_rh(idx_mask_cc_rh) = tmpalff;
                ts_rh_filt = boldts_cc_rh';
                ts_rh_filt(:,idx_mask_cc_rh) = tmpts_filt;
                clear tmpalff tmpts_filt
                fboldalff = [func_dir '/boldalff_cort.' ...
                    bands_label{idx_bd} '.mat'];
                save(fboldalff,'alff_cc_lh','alff_cc_rh');
                %vmhc computation
                disp(['vmhc estimation: ' bands_label{idx_bd} ' ...'])
                tmpvmhc = ccs_core_mhc(ts_lh_filt(:,idx_mask_hemi),...
                    ts_rh_filt(:,idx_mask_hemi));
                vmhc_cc_lh = mask_cc_lh;
                vmhc_cc_lh(idx_mask_hemi) = tmpvmhc;
                vmhc_cc_rh = vmhc_cc_lh;
                fboldvmhc = [func_dir '/boldvmhc_cort.' ...
                    bands_label{idx_bd} '.mat'];
                save(fboldvmhc,'vmhc_cc_lh','vmhc_cc_rh');
                %reho computation
                disp(['reho estimation: ' bands_label{idx_bd} ' ...'])
                reho_cc_lh = ccs_core_reho(ts_lh_filt,lh_nbrs(:,4));
                reho_cc_rh = ccs_core_reho(ts_rh_filt,rh_nbrs(:,4));
                fboldreho = [func_dir '/boldreho_cort.' ...
                    bands_label{idx_bd} '.mat'];
                save(fboldreho,'reho_cc_lh','reho_cc_rh');
            end
        end
    end
end
%end time
ecpu = cputime - tcpu
