%% setup major directories
clear; clc
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

%% load the RSNs and their Confidence Intervals
Conte32k_RSN_lh = gifti([ConteAtlas_dir ...
    '/RSN-networks.L.32k_fs_LR.label.gii']); %lh
RSN_cdata_lh = Conte32k_RSN_lh.cdata(:,1);
numVertices_lh = numel(RSN_cdata_lh);
RSN_names_lh = Conte32k_RSN_lh.labels.name;
RSN_keys_lh = Conte32k_RSN_lh.labels.key;
Conte32k_RSN_CI_lh = gifti([ConteAtlas_dir ...
    '/Yeo2011_7NetworksConfidence_N1000.L.32k_fs_LR.gii']);
DR_tmpt_lh = []; RSN_idx_lh = [];
for idxRSN=1:7
    nameRSN = ['7Networks_' num2str(idxRSN)];
    idxName = ccs_core_strfind(RSN_names_lh,nameRSN);
    idxKey = RSN_keys_lh(idxName);
    idxTMPT = find(RSN_cdata_lh==idxKey);
    DR_tmpt_lh{idxRSN} = Conte32k_RSN_CI_lh.cdata(idxTMPT);
    RSN_idx_lh{idxRSN} = idxTMPT;
end
%------------
Conte32k_RSN_rh = gifti([ConteAtlas_dir ...
    '/RSN-networks.R.32k_fs_LR.label.gii']); %rh
RSN_cdata_rh = Conte32k_RSN_rh.cdata(:,1);
numVertices_rh = numel(RSN_cdata_rh);
RSN_names_rh = Conte32k_RSN_rh.labels.name;
RSN_keys_rh = Conte32k_RSN_rh.labels.key;
Conte32k_RSN_CI_rh = gifti([ConteAtlas_dir ...
    '/Yeo2011_7NetworksConfidence_N1000.R.32k_fs_LR.gii']);
DR_tmpt_rh = []; RSN_idx_rh = [];
for idxRSN=1:7
    nameRSN = ['7Networks_' num2str(idxRSN)];
    idxName = ccs_core_strfind(RSN_names_rh,nameRSN);
    idxKey = RSN_keys_rh(idxName);
    idxTMPT = find(RSN_cdata_rh==idxKey);
    DR_tmpt_rh{idxRSN} = Conte32k_RSN_CI_rh.cdata(idxTMPT);
    RSN_idx_rh{idxRSN} = idxTMPT;
end

%% read the subject list
subjects_list = [proj_dir '/info/subjects_3t7t.xlsx'];
[status,sheets,xlFormat] = xlsfinfo(subjects_list);
[num,txt,raw] = xlsread(subjects_list,sheets{1});
subjects_id = num(:,1); nsubs = numel(subjects_id);

%% loop subjects from HCP1200 release
data_dir = [hcp_dir '/1200'];
rest_labels = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', ...
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};
bands_label = {'slow6','slow5','slow4','slow3','slow2','slow1'};
num_bands = numel(bands_label);
%start timing
tcpu = cputime;
for idx_sub=1:nsubs
    subject = num2str(subjects_id(idx_sub));
    sub_dir = [data_dir '/' subject '/MNINonLinear/Results'];
    for idx_rest=1:4
        disp(['subject ' subject ': ' rest_labels{idx_rest}])
        func_dir = [sub_dir '/' rest_labels{idx_rest}];
        frest = [func_dir '/' rest_labels{idx_rest} ...
            '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
        if ~exist(frest,'file')
            disp(['There is no ' rest_labels{idx_rest} ...
                'for subject ' subject])
        else
            boldts = ft_read_cifti(frest);     
            %time-freq info
            time = boldts.time; nsamples = length(time);
            rep_time = time(2)-time(1);
            fbands = ccs_core_lfobands(nsamples,rep_time);
            %extract cerebral cortex bold time series and clean up
            boldts_cc_lh = boldts.dtseries(boldts.brainstructure==1,:);
            boldts_cc_rh = boldts.dtseries(boldts.brainstructure==2,:);
            clear boldts
            %mask cereb cortex
            fboldmask = [func_dir '/boldmask_cort.mat'];
            load(fboldmask);
            idx_mask_cc_lh = find(mask_cc_lh==1);
            idx_mask_cc_rh = find(mask_cc_rh==1);
            %loop bands
            for idx_bd=3 %slow4
                %band-pass filtering
                disp(['band-pass filtering: ' ...
                    bands_label{idx_bd} ' ...']) 
                [~,~,tmpts_filt_lh] = ccs_core_ampfilt(...
                    boldts_cc_lh(idx_mask_cc_lh,:)',...
                    rep_time,fbands{idx_bd}(1),fbands{idx_bd}(2),'true');
                ts_filt_lh = zeros(nsamples,numVertices_lh);
                ts_filt_lh(:,idx_mask_cc_lh) = tmpts_filt_lh;
                [~,~,tmpts_filt_rh] = ccs_core_ampfilt(...
                    boldts_cc_rh(idx_mask_cc_rh,:)',...
                    rep_time,fbands{idx_bd}(1),fbands{idx_bd}(2),'true');
                ts_filt_rh = zeros(nsamples,numVertices_rh);
                ts_filt_rh(:,idx_mask_cc_rh) = tmpts_filt_rh;
                ts_filt_masked = [tmpts_filt_lh tmpts_filt_rh];
                clear tmpts_filt_lh tmpts_filt_rh
                %dual regression
                disp(['dual regression estimation: ' bands_label{idx_bd} ' ...'])
                dr_cc_lh = zeros(numVertices_lh,7);
                dr_cc_rh = zeros(numVertices_rh,7);
                dr_timeseries = zeros(nsamples,7);
                for idxTMPT=1:7
                    DR_tmpt = [DR_tmpt_lh{idxTMPT}; DR_tmpt_rh{idxTMPT}];
                    ts_filt_tmpt = [ts_filt_lh(:,RSN_idx_lh{idxTMPT}) ...
                        ts_filt_rh(:,RSN_idx_rh{idxTMPT})];
                    dr_timeseries(:,idxTMPT) = ccs_core_dualreg_space(ts_filt_tmpt',DR_tmpt);
                end
                tmp_maps = ccs_core_dualreg_time(ts_filt_masked',dr_timeseries);
                dr_cc_lh(idx_mask_cc_lh,:) = tmp_maps(1:numel(idx_mask_cc_lh),:);
                dr_cc_rh(idx_mask_cc_rh,:) = tmp_maps((1+numel(idx_mask_cc_lh)):end,:);
                %save the data
                fbolddualreg = [func_dir '/bolddualreg_cort.' ...
                    bands_label{idx_bd} '.mat'];
                save(fbolddualreg,'dr_cc_lh','dr_cc_rh','dr_timeseries');
            end
        end
    end
end
%end time
ecpu = cputime - tcpu

%% test the surface render
%cmax = max(abs(tmp_maps(:))); cmap_range = [-cmax cmax]*0.75;

%figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
%SurfStatView(dr_cc_lh(:,7), surfConte69_lh, ' ', 'white', 'true');
%SurfStatColLim(cmap_range); set(gcf, 'PaperPositionMode', 'auto');            

%figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
%SurfStatView(dr_cc_rh(:,7), surfConte69_rh, ' ', 'white', 'true'); 
%SurfStatColLim(cmap_range); set(gcf, 'PaperPositionMode', 'auto');
