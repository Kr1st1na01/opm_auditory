%% Reset all
clear
close all
clc
restoredefaultpath

%% Base paths
if contains(pwd,'/home/krimat')
    % Server:
    base_data_path = '/archive/24110_opm_auditory/';
    base_save_path = '/home/krimat/Documents/aud_data';
    base_matlab_path = '/home/krimat/Documents/MATLAB/';
    project_scripts_path = '/home/krimat/Documents/MATLAB/opm_auditory';
else
    % Laptop:
    base_data_path = 'C:/Users/Kristina/Documents/KTH dokument/MEX/24110_opm_auditory';
    base_save_path = 'C:/Users/Kristina/Documents/KTH dokument/MEX/Resultat';
    base_matlab_path = 'C:/Users/Kristina/Documents/MATLAB';
    project_scripts_path = 'C:/Users/Kristina/Documents/MATLAB/opm_auditory';
end

%% Set up fieldtrip
addpath ('/home/krimat/Documents/MATLAB/fieldtrip') % Fieldtrip path
%addpath(fullfile(base_matlab_path,'fieldtrip_private')) % Fieldtrip private functions
addpath(project_scripts_path)
ft_defaults

global ft_default
ft_default.showcallinfo = 'no';

%% Params
overwrite = [];
overwrite.preproc = true;
overwrite.coreg = false;
overwrite.mri = false;
overwrite.dip = true;
overwrite.mne = true;

params = [];
params.pre = 0.3; %sec
params.post = 0.9; %sec
params.pad = 0.2; %sec
params.filter = [];
params.filter.hp_freq = 1;
params.filter.lp_freq = 50;
params.filter.bp_freq = [];
params.filter.notch = sort([50 60]);
params.n_comp = 40;
params.ica_threshold = 0.8; % cutoff for EOG/ECG coherence
params.z_threshold = 20;
params.corr_threshold = 0.7; % correlation threshold for badchannel neighbors
params.opm_std_threshold = 2.5e-12;
params.eeg_std_threshold = 1e-4;
params.squidmag_std_threshold = 5e-12;
params.squidgrad_std_threshold = 5e-11;
params.hpi_freq = 33;

params.trigger_code = [1 3 5 11 13];
params.trigger_labels = {'Std' 'Low NG' 'Low Go' 'High NG' 'High Go'}; % Normal, No-Go, Go, High No-Go, High Go

params.oldtrigger_code = [1 18 20 10 12]; % The old trigger codes
params.oldtrigger_labels = ['Std' 'NoGo' 'Go' 'HighNoG0o' 'HighGo'];

%% Subjects + dates
subses = {'0005' '240208';
    '0905' '240229';
    '0916' '240320';
    '0953' '241104';
    '1096' '241022';
    '1153' '240321';
    '1167' '240425';
    '1186' '240925';
    '1190' '241023'; 
    '1191' '241024';
    '1193' '241029';
    '1194' '241029';
    '1195' '241030'};
mri_files = {'00000001.dcm' 
    '/mri/sub-15931_T1w.nii.gz'  
    '/nifti/anat/sub-15985_T1w.nii.gz'};

%% Loop over subjects
for i_sub = 1:size(subses,1)
    params.sub = ['sub_' num2str(i_sub,'%02d')];

    %% Paths
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    if ~exist(save_path, 'dir')
       mkdir(base_save_path)
    end
    if ~exist(fullfile(save_path,'figs'), 'dir')
       mkdir(fullfile(save_path,'figs'))
    end
    meg_file = fullfile(raw_path, 'meg', 'AudOddMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
    opm_file = fullfile(raw_path, 'osmeg', 'AudOddOPM_raw.fif');
    aux_file = fullfile(raw_path, 'meg', 'AudOddEEG.fif');
    hpi_file = fullfile(raw_path, 'osmeg', 'HPIpre_raw.fif');

    %% Read ICA files

    % OPM-MEG
    if exist(fullfile(save_path, [params.sub '_opmeeg_ica.mat']),'file') && overwrite.preproc==false % om filerna existerar så sker inte en overwrite och filerna laddas in
       load(fullfile(save_path, [params.sub '_opm_ica.mat']))
       opm_ica = data_ica_ds;
       load(fullfile(save_path, [params.sub '_opmeeg_ica.mat']))
       opmeeg_ica = data_ica_ds;
       %SQUID-MEG
       load(fullfile(save_path, [params.sub '_squid_ica.mat']))
       squidmag_ica = data_ica_ds;
       load(fullfile(save_path, [params.sub '_squideeg_ica.mat']))        
       squideeg_ica = data_ica_ds;
       clear data_ica_ds
    else
    
        ft_hastoolbox('mne', 1);
 
        % Read data
        [opm_cleaned, opmeeg_cleaned] = read_osMEG(opm_file, aux_file, save_path, params); % Read data
        [squid_cleaned, squideeg_cleaned] = read_cvMEG(meg_file, save_path, params); % Read data

        if i_sub <=3 % Change trigger codes in old recordings
            for i_trl = 2:length(params.trigger_code)
                opmeeg_cleaned.trialinfo(opmeeg_cleaned.trialinfo==params.oldtrigger_code(i_trl)) = params.trigger_code(i_trl); % A(A==yourvalue)=NewValue;
                opm_cleaned.trialinfo(opm_cleaned.trialinfo==params.oldtrigger_code(i_trl)) = params.trigger_code(i_trl);
                squideeg_cleaned.trialinfo(squideeg_cleaned.trialinfo==params.oldtrigger_code(i_trl)) = params.trigger_code(i_trl); 
                squid_cleaned.trialinfo(squid_cleaned.trialinfo==params.oldtrigger_code(i_trl)) = params.trigger_code(i_trl);
            end
        end
 clear i_trl

     %% ICA
        % ICA OPM
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        opm_ica = ica_MEG(opm_cleaned, save_path, params);

        cfg = [];
        cfg.elec = opmeeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_opmeeg_layout.mat']);
        opmeeg_layout = ft_prepare_layout(cfg);
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'opmeeg';
        opmeeg_ica = ica_MEG(opmeeg_cleaned, save_path, params);
    
        % ICA SQUID
        params.modality = 'squid';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'MEG*';
        squid_ica = ica_MEG(squid_cleaned, save_path, params); % Both mag and grad

        cfg = [];
        cfg.elec = squideeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_megeeg_layout.mat']);
        megeeg_layout = ft_prepare_layout(cfg);
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'squideeg';
        squideeg_ica = ica_MEG(squideeg_cleaned, save_path, params);

%% Dividing data for MMN och TFR
        % MMN data is generated first (= cropped data) and TFR data is generated after
        cfg = [];
        cfg.lpfilter  = 'yes';        % Apply lowpass filter
        cfg.lpfreq    = 20;      
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.100 0];
        cfg.toilim = [-0.100 0.5];
        
        params.ica = [opm_ica, opmeeg_ica, squid_ica, squideeg_ica];
        params.ica_labels = {'opm', 'opmeeg', 'squid', 'squideeg'};
        
        for i = 1:length(params.ica)
            cfg = [];
            cfg.lpfilter  = 'yes';        % Apply lowpass filter
            cfg.lpfreq    = 20;      
            cfg.demean          = 'yes';
            cfg.baselinewindow  = [-0.100 0];
            cfg.toilim = [-0.100 0.5];
            
            filt_downsampled_data = ft_preprocessing(cfg, params.ica(i)); % Data is processed with a lowpass filter of 40 Hz
            cropped_data = ft_redefinetrial(cfg, filt_downsampled_data); % Trials are cropped into time of interest (toi)
            save(fullfile(save_path, [params.sub '_' params.ica_labels{i} '_cropped_ica']), 'cropped_data',"-v7.3");

            cfg = [];
            cfg.hpfilter  = 'yes';        
            cfg.hpfreq    = 30;      % Apply highpass filter
            cfg.demean          = 'yes';
            cfg.baselinewindow  = [-0.100 0];

            TFR_data = ft_preprocessing(cfg, params.ica(i)); % Data is processed with a lowpass filter of 30 Hz
            save(fullfile(save_path, [params.sub '_' params.ica_labels{i} '_TFR_ica']), 'TFR_data',"-v7.3");
        end
        close all
        clear i
        clear -regexp ^TFR ^MMN
    end
end

%% Loop over subjects for timelocking

for i_sub = 1:size(subses,1)
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    
    %% Paths
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    if ~exist(save_path, 'dir')
       mkdir(base_save_path)
    end
    if ~exist(fullfile(save_path,'figs'), 'dir')
       mkdir(fullfile(save_path,'figs'))
    end
    meg_file = fullfile(raw_path, 'meg', 'AudOddMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
    opm_file = fullfile(raw_path, 'osmeg', 'AudOddOPM_raw.fif');
    aux_file = fullfile(raw_path, 'meg', 'AudOddEEG.fif');
    hpi_file = fullfile(raw_path, 'osmeg', 'HPIpre_raw.fif');

        %% Overwrite for timelock 
    if exist(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']),'file') && overwrite.preproc==false % om filerna existerar så sker inte en overwrite och filerna laddas in
        load(fullfile(save_path, [params.sub '_opm_timelocked.mat']))
        opm_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']))
        opmeeg_timelocked = timelocked;
        % SQUID-MEG
        load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat']))
        squidmag_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_squidgrad_timelocked.mat']))
        squidgrad_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_squideeg_timelocked.mat']))
        squideeg_timelocked = timelocked;
        clear timelocked
    else
        
        ft_hastoolbox('mne', 1);
       
        % Load data
        load(fullfile(save_path, [params.sub '_opm_cropped_ica.mat']))
        MMN_opm = cropped_data;
        load(fullfile(save_path, [params.sub '_opm_TFR_ica.mat']))
        TFR_opm = TFR_data;

        load(fullfile(save_path, [params.sub '_opmeeg_cropped_ica.mat']))
        MMN_opmeeg = cropped_data;
        load(fullfile(save_path, [params.sub '_opmeeg_TFR_ica.mat']))
        TFR_opmeeg = TFR_data;
    
        load(fullfile(save_path, [params.sub '_squid_cropped_ica.mat']))
        MMN_squid = cropped_data;
        load(fullfile(save_path, [params.sub '_squid_TFR_ica.mat']))
        TFR_squid = TFR_data;

        load(fullfile(save_path, [params.sub '_squideeg_cropped_ica.mat']))
        MMN_squideeg = cropped_data;
        load(fullfile(save_path, [params.sub '_squideeg_TFR_ica.mat']))
        TFR_squideeg = TFR_data;

        clear cropped_data
        clear TFR_data
        clear labels


        %% Average for OPM-MEG
        peak = [];
        peak.values = [];
        peak.labels = {};
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        params.pretimwin = 0.027;
        params.posttimwin = 0.033;
        params.freqsteps = 0.2;
        params.freqwin = [0 0.8];
        [opm_timelocked, peak] = timelock_MEG(MMN_opm, TFR_opm, params, save_path, peak); % Timelockar vanlig MMN och plottar för Std, Low och High och kör freqanalysis på TFR
        close all
        
        load(fullfile(save_path, [params.sub '_opmeeg_layout.mat']))
        opmeeg_layout = layout;
        params.modality = 'opmeeg';
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        [opmeeg_timelocked, peak] = timelock_MEG(MMN_opmeeg, TFR_opmeeg, params, save_path, peak); 
        close all

        %% Average SQUID-MEG

        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        [squidmag_timelocked, peak] = timelock_MEG(MMN_squid, TFR_squid, params, save_path, peak); 
        close all

        params.modality = 'squidgrad';
        params.layout = 'neuromag306planar.lay';
        params.chs = 'megplanar';
        params.amp_scaler = 1e15/100;
        params.amp_label = 'B [fT/cm]';
        [squidgrad_timelocked, peak] = timelock_MEG(MMN_squid, TFR_squid, params, save_path, peak); 
        close all

        load(fullfile(save_path, [params.sub '_megeeg_layout.mat']))
        megeeg_layout = layout;
        params.modality = 'squideeg';
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        [squideeg_timelocked, peak] = timelock_MEG(MMN_squideeg, TFR_squideeg, params, save_path, peak); 
        save(fullfile(save_path, [params.sub '_peaks']), 'peak' ,"-v7.3");
        close all
        clear -regexp ^TFR ^MMN ^opm ^opmeeg ^squid ^squidgrad ^squidmag ^squideeg layout

    end
% %%
%     params = rmfield(params,{'modality', 'layout', 'chs', 'amp_scaler', 'amp_label'}); % remove fields used for picking modality
%     create_bads_reports(base_save_path, i_sub, params);
%     close all
end

%% --- Group sensor level -------------------------------------------------
if ~exist(fullfile(base_save_path,'figs'), 'dir')
       mkdir(fullfile(base_save_path,'figs'))
end
subs = 1:size(subses,1);
sensor_results_goup(base_save_path,subs, params)

%% Prepare MRIs
for i_sub = 2:size(subses,1) % Skippar första sub ?
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}],'mri');
    if exist(fullfile(save_path, 'headmodels.mat'),'file') && overwrite.mri==false
        load(fullfile(save_path, 'headmodels.mat'));
        load(fullfile(save_path, 'meshes.mat'));
        load(fullfile(save_path, 'mri_resliced.mat'));
    else
        meg_file = fullfile(raw_path, 'meg', 'MEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'MEG_proc-tsss+corr98.fif');
        end
        mri_file = fullfile(mri_path, 'orig','001.mgz');
        [headmodels, meshes] = prepare_mri(mri_file,meg_file,save_path);
        close all
    end
end

%% HPI localization
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    hpi_file = fullfile(raw_path, 'osmeg', 'HPIpre_raw.fif');

    if exist(fullfile(save_path, 'opm_trans.mat'),'file') && overwrite.coreg==false
        load(fullfile(save_path, 'hpi_fit.mat'));
        load(fullfile(save_path, 'opm_trans.mat'));
    else
        meg_file = fullfile(raw_path, 'meg', 'MEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'MEG_proc-tsss+corr98.fif');
        end
        ft_hastoolbox('mne', 1);
        load(fullfile(save_path, [params.sub '_opm_ica_ds']));
        params.include_chs = data_ica_ds.label(find(contains(data_ica_ds.label,'bz')));
        params.include_chs = params.include_chs([1:75 77:end]);
        [hpi_fit, opm_trans, hpi_fit_tf] = fit_hpi(hpi_file, meg_file, save_path, params);
        close all
    end
end

% Transform for OPM data
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    if exist(fullfile(save_path, 'headmodels.mat'),'file') && exist(fullfile(save_path, 'opm_trans.mat'),'file') && overwrite.coreg==true
        load(fullfile(save_path, 'opm_trans.mat'));
        load(fullfile(save_path, 'headmodels.mat'));
        load(fullfile(save_path, 'meshes.mat'));
        load(fullfile(save_path, [params.sub '_opm_timelocked.mat']))
        opm_timelocked = timelocked;
        opm_timelockedT = opm_timelocked;
        load(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']))
        opmeeg_timelocked = timelocked;
        opmeeg_timelockedT = opmeeg_timelocked;
        load(fullfile(save_path, [params.sub '_squideeg_timelocked.mat']))
        squideeg_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat']))
        squidmag_timelocked = timelocked;
        clear timelocked;
        meg_file = fullfile(raw_path, 'meg', 'MEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'MEG_proc-tsss+corr98.fif');
        end
        headshape = ft_read_headshape(meg_file);
        for i = 1:5
            opm_timelockedT{i}.grad.chanpos = opm_trans.transformPointsForward(opm_timelocked{i}.grad.chanpos*1e2)*1e-2;
            opm_timelockedT{i}.grad.coilpos = opm_trans.transformPointsForward(opm_timelocked{i}.grad.coilpos*1e2)*1e-2;
            opmeeg_timelockedT{i}.elec.chanpos = squideeg_timelocked{i}.elec.chanpos;
            opmeeg_timelockedT{i}.elec.elecpos = squideeg_timelocked{i}.elec.elecpos;
        end
        
        h = figure; 
        hold on; 
        ft_plot_sens(opm_timelockedT{1}.grad,'unit','cm')
        ft_plot_sens(opmeeg_timelockedT{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
        ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.7,'FaceColor',[229 194 152]/256,'unit','cm')
        ft_plot_headmodel(headmodels.headmodel_meg)
        ft_plot_headshape(headshape)
        hold off;
        title('OPM-MEG')
        view([-140 10])
        savefig(h, fullfile(save_path, 'figs', 'opm_layout.fig'))
        saveas(h, fullfile(save_path, 'figs', 'opm_layout.jpg'))
    
        h = figure; 
        hold on
        ft_plot_sens(squidmag_timelocked{1}.grad,'unit','cm')
        ft_plot_sens(squideeg_timelocked{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
        ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.7,'FaceColor',[229 194 152]/256,'unit','cm')
        ft_plot_headmodel(headmodels.headmodel_meg)
        ft_plot_headshape(headshape)
        hold off;
        title('SQUID-MEG')
        view([-140 10])
        savefig(h, fullfile(save_path, 'figs', 'squid_layout.fig'))
        saveas(h, fullfile(save_path, 'figs', 'squid_layout.jpg'))
        close all

        %% Save
        save(fullfile(save_path, [params.sub '_opm_timelockedT']), 'opm_timelockedT', '-v7.3');
        save(fullfile(save_path, [params.sub '_opmeeg_timelockedT']), 'opmeeg_timelockedT', '-v7.3');

    else
        disp('Required files not found. No transformed OPM data was saved.')
    end
end

%% Dipole fits
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);

    if exist(fullfile(save_path, 'dipoles.mat'),'file') && overwrite.dip==false
        load(fullfile(save_path, 'dipoles.mat'));
    else
        clear headmodels mri_resliced
        clear squidmag_timelocked squidgrad_timelocked squideeg_timelocked
        clear opm_timelockedT opmeeg_timelockedT
        load(fullfile(save_path, [params.sub '_opm_timelockedT.mat']))
        load(fullfile(save_path, [params.sub '_opmeeg_timelockedT.mat']))
        load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat']))
        squidmag_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_squidgrad_timelocked.mat']))
        squidgrad_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_squideeg_timelocked.mat']))
        squideeg_timelocked = timelocked;
        clear timelocked
        m100_latency = cell(5,1);
        for i_ph = 1:5
            m100_latency{i_ph} = [];
            load(fullfile(save_path, [params.sub '_opm_M100'])); 
            m100_latency{i_ph}.opm = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_opmeeg_M100'])); 
            m100_latency{i_ph}.opmeeg = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_squidmag_M100'])); 
            m100_latency{i_ph}.megmag = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_squidgrad_M100'])); 
            m100_latency{i_ph}.megplanar = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_squideeg_M100']));
            m100_latency{i_ph}.megeeg = M100{i_ph}.peak_latency;
        end
        load(fullfile(save_path, 'headmodels.mat'));
        load(fullfile(save_path, 'mri_resliced.mat'));
        [squidmag_dipole, squidgrad_dipole, opm_dipole, squideeg_dipole, opmeeg_dipole] = fit_dipoles(save_path,squidmag_timelocked,squidgrad_timelocked,squideeg_timelocked,opm_timelockedT,opmeeg_timelockedT,headmodels,mri_resliced,m100_latency,params);
    end
end

%% --- Group source level -------------------------------------------------
if ~exist(fullfile(base_save_path,'figs'), 'dir')
       mkdir(fullfile(base_save_path,'figs'))
end
subs = 2:13;
dipole_results_goup(base_save_path,subs, params)

%% Prepare MNE sourcemodel 
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);

    clear headmodels megmag_timelocked opm_timelockedT meshes filename
    clear headmodels meshes filename
    load(fullfile(save_path,'headmodels.mat'));
    load(fullfile(save_path,'meshes.mat'));    
    load(fullfile(save_path, 'mri_resliced.mat'));
    load(fullfile(save_path, [params.sub '_opm_timelockedT.mat']))
    load(fullfile(save_path, [params.sub '_opmeeg_timelockedT.mat']))
    load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat']))
    squidmag_timelocked = timelocked;
    load(fullfile(save_path, [params.sub '_squideeg_timelocked.mat']))
    squideeg_timelocked = timelocked;
    clear timelocked

    % Read and transform cortical restrained source model
    files = dir(fullfile(mri_path,'workbench'));
    if i_sub ==5
        files = dir(fullfile(save_path,'wb'));
    end
    for i = 1:length(files)
        if endsWith(files(i).name,'.L.midthickness.8k_fs_LR.surf.gii')
            filename = fullfile(mri_path,'workbench',files(i).name);
            %return;
        end
    end
    sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});
    T = mri_resliced.transform/mri_resliced.hdr.vox2ras;
    sourcemodelT = ft_transform_geometry(T, sourcemodel);
    sourcemodelT.inside = true(size(sourcemodelT.pos,1),1);
    %sourcemodelOPM = sourcemodelT;
    %sourcemodelOPM.pos = opm_trans.transformPointsForward(sourcemodelOPM.pos);

    % Plot source and head models
    h=figure; 
    ft_plot_mesh(sourcemodelT, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
    hold on; 
    ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
    ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
    ft_plot_sens(opm_timelockedT{1}.grad,'unit','cm')
    ft_plot_sens(opm_timelockedT{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
    hold off;
    title('OPM-MEG')
    view([-140 10])
    savefig(h, fullfile(save_path, 'figs', 'opm_layout2.fig'))
    saveas(h, fullfile(save_path, 'figs', 'opm_layout2.jpg'))

    h=figure; 
    ft_plot_mesh(sourcemodelT, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
    hold on; 
    ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
    ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
    ft_plot_sens(squidmag_timelocked{1}.grad,'unit','cm')
    ft_plot_sens(squideeg_timelocked{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
    hold off;
    title('SQUID-MEG')
    view([-140 10])
    savefig(h, fullfile(save_path, 'figs', 'meg_layout2.fig'))
    saveas(h, fullfile(save_path, 'figs', 'meg_layout2.jpg'))
    
    close all
    % Save
    save(fullfile(save_path, [params.sub '_sourcemodelT']), 'sourcemodelT', '-v7.3');
end

%% MNE
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);

    %% MNE fit
    if exist(fullfile(save_path, 'mne_fits.mat'),'file') && overwrite.mne==false
        load(fullfile(save_path, 'mne_fits.mat'));
    else
        clear headmodels sourcemodelT
        load(fullfile(save_path, [params.sub '_sourcemodelT']));
        load(fullfile(save_path,'headmodels.mat'));
        clear megmag_timelocked magplanar_timelocked squideeg_timelocked
        clear opm_timelockedT opmeeg_timelockedT
        load(fullfile(save_path, [params.sub '_opm_timelockedT.mat']))
        load(fullfile(save_path, [params.sub '_opmeeg_timelockedT.mat']))
        load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat']))
        squidmag_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_squidgrad_timelocked.mat']))
        squidgrad_timelocked = timelocked;
        load(fullfile(save_path, [params.sub '_squideeg_timelocked.mat']))
        squideeg_timelocked = timelocked;
        clear timelocked
        m100_latency = cell(5,1);
        for i_ph = 1:5
            m100_latency{i_ph} = [];
            load(fullfile(save_path, [params.sub '_opm_M100'])); 
            m100_latency{i_ph}.opm = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_opmeeg_M100'])); 
            m100_latency{i_ph}.opmeeg = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_squidmag_M100'])); 
            m100_latency{i_ph}.squidmag = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_squidgrad_M100'])); 
            m100_latency{i_ph}.squidgrad = M100{i_ph}.peak_latency;
            load(fullfile(save_path, [params.sub '_squideeg_M100']));
            m100_latency{i_ph}.squideeg = M100{i_ph}.peak_latency;
        end
        [megmag_mne, megplanaer_mne, opm_mne, eeg_mne] = fit_mne(save_path,squidmag_timelocked,squidgrad_timelocked,squideeg_timelocked,opm_timelockedT,opmeeg_timelockedT,headmodels,sourcemodelT,m100_latency,params);
    end
end
close all
