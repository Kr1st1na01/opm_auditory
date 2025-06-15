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
params.hpi_gof = 0.9;

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
       mkdir(save_path)
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

        % Dividing data for Evoked analysis och Frequency tag analysis
        params.ica = [opm_ica, opmeeg_ica, squid_ica, squideeg_ica];
        params.ica_labels = {'opm', 'opmeeg', 'squid', 'squideeg'};
        
        for i = 1:length(params.ica) 
            cfg = [];
            cfg.lpfilter  = 'yes';        % Apply lowpass filter
            cfg.lpfreq    = 20;      
            cfg.demean          = 'yes';
            cfg.baselinewindow  = [-0.100 0];
            cfg.toilim = [-0.100 0.5];
            
            % Evoked dataset
            filt_downsampled_data = ft_preprocessing(cfg, params.ica(i)); % Data is processed with a lowpass filter of 40 Hz
            evoked_data = ft_redefinetrial(cfg, filt_downsampled_data); % Trials are cropped into time of interest (toi)
            save(fullfile(save_path, [params.sub '_' params.ica_labels{i} '_evoked_ica']), 'evoked_data',"-v7.3");

            cfg = [];
            cfg.hpfilter  = 'yes';        
            cfg.hpfreq    = 30;      % Apply highpass filter
            cfg.demean          = 'yes';
            cfg.baselinewindow  = [-0.100 0];
            
            % TFR dataset
            frequency_tag_data = ft_preprocessing(cfg, params.ica(i)); % Data is processed with a lowpass filter of 30 Hz
            save(fullfile(save_path, [params.sub '_' params.ica_labels{i} '_frequency_tag_ica']), 'frequency_tag_data',"-v7.3");
        end
        close all
        clear i
        clear -regexp ^long ^cropped
    end
end

%% Loop over subjects 
for i_sub = 1:size(subses,1)
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    
    %% Paths
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
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
       
        %% Load data
        
        % Load data
        Evoked_opm = load(fullfile(save_path, [params.sub '_opm_evoked_ica.mat'])).evoked_data;
        Freqtag_opm = load(fullfile(save_path, [params.sub '_opm_frequency_tag_ica.mat'])).frequency_tag_data;

        Evoked_opmeeg = load(fullfile(save_path, [params.sub '_opmeeg_evoked_ica.mat'])).evoked_data;
        Freqtag_opmeeg = load(fullfile(save_path, [params.sub '_opmeeg_frequency_tag_ica.mat'])).frequency_tag_data;
    
        Evoked_squid = load(fullfile(save_path, [params.sub '_squid_evoked_ica.mat'])).evoked_data;
        Freqtag_squid = load(fullfile(save_path, [params.sub '_squid_frequency_tag_ica.mat'])).frequency_tag_data;

        Evoked_squideeg = load(fullfile(save_path, [params.sub '_squideeg_evoked_ica.mat'])).evoked_data;
        Freqtag_squideeg = load(fullfile(save_path, [params.sub '_squideeg_frequency_tag_ica.mat'])).frequency_tag_data;

        clear evoked_data frequency_tag_data 

        %% Average for OPM-MEG
        ylimit = {};
        
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        params.pow_label = 'Power [T^2]';
        %[opm_timelocked, peak, ylimit] = evoked_analysis(Evoked_opm, params, save_path, ylimit); % Timelockar vanlig MMN och plottar för Std, Low och High och kör freqanalysis på TFR
        peak = freqtag_analysis(Freqtag_opm, params, save_path, peak);
        save(fullfile(save_path, [params.sub '_opm_peaks']), 'peak' ,"-v7.3"); % Peak values saved for each modality
        close all
        
        load(fullfile(save_path, [params.sub '_opmeeg_layout.mat']))
        opmeeg_layout = layout;
        params.modality = 'opmeeg';
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        params.pow_label = 'Power [V^2]';        
        [opmeeg_timelocked, peak, ylimit] = evoked_analysis(Evoked_opmeeg, params, save_path, ylimit); 
        peak = freqtag_analysis(Freqtag_opmeeg, params, save_path, peak);
        save(fullfile(save_path, [params.sub '_opmeeg_peaks']), 'peak' ,"-v7.3");
        close all

        %% Average SQUID-MEG
        params.modality = 'squidmag';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        params.pow_label = 'Power [T^2]';
        [squidmag_timelocked, peak, ylimit] = evoked_analysis(Evoked_squid, params, save_path, ylimit); 
        peak = freqtag_analysis(Freqtag_squid, params, save_path, peak);
        save(fullfile(save_path, [params.sub '_squidmag_peaks']), 'peak' ,"-v7.3");
        close all

        params.modality = 'squidgrad';
        params.layout = 'neuromag306planar.lay';
        params.chs = 'megplanar';
        params.amp_scaler = 1e15/100;
        params.amp_label = 'B [fT/cm]';
        params.pow_label = 'Power [T^2/cm^2]';
        [squidgrad_timelocked, peak, ylimit] = evoked_analysis(Evoked_squid, params, save_path, ylimit); 
        peak = freqtag_analysis(Freqtag_squid, params, save_path, peak);
        save(fullfile(save_path, [params.sub '_squidgrad_peaks']), 'peak' ,"-v7.3");
        close all

        load(fullfile(save_path, [params.sub '_megeeg_layout.mat']))
        megeeg_layout = layout;
        params.modality = 'squideeg';
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        params.pow_label = 'Power [V^2]';
        [squideeg_timelocked, peak, ylimit] = evoked_analysis(Evoked_squideeg, params, save_path, ylimit); 
        peak = freqtag_analysis(Freqtag_squideeg, params, save_path, peak);
        save(fullfile(save_path, [params.sub '_squideeg_peaks']), 'peak' ,"-v7.3");
        close all

        clear -regexp ^TFR ^MMN ^opm ^opmeeg ^squid ^squidgrad ^squidmag ^squideeg layout peak ^Evoked_ ^Freqtag_

    end
% %%
%     params = rmfield(params,{'modality', 'layout', 'chs', 'amp_scaler', 'amp_label'}); % remove fields used for picking modality
%     create_bads_reports(base_save_path, i_sub, params);
%     close all
end

%% --- Statistical analysis -----------------------------------------------
stats.opm_val = [];
stats.opmeeg_val = [];
stats.squidmag_val = [];
stats.squidgrad_val = [];
stats.squideeg_val = [];

for i_sub = 1:size(subses,1)
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    
    % Paths
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);

    % Load the data and save them into one matrix each in struct stats
    % Each matrix will be 13 columns, one for each subject
    load(fullfile(save_path, [params.sub '_opm_peaks.mat']));
    stats.opm_val = [stats.opm_val peak.values];

    load(fullfile(save_path, [params.sub '_opmeeg_peaks.mat']));
    stats.opmeeg_val = [stats.opmeeg_val peak.values];

    load(fullfile(save_path, [params.sub '_squidmag_peaks.mat']));
    stats.squidmag_val = [stats.squidmag_val peak.values];

    load(fullfile(save_path, [params.sub '_squidgrad_peaks.mat']));
    stats.squidgrad_val = [stats.squidgrad_val peak.values];

    load(fullfile(save_path, [params.sub '_squideeg_peaks.mat'])).peak;
    stats.squideeg_val = [stats.squideeg_val peak.values];
end

if ~exist(fullfile(base_save_path,'statistics'), 'dir') % All statistics are saved in their own folder
   mkdir(fullfile(base_save_path,'statistics'))
end

stats.labels = peak.labels;
stats.labels = replace(erase(stats.labels, 'squideeg_'),'_', ', ');

clear peak

% Saving indexes and tables
i_amp = find(contains(stats.labels, 'amplitude'));
i_pow = find(contains(stats.labels, 'FFT'));
i_lat = find(contains(stats.labels, 'latency'));
% statistics_amp = table;
% statistics_pow = table;
% statistics_lat = table;
i_ap = {i_amp, i_pow}; % The indices of amplitudes and power

% First row of tables is the labels
statistics_amp = cell2table(stats.labels(i_amp), 'VariableNames',{'Labels'});
statistics_pow = cell2table(stats.labels(i_pow), 'VariableNames',{'Labels'});
statistics_lat = cell2table(stats.labels(i_lat), 'VariableNames',{'Labels'});

list = {stats.opm_val stats.squidgrad_val; stats.opm_val stats.squidmag_val; stats.squidmag_val stats.squidgrad_val; stats.opmeeg_val stats.squideeg_val}; % The matrices that will be compared are put beside each other in list
comparisons = {'opm vs squidgrad', 'opm vs squidmag', 'squidmag vs squidgrad', 'opmeeg vs squideeg'}; % The different comparisons
names = {'Amplitude', 'Power', 'Latency'; '#0072BD', '#D95319', '#77AC30'}; % For the histogram, the different comparisons and their colours

ratio = [];
diff = [];
log_ratio = [];
for i = 1:size(list,1)
    ratio = abs(list{i,1})./abs(list{i,2}); % The ratio for opm/squidgrad, opm/squidmag, squidmag/squidgrad, opmeeg/squideeg
    log_ratio = log10(ratio); % The ratios are logged
    diff = abs(list{i,1})-abs(list{i,2}); % The latency difference is calculated
      
    h = figure;
    for l = 1:numel(i_lat) % Latency plot
        histogram(diff(i_lat(l),:), 'FaceColor', names{2,end}); % Pick out the latency values based on the index in i_apl
        ylabel([names{1,end} ' Difference'])
    end
    sgtitle(h, comparisons{i}) % Big plot title
    saveas(h, fullfile(base_save_path, 'statistics', [comparisons{i} '_histogram_latency difference.jpg']))

    h = figure;
    for ap = 1:numel(i_ap) % Plot histogram, first subplot is ratio
        subplot(2,length(i_ap),ap)
        histogram(ratio(i_ap{ap},:), 'DisplayName', names{1,ap}, 'FaceColor', names{2,ap}); % Pick out the ratio values based on the index in i_apl
        title(names{1,ap})
        if ap==1
        ylabel('Ratio')
        end
       
        % Second subplot is log ratios
        subplot(2,length(i_ap),ap+length(i_ap))
        histogram(log_ratio(i_ap{ap},:), 'FaceColor', names{2,ap});
        if ap==1
        ylabel('Log ratio')
        end
    end
    sgtitle(h, comparisons{i}) % Big plot title
    saveas(h, fullfile(base_save_path, 'statistics', [comparisons{i} '_histogram_ratio vs log ratio' '.jpg']))

    save = {};
    for k = 1:size(i_amp) % Amplitudes
        [~, p] = ttest(abs(list{i,1}(i_amp(k),:)), abs(list{i,2}(i_amp(k),:))); % Calculating h and p-value
        save{k,1} = p;
        save{k,2} = mean(ratio(i_amp(k),:));
        save{k,3} = mean(log_ratio(i_amp(k),:))*20;
        [h, p, ~, stats] = ttest(log_ratio(i_amp(k),:)); % Calculating h and p-value for ratio
        save{k,4} = [h,p];
        save{k,5} = stats.sd;
    end
    statistics_pre = cell2table(save); % Putting the values in the table
    statistics_pre.Properties.VariableNames = {comparisons{i}, ['mean ratio ' int2str(i)], ['mean log ratio ' int2str(i) ' (dB)'], ['ttest on log ratio ' int2str(i)], ['standard deviation' int2str(i)]};
    statistics_amp = [statistics_amp statistics_pre]; % Adding the table to the big table

    save = {};
    for k = 1:size(i_pow) % Power
        [~, p] = ttest(abs(list{i,1}(i_pow(k),:)), abs(list{i,2}(i_pow(k),:))); % Calculating h and p-value
        save{k,1} = p;
        save{k,2} = mean(ratio(i_pow(k),:));
        save{k,3} = mean(log_ratio(i_pow(k),:))*10;
        [h, p, ~, stats] = ttest(log_ratio(i_pow(k),:)); % Calculating h and p-value for ratio
        save{k,4} = [h,p];
        save{k,5} = stats.sd;
    end
    statistics_pre = cell2table(save); % Putting the values in the table
    statistics_pre.Properties.VariableNames = {comparisons{i}, ['mean ratio ' int2str(i)], ['mean log ratio ' int2str(i) ' (dB)'], ['ttest on log ratio' int2str(i)], ['standard deviation' int2str(i)]}; % Changing title
    statistics_pow = [statistics_pow statistics_pre]; % Adding the pre-table to the table
    save = {};

    for k = 1:size(i_lat) % Latency
        [~, p] = ttest(list{i,1}(i_lat(k),:), list{i,2}(i_lat(k),:)); % Calculating h and p-value
        save{k,1} = p;
        save{k,2} = mean(diff(i_lat(k),:));
        [h, p] = ttest(diff(i_lat(k),:)); % Calculating h and p-value for ratio
        save{k,3} = [h,p];
        save{k,4} = stats.sd;
    end
    statistics_pre = cell2table(save); % Putting the values in the table
    statistics_pre.Properties.VariableNames = {comparisons{i}, ['mean difference ' int2str(i)], ['ttest on difference ' int2str(i)], ['standard deviation' int2str(i)]};
    statistics_lat = [statistics_lat statistics_pre]; % Adding the table to the big table
end  

% Row names, when saved as cvs or xlsx it does not show up properly so they
% are added to the table from cells
% statistics_amp.Properties.RowNames = replace(stats.lab(i_amp),'_', ', ');
% statistics_pow.Properties.RowNames = replace(stats.lab(i_pow),'_', ', ');
% statistichrcs_lat.Properties.RowNames = replace(stats.lab(i_lat),'_', ', ');

writetable(statistics_amp, fullfile(base_save_path, 'statistics', 'table_amplitudes.csv')) % The tables are saved
writetable(statistics_pow, fullfile(base_save_path, 'statistics', 'table_power.csv'))
writetable(statistics_lat, fullfile(base_save_path, 'statistics', 'table_latency.csv'))

close all
clear k save ratio log_ratio list i_amp i_pow h l names p diff comparisons i ap
clear -regexp i_ statistics_

%% --- Group sensor level -------------------------------------------------

%% Prepare MRIs
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}],'mri');
    if exist(fullfile(save_path, 'headmodels.mat'),'file') && overwrite.mri==false
        disp(['Not overwriting MRI for ' params.sub]);
    else
        meg_file = fullfile(raw_path, 'meg', 'AudOddMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        mri_file = fullfile(mri_path, 'orig','001.mgz');
        prepare_mri(mri_file,meg_file,save_path);
        close all
    end
end

%% HPI localization
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    hpi_path = fullfile(raw_path, 'osmeg'); %hpi_file = fullfile(raw_path, 'osmeg', 'HPIpre_raw.fif');

    if exist(fullfile(save_path, 'opm_trans.mat'),'file') && overwrite.coreg==false
        disp(['Not overwriting OPM transform for ' params.sub]);
    else
        meg_file = fullfile(raw_path, 'meg', 'AudOddMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        ft_hastoolbox('mne', 1);
        load(fullfile(save_path, [params.sub '_opm_long_ica']));
        params.include_chs = TFR_data.label(find(contains(TFR_data.label,'bz')));
        fit_hpi(hpi_path, meg_file, save_path, params);
        close all
    end
end

%% Transform for OPM data
for i_sub = 2:size(subses,1) % Skip first subject
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    save_path = fullfile(base_save_path, params.sub);
    shared_folder_path = fullfile('/home/share/opm_auditory', params.sub);
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    if ~exist(fullfile(save_path,'source analysis'), 'dir')
       mkdir(fullfile(save_path,'source analysis'))
    end

    if exist(fullfile(save_path, 'headmodels.mat'),'file') && exist(fullfile(save_path, 'opm_trans.mat'),'file') && or(overwrite.preproc,or(overwrite.coreg,overwrite.mri))
        clear headmodels meshes filename headshape 
        headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
        meshes = load(fullfile(save_path,'meshes.mat')).meshes;    
        mri_resliced = load(fullfile(save_path, 'mri_resliced.mat')).mri_resliced;
        opm_trans = load(fullfile(save_path, 'opm_trans.mat')).opm_trans;
        
        clear opm_timelockedT opmeeg_timelcokedT 
        clear squideeg_timelocked squidmag_timelocked
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelocked_Std.mat'])).timelocked_data;
        opmeeg_timelockedT = load(fullfile(save_path, [params.sub '_opmeeg_timelocked_Std.mat'])).timelocked_data;
        squideeg_timelocked = load(fullfile(save_path, [params.sub '_squideeg_timelocked_Std.mat'])).timelocked_data;
        squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked_Std.mat'])).timelocked_data;

        opm_timelockedT.grad = ft_convert_units(opm_timelockedT.grad,'cm');

        % Transform opm & opmeeg data 
        meg_file = (fullfile(raw_path, 'meg', 'AudOddMEG_proc-tsss+corr98+mc+avgHead_meg.fif'));
        headshape = ft_read_headshape(meg_file);

        opm_timelockedT.grad.chanpos = opm_trans.transformPointsForward(opm_timelockedT.grad.chanpos);
        opm_timelockedT.grad.coilpos = opm_trans.transformPointsForward(opm_timelockedT.grad.coilpos);
        opm_timelockedT.grad.chanori = (opm_trans.Rotation'*opm_timelockedT.grad.chanori')';
        opm_timelockedT.grad.coilori = (opm_trans.Rotation'*opm_timelockedT.grad.coilori')';
        opmeeg_timelockedT.elec.chanpos = squideeg_timelocked.elec.chanpos;
        opmeeg_timelockedT.elec.elecpos = squideeg_timelocked.elec.elecpos;

        % Read and transform cortical restrained source model
        clear sourcemodel sourcemodel_inflated
        files = dir(fullfile(mri_path,'workbench'));
        for i = 1:length(files)
            if endsWith(files(i).name,'.L.midthickness.8k_fs_LR.surf.gii')
                filename = fullfile(mri_path,'workbench',files(i).name);
            elseif endsWith(files(i).name,'.L.aparc.8k_fs_LR.label.gii')
                filename2 = fullfile(mri_path,'workbench',files(i).name);
            end
        end
        sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});

        aparc_L = ft_read_atlas({filename2,filename});
        aparc_R = ft_read_atlas({strrep(filename2,'.L.','.R.'),strrep(filename,'.L.','.R.')});
        tmp = ft_read_atlas(strrep(filename2, '.L.', '.R.'),'format','caret_label');
        n_labels = length(aparc_L.parcellationlabel);
        atlas = [];
        atlas.parcellationlabel = [aparc_L.parcellationlabel; aparc_R.parcellationlabel];
        atlas.parcellation = [aparc_L.parcellation; aparc_R.parcellation + n_labels];
        atlas.rgba = [aparc_L.rgba; aparc_R.rgba; [0 0 0 1]];
        n_labels = length(atlas.parcellationlabel);
        atlas.parcellation(isnan(atlas.parcellation))=n_labels+1;
        sourcemodel.brainstructure = atlas.parcellation;
        sourcemodel.brainstructurelabel = atlas.parcellationlabel;
        sourcemodel.brainstructurecolor = atlas.rgba;
        clear atlas aparc_L aparc_R

        T = mri_resliced.transform/mri_resliced.hdr.vox2ras;
        sourcemodel = ft_transform_geometry(T, sourcemodel);
        sourcemodel.inside = true(size(sourcemodel.pos,1),1);

        for i = 1:length(files)
            if endsWith(files(i).name,'.L.inflated.8k_fs_LR.surf.gii')
                filename = fullfile(mri_path,'workbench',files(i).name);
            end
        end
        sourcemodel_inflated = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});
        sourcemodel_inflated = ft_transform_geometry(T, sourcemodel_inflated);
        sourcemodel_inflated.inside = true(size(sourcemodel_inflated.pos,1),1);
        sourcemodel_inflated.brainstructure = sourcemodel.brainstructure;
        sourcemodel_inflated.brainstructurelabel = sourcemodel.brainstructurelabel;
        sourcemodel_inflated.brainstructurecolor = sourcemodel.brainstructurecolor;

        % Plot source and head models
        h=figure; 
        ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
        hold on; 
        ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
        ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
        ft_plot_sens(opm_timelockedT.grad,'unit','cm')
        ft_plot_sens(opmeeg_timelockedT.elec,'unit','cm', 'style', '.r','elecsize',20)
        hold off;
        title('OPM-MEG')
        view([-140 10])
        savefig(h, fullfile(save_path, 'source analysis', 'opm_layout2.fig'))
        saveas(h, fullfile(save_path, 'source analysis', 'opm_layout2.jpg'))
    
        h=figure; 
        ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
        hold on; 
        ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
        ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
        ft_plot_sens(squidmag_timelocked.grad,'unit','cm')
        ft_plot_sens(squideeg_timelocked.elec,'unit','cm', 'style', '.r','elecsize',20)
        hold off;
        title('SQUID-MEG')
        view([-140 10])
        savefig(h, fullfile(save_path, 'source analysis', 'meg_layout2.fig'))
        saveas(h, fullfile(save_path, 'source analysis', 'meg_layout2.jpg'))

        close all

        %% Save
        save(fullfile(save_path, [params.sub '_opm_timelockedT']), 'opm_timelockedT', '-v7.3');
        save(fullfile(save_path, [params.sub '_opmeeg_timelockedT']), 'opmeeg_timelockedT', '-v7.3');
        save(fullfile(save_path, [params.sub '_sourcemodel']), 'sourcemodel', '-v7.3');
        save(fullfile(save_path, [params.sub '_sourcemodel_inflated']), 'sourcemodel_inflated', '-v7.3');
        
        clear opm_timelockedT opmeeg_timelockedT sourcemodel_inflated sourcemodel
    else
        disp(['Required files not found. No transformed OPM/sourcemodel data was saved for ' params.sub])
    end
end


%% Load files with peak values

stats.opm_val = [];
stats.opmeeg_val = [];
stats.squidmag_val = [];
stats.squidgrad_val = [];
stats.squideeg_val = [];

if exist('stats', 'var')
    disp('stats is in the workspace')
else
    for i_sub = 1:size(subses,1)
        params.sub = ['sub_' num2str(i_sub,'%02d')];
        
        % Paths
        raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
        save_path = fullfile(base_save_path,params.sub);
    
        % Load the data and save them into one matrix each in struct stats
        % Each matrix will be 13 columns, one for each subject
        load(fullfile(save_path, [params.sub '_opm_peaks.mat']));
        stats.opm_val = [stats.opm_val peak.values];
    
        load(fullfile(save_path, [params.sub '_opmeeg_peaks.mat']));
        stats.opmeeg_val = [stats.opmeeg_val peak.values];
    
        load(fullfile(save_path, [params.sub '_squidmag_peaks.mat']));
        stats.squidmag_val = [stats.squidmag_val peak.values];
    
        load(fullfile(save_path, [params.sub '_squidgrad_peaks.mat']));
        stats.squidgrad_val = [stats.squidgrad_val peak.values];
    
        load(fullfile(save_path, [params.sub '_squideeg_peaks.mat'])).peak;
        stats.squideeg_val = [stats.squideeg_val peak.values];
    end
end

if ~exist(fullfile(base_save_path,'statistics'), 'dir') % All statistics are saved in their own folder
   mkdir(fullfile(base_save_path,'statistics'))
end

stats.labels = peak.labels;
stats.labels = replace(erase(stats.labels, 'squideeg_'),'_', ', ');

clear peak

    %% Dipole fits
for i_sub = 2:size(subses,1) % Skip first subject
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    save_path = fullfile(base_save_path, params.sub);
    raw_path = fullfile(base_data_path,'MEG',['NatMEG' subses{i_sub,1}], subses{i_sub,2});
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);    
    if exist(fullfile(save_path, 'source analysis', 'dipoles.mat'),'file') && overwrite.dip==false
        disp(['Not overwriting dipole source reconstruction for ' params.sub]);
        else
            clear headmodels mri_resliced
            headmodels = load(fullfile(save_path, 'headmodels.mat')).headmodels;
            mri_resliced = load(fullfile(save_path, 'mri_resliced.mat')).mri_resliced;
            
            clear squimdag_timelocked squidgrad_timelocked opm_timelockedT
            opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
            opmeeg_timelockedT = load(fullfile(save_path, [params.sub '_opmeeg_timelockedT.mat'])).opmeeg_timelockedT;
            squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked_Std.mat'])).timelocked_data;
            squidgrad_timelocked = load(fullfile(save_path, [params.sub '_squidgrad_timelocked_Std.mat'])).timelocked_data;
            squideeg_timelocked = load(fullfile(save_path, [params.sub '_squideeg_timelocked_Std.mat'])).timelocked_data;

            clear M100_opm M100_squidmag M100_squidgrad
            row_idx = find(contains(stats.lab, 'M100, Std')); % Row index for latency
            [squidmag_dipole, squidgrad_dipole, opm_dipole] = fit_dipoles(save_path, squidmag_timelocked, squidgrad_timelocked, squideeg_timelocked, opm_timelockedT, opmeeg_timelockedT, headmodels, mri_resliced, params, stats, row_idx);
    end
end
%%
% Dipole group analysis
subs = 2:13;
dipole_results_group(base_save_path,subs, params)

%% Prepare Leadfields
for i_sub = 2:size(subses,1)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    if exist(fullfile(save_path, 'source analysis', [params.sub '_squidgrad_leadfield.mat']),'file') %&& overwrite.dip==false
        disp(['Not overwriting dipole source reconstruction for ' params.sub]);
   % Load files
    else
        clear headmodels sourcemodel sourcemodel_inflated
        sourcemodel = load(fullfile(save_path, 'source analysis', [params.sub '_sourcemodel.mat'])).sourcemodel;
        headmodels = load(fullfile(save_path, 'headmodels.mat')).headmodels;
        
        clear squimdag_timelocked squidgrad_timelocked opm_timelockedT
        opm_timelockedT = load(fullfile(save_path, 'source analysis', [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
        squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked_Std.mat'])).timelocked_data;
        squidgrad_timelocked = load(fullfile(save_path, [params.sub '_squidgrad_timelocked_Std.mat'])).timelocked_data;
        
        % Prepare leadfields
        prep_leadfield(headmodels,sourcemodel, squidmag_timelocked, [params.sub '_squidmag'], save_path)
        prep_leadfield(headmodels,sourcemodel, squidgrad_timelocked, [params.sub '_squidgrad'], save_path)
        prep_leadfield(headmodels,sourcemodel, opm_timelockedT, [params.sub '_opm'], save_path)
    end


%% MNE fit (did not finish)
    if exist(fullfile(save_path, 'mne_fits.mat'),'file') && overwrite.mne==false
        disp(['Not overwriting MNE source reconstruction for ' params.sub]);
    else
        %LOAD LEADFIELDS AND SOURCEMODELS
        clear headmodels sourcemodel sourcemodel_inflated
        sourcemodel = load(fullfile(save_path, [params.sub '_sourcemodel.mat'])).sourcemodel;
        sourcemodel_inflated = load(fullfile(save_path, [params.sub '_sourcemodel_inflated.mat'])).sourcemodel_inflated;
        headmodels = load(fullfile(save_path, 'headmodels.mat')).headmodels;
        
        clear squimdag_timelocked squidgrad_timelocked opm_timelockedT
        opm_timelockedT = load(fullfile(save_path, 'source analysis', [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
        squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked_Std.mat'])).timelocked_data;
        squidgrad_timelocked = load(fullfile(save_path, [params.sub '_squidgrad_timelocked_Std.mat'])).timelocked_data;
        
        squidmag_timelocked.label = squidmag_timelocked.label(contains(squidmag_timelocked.label, 'MEG'));
        squidgrad_timelocked.label = squidgrad_timelocked.label(contains(squidgrad_timelocked.label, 'MEG'));
        
        leadfield = load(fullfile(save_path, [params.sub '_squidmag_leadfield.mat'])).leadfield;
        params.modality = 'squidmag';
        params.scalesourcecov = 'no';
        row_idx = find(contains(stats.labels, 'M100, Std')); % Row index for latency 
        val = stats.squidmag_val(row_idx,i_sub);
        fit_mne(save_path, squidmag_timelocked, headmodels, sourcemodel, sourcemodel_inflated, leadfield, params, val);
        
        leadfield = load(fullfile(save_path, [params.sub '_squidgrad_leadfield.mat'])).leadfield;
        params.modality = 'squidgrad';
        params.scalesourcecov = 'no';
        row_idx = find(contains(stats.labels, 'M100, Std')); % Row index for latency
        val = stats.squidgrad_val(row_idx,i_sub);
        fit_mne(save_path, squidgrad_timelocked, headmodels, sourcemodel, sourcemodel_inflated, leadfield, params, val);
               
        leadfield = load(fullfile(save_path, [params.sub '_opm_leadfield.mat'])).leadfield;
        params.modality = 'opm';
        params.scalesourcecov = 'yes';
        row_idx = find(contains(stats.labels, 'M100, Std')); % Row index for latency
        val = stats.opm_val(row_idx,i_sub);
        fit_mne(save_path, opm_timelockedT, headmodels, sourcemodel, sourcemodel_inflated, leadfield, params, val);

        clear leadfield
    end
end
close all

%% MNE group analysis (använd som inspo, ändra)
subs = [2:10 12:13];
mne_results_goup(base_save_path, subs, params);
