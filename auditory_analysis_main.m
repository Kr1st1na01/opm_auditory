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

%% Dividing data for MMN och TFR
        % MMN data is generated first (= cropped data) and TFR data is generated after
        
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

            long_data = ft_preprocessing(cfg, params.ica(i)); % Data is processed with a lowpass filter of 30 Hz
            save(fullfile(save_path, [params.sub '_' params.ica_labels{i} '_long_ica']), 'long_data',"-v7.3");
        end
        close all
        clear i
        clear -regexp ^long ^cropped
    end
end

%% Loop over subjects 

for i_sub = 2%:size(subses,1)
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
        Evoked_opm = load(fullfile(save_path, [params.sub '_opm_cropped_ica.mat'])).cropped_data;
        Freqtag_opm = load(fullfile(save_path, [params.sub '_opm_long_ica.mat'])).long_data;

        Evoked_opmeeg = load(fullfile(save_path, [params.sub '_opmeeg_cropped_ica.mat'])).cropped_data;
        Freqtag_opmeeg = load(fullfile(save_path, [params.sub '_opmeeg_long_ica.mat'])).long_data;
    
        Evoked_squid = load(fullfile(save_path, [params.sub '_squid_cropped_ica.mat'])).cropped_data;
        Freqtag_squid = load(fullfile(save_path, [params.sub '_squid_long_ica.mat'])).long_data;

        Evoked_squideeg = load(fullfile(save_path, [params.sub '_squideeg_cropped_ica.mat'])).cropped_data;
        Freqtag_squideeg = load(fullfile(save_path, [params.sub '_squideeg_long_ica.mat'])).long_data;

        clear cropped_data long_data labels 

        %% Average for OPM-MEG
        peak = [];
        peak.values = [];
        peak.labels = {};
        ylimit = {};
        
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        params.pow_label = 'Power [T^2]';
        [opm_timelocked, peak, ylimit] = evoked_analysis(Evoked_opm, params, save_path, peak, ylimit); % Timelockar vanlig MMN och plottar för Std, Low och High och kör freqanalysis på TFR
        peak = freqtag_analysis(Freqtag_opm, params, save_path, peak);
        close all
        
        load(fullfile(save_path, [params.sub '_opmeeg_layout.mat']))
        opmeeg_layout = layout;
        params.modality = 'opmeeg';
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        params.pow_label = 'Power [V^2]';        
        %[opmeeg_timelocked, peak, ylimit] = evoked_analysis(Evoked_opmeeg, params, save_path, peak, ylimit); 
        %peak = freqtag_analysis(Freqtag_opmeeg, params, save_path, peak);
        close all

        %% Average SQUID-MEG
        params.modality = 'squidmag';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        params.pow_label = 'Power [T^2]';
        [squidmag_timelocked, peak, ylimit] = evoked_analysis(Evoked_squid, params, save_path, peak, ylimit); 
        peak = freqtag_analysis(Freqtag_squid, params, save_path, peak);
        close all

        params.modality = 'squidgrad';
        params.layout = 'neuromag306planar.lay';
        params.chs = 'megplanar';
        params.amp_scaler = 1e15/100;
        params.amp_label = 'B [fT/cm]';
        params.pow_label = 'Power [T^2/cm^2]';
        [squidgrad_timelocked, peak, ylimit] = evoked_analysis(Evoked_squid, params, save_path, peak, ylimit); 
        peak = freqtag_analysis(Freqtag_squid, params, save_path, peak);
        close all

        load(fullfile(save_path, [params.sub '_megeeg_layout.mat']))
        megeeg_layout = layout;
        params.modality = 'squideeg';
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        params.pow_label = 'Power [V^2]';
        [squideeg_timelocked, peak, ylimit] = evoked_analysis(Evoked_squideeg, params, save_path, peak, ylimit); 
        peak = freqtag_analysis(Freqtag_squideeg, params, save_path, peak);
        
        save(fullfile(save_path, [params.sub '_peaks']), 'peak' ,"-v7.3");
        close all
        clear -regexp ^TFR ^MMN ^opm ^opmeeg ^squid ^squidgrad ^squidmag ^squideeg layout peak ^Evoked_ ^Freqtag_

    end
% %%
%     params = rmfield(params,{'modality', 'layout', 'chs', 'amp_scaler', 'amp_label'}); % remove fields used for picking modality
%     create_bads_reports(base_save_path, i_sub, params);
%     close all
end

%% --- Statistical analysis -----------------------------------------------
stats.values = [];

for i_sub = 1:size(subses,1)
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    
    % Paths
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);

    % Load the data
    load(fullfile(save_path, [params.sub '_peaks.mat']))
    stats.values = [stats.values peak.values];
end

if ~exist(fullfile(base_save_path,'statistics'), 'dir')
   mkdir(fullfile(base_save_path,'statistics'))
end

stats.labels = peak.labels;
clear peak

stats.val_opm = stats.values(startsWith(stats.labels, 'opm_'), :);
stats.val_opmeeg = stats.values(startsWith(stats.labels, 'opmeeg'), :);
stats.val_squidmag = stats.values(startsWith(stats.labels, 'squidmag'), :);
stats.val_squidgrad = stats.values(startsWith(stats.labels, 'squidgrad'), :);
stats.val_squideeg = stats.values(startsWith(stats.labels, 'squideeg'), :);
stats.lab = replace(erase(stats.labels(startsWith(stats.labels, 'opm_')), 'opm_'),'_', ', ');

% Saving indexes and tables
i_amp = find(contains(stats.lab, 'amplitude'));
i_pow = find(contains(stats.lab, 'FFT'));
i_lat = find(contains(stats.lab, 'latency'));
% statistics_amp = table;
% statistics_pow = table;
% statistics_lat = table;
i_ap = {i_amp, i_pow};

% First row of tables is the labels
statistics_amp = cell2table(stats.lab(i_amp), 'VariableNames',{'Labels'});
statistics_pow = cell2table(stats.lab(i_pow), 'VariableNames',{'Labels'});
statistics_lat = cell2table(stats.lab(i_lat), 'VariableNames',{'Labels'});

list = {stats.val_opm stats.val_squidgrad; stats.val_opm stats.val_squidmag; stats.val_squidmag stats.val_squidgrad; stats.val_opmeeg stats.val_squideeg};
comparisons = {'opm vs squidgrad', 'opm vs squidmag', 'squidmag vs squidgrad', 'opmeeg vs squideeg'};
names = {'Amplitude', 'Power', 'Latency'; '#0072BD', '#D95319', '#77AC30'};

ratio = [];
diff = [];
log_ratio = [];
for i = 1:size(list,1)
    ratio = abs(list{i,1})./abs(list{i,2}); % Now the list will contain one comparision at the moment, opm/squidgrad, opm/squidmag, squidmag/squidgrad, opmeeg/squideeg
    log_ratio = log10(ratio);
    diff = abs(list{i,1})-abs(list{i,2});
      
    h = figure;
    for l = 1:numel(i_lat) % Latency plot
        histogram(diff(i_lat(l),:), 'FaceColor', names{2,end}); % Pick out the latency values based on the index in i_apl
        ylabel([names{1,end} ' Difference'])
    end
    sgtitle(h, comparisons{i}) % Big plot title
    saveas(h, fullfile(base_save_path, 'statistics', [comparisons{i} '_histogram_latency difference.jpg']))

    h = figure;
    for ap = 1:numel(i_ap)    
        subplot(2,length(i_ap),ap) % Plot histogram
        histogram(ratio(i_ap{ap},:), 'DisplayName', names{1,ap}, 'FaceColor', names{2,ap}); % Pick out the ratio values based on the index in i_apl
        title(names{1,ap})
        if ap==1
        ylabel('Ratio')
        end
       
        % Plotting all log ratios       
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
    statistics_pre.Properties.VariableNames = {comparisons{i}, ['mean ratio ' int2str(i)], ['mean log ratio ' int2str(i) ' (dB)'], ['ttest on log ratio' int2str(i)], ['standard deviation' int2str(i)]};
    statistics_pow = [statistics_pow statistics_pre]; % Adding the table to the big table
    save = {};

    for k = 1:size(i_lat) % Latency
        [~, p] = ttest(abs(list{i,1}(i_lat(k),:)), abs(list{i,2}(i_lat(k),:))); % Calculating h and p-value
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

% Row names but when saved as cvs or xlsx it does not show
% statistics_amp.Properties.RowNames = replace(stats.lab(i_amp),'_', ', ');
% statistics_pow.Properties.RowNames = replace(stats.lab(i_pow),'_', ', ');
% statistichrcs_lat.Properties.RowNames = replace(stats.lab(i_lat),'_', ', ');

writetable(statistics_amp, fullfile(base_save_path, 'statistics', 'table_amplitudes.csv'))
writetable(statistics_pow, fullfile(base_save_path, 'statistics', 'table_power.csv'))
writetable(statistics_lat, fullfile(base_save_path, 'statistics', 'table_latency.csv'))
% 
% h = figure;
% x1=table2array(statistics_amp(:,9));
% x2=table2array(statistics_lat(:,7));
% x3=table2array(statistics_pow(:,9));
% matris = [x3];
% boxplot(matris)
% hold on
% errorbar(table2array(statistics_amp(:,11)))
% legend
% hold off


close all
clear k save ratio log_ratio list i_amp i_pow h l names p stats diff comparisons i ap
clear -regexp i_ statistics_

%% --- Group sensor level -------------------------------------------------
% Load files with peak values
stats.values = [];

for i_sub = 1:size(subses,1)
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    
    % Paths
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);

    % Load the data
    load(fullfile(save_path, [params.sub '_peaks.mat']))
    stats.values = [stats.values peak.values];
end

stats.labels = peak.labels;
clear peak

stats.lab = replace(erase(stats.labels(startsWith(stats.labels, 'opm_')), 'opm_'),'_', ', ');
i_lat = find(contains(stats.lab, 'latency'));
stats.lab = stats.lab(i_lat);

stats.val_opm = stats.values(startsWith(stats.labels, 'opm_'), :);
stats.val_opm = stats.val_opm(i_lat,:);

stats.val_opmeeg = stats.values(startsWith(stats.labels, 'opmeeg'), :);
stats.val_opmeeg = stats.val_opmeeg(i_lat,:);

stats.val_squidmag = stats.values(startsWith(stats.labels, 'squidmag'), :);
stats.val_squidmag = stats.val_squidmag(i_lat,:);

stats.val_squidgrad = stats.values(startsWith(stats.labels, 'squidgrad'), :);
stats.val_squidgrad = stats.val_squidgrad(i_lat,:);

stats.val_squideeg = stats.values(startsWith(stats.labels, 'squideeg'), :);
stats.val_squideeg = stats.val_squideeg(i_lat,:);


% Aligning head- and sourcemodel
for i_sub = 2 %2:size(subses,1) % Skip first subject
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    save_path = fullfile(base_save_path, params.sub);
    shared_folder_path = fullfile('/home/share/opm_auditory', params.sub);
    raw_path = fullfile(base_data_path,'MEG',['NatMEG' subses{i_sub,1}], subses{i_sub,2});
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    if ~exist(fullfile(save_path,'source analysis'), 'dir')
       mkdir(fullfile(save_path, params.sub,'source analysis'))
    end

    % Load Headmodels
    headmodels = load(fullfile(shared_folder_path, 'headmodels.mat')).headmodels;

    % Load sourcemodels
    sourcemodel = load(fullfile(shared_folder_path, [params.sub '_sourcemodel.mat'])).sourcemodel;
    close all

    % Transform for OPM data
    meshes = load(fullfile(shared_folder_path,'meshes.mat')).meshes;    
    mri_resliced = load(fullfile(shared_folder_path, 'mri_resliced_cm.mat')).mri_resliced_cm;
    opm_trans = load(fullfile(shared_folder_path, 'opm_trans.mat')).opm_trans;
     
    opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelocked_Std.mat'])).timelocked_data;
    opmeeg_timelockedT = load(fullfile(save_path, [params.sub '_opmeeg_timelocked_Std.mat'])).timelocked_data;
    squideeg_timelocked = load(fullfile(save_path, [params.sub '_squideeg_timelocked_Std.mat'])).timelocked_data;
    squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked_Std.mat'])).timelocked_data;
    
    opm_timelockedT.grad = ft_convert_units(opm_timelockedT.grad,'cm');

    % Transform opm & opmeeg data 
    opm_timelockedT.grad.chanpos = opm_trans.transformPointsForward(opm_timelockedT.grad.chanpos);
    opm_timelockedT.grad.coilpos = opm_trans.transformPointsForward(opm_timelockedT.grad.coilpos);
    opm_timelockedT.grad.chanori = (opm_trans.Rotation'*opm_timelockedT.grad.chanori')';
    opm_timelockedT.grad.coilori = (opm_trans.Rotation'*opm_timelockedT.grad.coilori')';
    opmeeg_timelockedT.elec.chanpos = squideeg_timelocked.elec.chanpos;
    opmeeg_timelockedT.elec.elecpos = squideeg_timelocked.elec.elecpos;

    % Plot aligned source- and head models
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

    % MEG
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
    save(fullfile(save_path, 'source analysis', [params.sub '_opm_timelockedT']), 'opm_timelockedT', '-v7.3');
    save(fullfile(save_path, 'source analysis', [params.sub '_opmeeg_timelockedT']), 'opmeeg_timelockedT', '-v7.3');
    save(fullfile(save_path, 'source analysis', [params.sub '_sourcemodel']), 'sourcemodel', '-v7.3');
    %save(fullfile(save_path, [params.sub '_sourcemodel_inflated']), 'sourcemodel_inflated', '-v7.3');

end

            % What do I want to look at? MMN,ososv
    %         list = ['Std','MMN'];
    %         for i = 1:length(list)
    %             ind = find(contains(list(i)));
    %         end
    %% Dipole fits
for i_sub = 2 %2:size(subses,1) % Skip first subject
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    save_path = fullfile(base_save_path, params.sub);
    shared_folder_path = fullfile('/home/share/opm_auditory', params.sub);
    raw_path = fullfile(base_data_path,'MEG',['NatMEG' subses{i_sub,1}], subses{i_sub,2});
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);    
    if exist(fullfile(save_path, 'source analysis', 'dipoles.mat'),'file') && overwrite.dip==false
        disp(['Not overwriting dipole source reconstruction for ' params.sub]);
        else
            clear headmodels mri_resliced
            headmodels = load(fullfile(shared_folder_path, 'headmodels.mat')).headmodels;
            mri_resliced = load(fullfile(shared_folder_path, 'mri_resliced_cm.mat')).mri_resliced_cm;
            
            clear squimdag_timelocked squidgrad_timelocked opm_timelockedT
            opm_timelockedT = load(fullfile(save_path, 'source analysis', [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
            opmeeg_timelockedT = load(fullfile(save_path, 'source analysis', [params.sub '_opmeeg_timelockedT.mat'])).opmeeg_timelockedT;
            squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked_Std.mat'])).timelocked_data;
            squidgrad_timelocked = load(fullfile(save_path, [params.sub '_squidgrad_timelocked_Std.mat'])).timelocked_data;
            squideeg_timelocked = load(fullfile(save_path, [params.sub '_squideeg_timelocked_Std.mat'])).timelocked_data;

            clear M100_opm M100_squidmag M100_squidgrad
            M100_opm = load(fullfile(save_path, [params.sub '_opm_Evoked data_M300_MMN_max sensor.mat'])).timelocked; 
            M100_squidmag = load(fullfile(save_path, [params.sub '_squidmag_Evoked data_M100_Std_max sensor.mat'])).timelocked; 
            M100_squidgrad = load(fullfile(save_path, [params.sub '_squidgrad_Evoked data_M100_Std_max sensor.mat'])).timelocked;
            row_idx = find(contains(stats.lab, 'MMN')); % Row index for latency
            [squidmag_dipole, squidgrad_dipole, opm_dipole] = fit_dipoles(save_path, squidmag_timelocked, squidgrad_timelocked, squideeg_timelocked, opm_timelockedT, opmeeg_timelockedT, headmodels, mri_resliced, M100_squidmag, MMN_squidgrad, MMN_opm, params, stats, row_idx);
    end
end
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
    if exist(fullfile(save_path, 'source analysis', [params.sub '_opm_leadfield.mat']),'file') && overwrite.dip==false
        disp(['Not overwriting dipole source reconstruction for ' params.sub]);
   % Load files
    else
        clear headmodels sourcemodel sourcemodel_inflated
        sourcemodel = load(fullfile(shared_folder_path, [params.sub '_sourcemodel.mat'])).sourcemodel;
        sourcemodel_inflated = load(fullfile(shared_folder_path, [params.sub '_sourcemodel_inflated.mat'])).sourcemodel_inflated;
        headmodels = load(fullfile(shared_folder_path, 'headmodels.mat')).headmodels;
        
        clear squimdag_timelocked squidgrad_timelocked opm_timelockedT
        opm_timelockedT = load(fullfile(save_path, 'source analysis', [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
        squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked_Std.mat'])).timelocked_data;
        squidgrad_timelocked = load(fullfile(save_path, [params.sub '_squidgrad_timelocked_Std.mat'])).timelocked_data;
        
%         % Prepare
%         prep_leadfield(headmodels,sourcemodel, squidmag_timelocked, [params.sub '_squidmag'], save_path)
%         prep_leadfield(headmodels,sourcemodel, squidgrad_timelocked, [params.sub '_squidgrad'], save_path)
%         prep_leadfield(headmodels,sourcemodel, opm_timelockedT, [params.sub '_opm'], save_path)
    end


%% MNE fit
    if exist(fullfile(save_path, 'mne_fits.mat'),'file') && overwrite.mne==false
        disp(['Not overwriting MNE source reconstruction for ' params.sub]);
    else
        %LOAD LEADFIELDS AND SOURCEMODELS
        clear headmodels sourcemodel sourcemodel_inflated
        sourcemodel = load(fullfile(shared_folder_path, [params.sub '_sourcemodel.mat'])).sourcemodel;
        sourcemodel_inflated = load(fullfile(shared_folder_path, [params.sub '_sourcemodel_inflated.mat'])).sourcemodel_inflated;
        headmodels = load(fullfile(shared_folder_path, 'headmodels.mat')).headmodels;
        
        clear squimdag_timelocked squidgrad_timelocked opm_timelockedT
        opm_timelockedT = load(fullfile(save_path, 'source analysis', [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
        squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked_Std.mat'])).timelocked_data;
        squidgrad_timelocked = load(fullfile(save_path, [params.sub '_squidgrad_timelocked_Std.mat'])).timelocked_data;

        params.use_cov_all = false;
        
        leadfield = load(fullfile(save_path, 'source analysis', [params.sub '_squidmag_leadfield.mat'])).leadfield;
        params.modality = 'squidmag';
        row_idx = find(contains(stats.lab, 'M300, No Go')); % Row index for latency
        val = stats.val_squidmag(row_idx,i_sub);
        fit_mne(save_path, squidmag_timelocked, headmodels, sourcemodel, sourcemodel_inflated, leadfield, params, val);
        
        leadfield = load(fullfile(save_path, 'source analysis', [params.sub '_squidgrad_leadfield.mat'])).leadfield;
        params.modality = 'squidgrad';
        row_idx = find(contains(stats.lab, 'M300 No Go')); % Row index for latency
        val = stats.val_squidgrad(row_idx,i_sub);
        fit_mne(save_path, squidgrad_timelocked, headmodels, sourcemodel, sourcemodel_inflated, leadfield, params, val);
               
        leadfield = load(fullfile(save_path, 'source analysis', [params.sub '_opm_leadfield.mat'])).leadfield;
        params.modality = 'opm';
        row_idx = find(contains(stats.lab, 'M300 No Go')); % Row index for latency
        val = stats.val_opm(row_idx,i_sub);
        fit_mne(save_path, opm_timelockedT, headmodels, sourcemodel, sourcemodel_inflated, leadfield, params, val);

        clear leadfield
    end
end
close all

%% MNE group analysis (använd som inspo, ändra)
subs = [2:10 12:13];
mne_results_goup(base_save_path, subs, params);
