function [] = freqanalysis(data, params, save_path)
cfg = [];
cfg.trials = params.trials;
epochs = ft_selectdata(cfg, data);

% cfg = [];
% cfg.output          = 'pow';            % Return PSD
% cfg.channel         = 'all';            % Calculate for MEG and EEG
% cfg.method          = 'mtmfft';
% cfg.pad             = 'nextpow2';       % Trial padding
% cfg.parameter       = 'powspctrm';
% cfg.showlabels      = 'yes';
% cfg.layout          = params.layout; % Layout for MEG/EEG sensors
% cfg.foilim          = [1 95];           % Frequency range
% cfg.xlim            = [1 45];           % Frequencies to plot
% 
% %% PSD: Singletaper
% cfg.taper           = 'hanning';        % Hann window as taper, lets us define the windowing/tapering function
% 
% psd_hann = ft_freqanalysis(cfg, epochs);
% 
% % Plot PSD
% h = figure;
% ft_multiplotER(cfg, psd_hann);
% saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_single taper_' params.condition '.jpg']))
% 
% %% PSD: Multitaper
% cfg.taper           = 'dpss';         % Multitapers based on Slepian sequences
% cfg.tapsmofrq       = 2;              % Smoothing +/- 2 Hz, change smoothing and try maybe 10
% 
% psd_dpss = ft_freqanalysis(cfg, epochs);
% 
% % Plot PSD
% h = figure;
% ft_multiplotER(cfg, psd_dpss);
% saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_multi taper_' params.condition '.jpg']))
% 
% %% TFR: Singletaper
% cfg = [];
% cfg.output      = 'pow';      
% cfg.channel     = params.chs;
% cfg.method      = 'mtmconvol';
% cfg.taper       = 'hanning';    % Hann window as taper
% cfg.foi         = 1:1:45;       % Frequencies we want to estimate from 1 Hz to 45 Hz in steps of 1HZ
% cfg.toi         = -0.5:0.01:0.9;            % Timepoints to center on
% cfg.t_ftimwin   = 0.5*ones(length(cfg.foi),1);  % length of time window
% 
% tfr_hann = ft_freqanalysis(cfg, epochs);
% 
% % Plot TFR
% cfg = [];
% cfg.parameter       = 'powspctrm';
% cfg.layout          = params.layout;
% cfg.showlabels      = 'yes';
% 
% h = figure;
% ft_multiplotTFR(cfg, tfr_hann);
% saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_TFR single_' params.condition '.jpg']))
% 
% % Plot relative TFR with baseline correction
% cfg.baselinetype    = 'relative';  % Type of baseline, see help ft_multiplotTFR
% cfg.baseline        = [0 inf];    % Time of baseline. Here the entire epoch is used as baseline.
% 
% h = figure;
% ft_multiplotTFR(cfg, tfr_hann);
% saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_TFR baseline single_' params.condition '.jpg']))
% 
% %% TFR: Single taper, varying taper length
% % The time windows here vary with the individual frequencies of interest.
% cfg = [];
% cfg.output      = 'pow';
% cfg.channel     = params.chs;
% cfg.method      = 'mtmconvol';
% cfg.taper       = 'hanning';         % Hann window as taper
% cfg.foi         = 1:1:45;            % Frequencies we want to estimate from 1 Hz to 45 Hz in steps of 1HZ
% cfg.toi         = -0.5:0.01:0.9; % Times to center on
% cfg.t_ftimwin    = 5./cfg.foi; 
% 
% tfr_hann5 = ft_freqanalysis(cfg, epochs);
% 
% % Plot TFR
% cfg.parameter       = 'powspctrm';
% cfg.layout          = params.layout;
% cfg.showlabels      = 'yes';
% cfg.baselinetype    = 'relative';
% cfg.baseline        = [0 inf];
% 
% h = figure;
% ft_multiplotTFR(cfg, tfr_hann5);
% saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_TFR single5_' params.condition '.jpg']))

%% TFR: multitaper
cfg = [];
cfg.output      = 'pow';
cfg.channel     = 'all';
cfg.method      = 'mtmconvol';
cfg.taper       = 'dpss';             % Slepian sequence as tapers
cfg.foi         = 1:1:45;             % Frequencies we want to estimate from 1 Hz to 45 Hz in steps of 1HZ
cfg.toi         = -0:0.01:0.9;  % Times to center on
cfg.t_ftimwin   = 5./cfg.foi;         % length of time window
cfg.tapsmofrq   = 0.5*cfg.foi;        % Smoothing

tfr_dpss = ft_freqanalysis(cfg, epochs);

% Plot TFR
cfg.parameter       = 'powspctrm';
cfg.layout          = params.layout;
cfg.showlabels      = 'yes';
cfg.baselinetype    = 'relative';
cfg.baseline        = [0 inf];

h = figure;
ft_multiplotTFR(cfg, tfr_dpss);
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_TFR multi_' params.condition '.jpg']))

%% TFR: Wavelet analysis
% The "tapering" is done on a wave-function that is fitted to the signal
% centred on the timepoints of interest.
cfg = [];
cfg.channel     = 'all';
cfg.method      = 'wavelet';
cfg.foi         = 1:1:45;
cfg.toi         = -params.pre:0.01:params.post;
cfg.width       = 5;                        % Number of cycles
cfg.pad         = 'nextpow2';

tfr_wavelet = ft_freqanalysis(cfg, epochs);

% Plot TFR
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = params.layout;
cfg.showlabels      = 'yes';
cfg.baselinetype    = 'relative';
cfg.baseline        = [-0.1 0];

h = figure;
ft_multiplotTFR(cfg, tfr_wavelet);
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_TFR wave_' params.condition '.jpg']))


end