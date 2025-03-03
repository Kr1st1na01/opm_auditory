function [] = freqanalysis(data, params, pks, save_path)
% Här gör vi TFR

cfg = [];
cfg.trials = params.trials;
epochs = ft_selectdata(cfg, data);

%% TFR: multitaper
cfg = [];
cfg.output      = 'pow';
cfg.channel     = 'all';
cfg.method      = 'mtmconvol';
cfg.taper       = 'dpss';             % Slepian sequence as tapers
cfg.foi         = 1:1:45;             % Frequencies we want to estimate from 1 Hz to 45 Hz in steps of 1HZ
cfg.toi         = -params.pre:0.01:params.post; % Times to center on
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

% Channel with highest peak
h = figure;
cfg.channel = pks;
ft_singleplotTFR(cfg, tfr_dpss)
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_single TFR multi_' params.condition '.jpg']))

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

% Channel with highest peak
h = figure;
cfg.channel = pks;
ft_singleplotTFR(cfg, tfr_dpss)
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_single TFR wave_' params.condition '.jpg']))

end