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
cfg.foi         = 30:1:50;             % Frequencies we want to estimate, total 21 frequencies
cfg.toi         = -params.pre:0.05:params.post; % Times to center on, the windows are 500 ms, from -0.3-0.9=1.2s
cfg.t_ftimwin   = ones(length(cfg.foi), 1)*0.5;         % length of time window, 0.5s.
cfg.tapsmofrq   = zeros(length(cfg.foi), 1)+2;        % Smoothing. dF = 1/dT (dT = 0.5 --> dF = 2 --> smoothing needs to be 2 Hz).

TFRhann_multi = ft_freqanalysis(cfg, epochs); % The power spectrum is calculated

% Plot TFR
cfg.parameter       = 'powspctrm';
cfg.layout          = params.layout;
cfg.showlabels      = 'yes';
cfg.baselinetype    = 'relative';
cfg.baseline        = [-params.pre 0];

h = figure;
ft_multiplotTFR(cfg, TFRhann_multi);
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_TFR multi_' params.condition '.jpg']))

% Topoplot
h = figure;
cfg.channel = pks;
ft_topoplotTFR(cfg, TFRhann_multi)
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality 'Topoplot TFR multi_' params.condition '.jpg']))

% Channel with highest peak
h = figure;
cfg.channel = pks;
ft_singleplotTFR(cfg, TFRhann_multi)
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_single TFR multi_' params.condition '.jpg']))

%% TFR: Wavelet analysis
% The "tapering" is done on a wave-function that is fitted to the signal
% centred on the timepoints of interest.
cfg = [];
cfg.channel     = 'all';
cfg.method      = 'wavelet';
cfg.foi         = 30:1:50;
cfg.toi         = -params.pre:0.01:params.post;
cfg.width       = 5;                        % Number of cycles
cfg.pad         = 'nextpow2';

TFR_wavelet = ft_freqanalysis(cfg, epochs);

% Plot TFR
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = params.layout;
cfg.showlabels      = 'yes';
cfg.baselinetype    = 'relative';
cfg.baseline        = [-params.pre 0];

h = figure;
ft_multiplotTFR(cfg, TFR_wavelet);
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_TFR wave_' params.condition '.jpg']))

% Channel with highest peak
h = figure;
cfg.channel = pks;
ft_singleplotTFR(cfg, TFRhann_multi)
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_single TFR wave_' params.condition '.jpg']))

end