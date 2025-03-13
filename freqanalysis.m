function [] = freqanalysis(data, params, save_path, peak)
% Här gör vi TFR
cfg = [];
cfg.trials = params.trials;
epochs = ft_selectdata(cfg, data);

%% FFT
cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'prestim';
timelock = ft_timelockanalysis(cfg, epochs); % labelxtime

cfg = [];
cfg.output      = 'pow';
cfg.channel     = 'all';
cfg.method      = 'mtmfft';
cfg.taper       = 'hanning';             % Slepian sequence as tapers
cfg.foi         = 30:1:50;             % Frequencies we want to estimate, total 21 frequencies
cfg.pad = 2;
FFT_timelocked = ft_freqanalysis(cfg, timelock); % LabelxFreq

% Plot TFR
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = params.layout;
cfg.showlabels      = 'yes';

h = figure;
ft_multiplotER(cfg, FFT_timelocked);
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_FFT multi_' params.condition '.jpg']))

% Channel with highest peak
h = figure;
FFT_mean = mean(FFT_timelocked.powspctrm, 2);
[val, pks] = max(FFT_mean);
peak.values(end+1, 1) = val;
peak.labels = {peak.labels;[params.modality '_' params.condition '_FFT_multi']};
cfg.channel = pks;
ft_singleplotER(cfg, FFT_timelocked)
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_FFT multi_' params.condition '.jpg']))


% Topoplot
h = figure;
cfg = [];
cfg.baselinetype    = 'relative';
cfg.baseline        = [-params.pre 0];
cfg.layout = params.layout;
cfg.xlim = [params.pretimwin params.posttimwin]; %FIX TIME 0.4-0.6 ish
ft_topoplotTFR(cfg, FFT_timelocked)
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Topoplot FFT multi_' params.condition '.jpg']))


%% TFR: multitaper
cfg = [];
cfg.output      = 'pow';
cfg.channel     = 'all';
cfg.method      = 'mtmconvol';
cfg.taper       = 'dpss';             % Slepian sequence as tapers
cfg.foi         = 30:1:50;             % Frequencies we want to estimate, total 21 frequencies
cfg.toi         = -params.pre:0.05:params.post; % Times to center on, the windows are 500 ms, from -0.3-0.9=1.2s
cfg.t_ftimwin   = ones(length(cfg.foi), 1)*0.25;         % length of time window, 0.5s.
cfg.tapsmofrq   = zeros(length(cfg.foi), 1)+4;   % Smoothing. dF = 1/dT (dT = 0.5 --> dF = 2 --> smoothing needs to be 2 Hz).
cfg.pad = 2;

TFRhann_multi = ft_freqanalysis(cfg, epochs); % The power spectrum is calculated

% Plot TFR
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = params.layout;
cfg.showlabels      = 'yes';
cfg.baselinetype    = 'relative';
cfg.baseline        = [-params.pre 0];

h = figure;
ft_multiplotTFR(cfg, TFRhann_multi);
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_TFR multi_' params.condition '.jpg']))

% Channel with highest peak
h = figure;
val = nanmean(TFRhann_multi.powspctrm(pks, 10, :));
peak.values(end+1, 1) = val;
peak.labels = {peak.labels;[params.modality '_' params.condition '_TFR_multi']};
cfg.channel = pks;
ft_singleplotTFR(cfg, TFRhann_multi)
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_TFR multi_' params.condition '.jpg']))

% Topoplot
h = figure;
cfg = [];
cfg.baselinetype    = 'relative';
cfg.baseline        = [-params.pre 0];
cfg.layout = params.layout;
cfg.xlim = [params.pretimwin params.posttimwin]; %FIX TIME 0.4-0.6 ish
ft_topoplotTFR(cfg, TFRhann_multi)
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Topoplot TFR multi_' params.condition '.jpg']))

end