function [peak] = freqanalysis(data, params, save_path, peak)
% FFT and TFR is done here

cfg = [];
cfg.trials = params.trials;
epochs = ft_selectdata(cfg, data);

%% FFT
cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'prestim';
timelock = ft_timelockanalysis(cfg, epochs);

cfg = [];
cfg.output      = 'pow';
cfg.channel     = 'all';
cfg.method      = 'mtmfft';
cfg.taper       = 'hanning';           
cfg.foi         = 30:1:50;             % Frequencies we want to estimate, total 21 frequencies
cfg.pad = 2;
FFT_timelocked = ft_freqanalysis(cfg, timelock);

% Plot butterfly on the timelocked data
h = figure;
plot(FFT_timelocked.freq, FFT_timelocked.powspctrm);
xlabel('frequencies [Hz]')
ylabel(params.pow_label)
title(['Frequency analysis ' params.modality ' - ' params.condition ' (n_{chs}=' num2str(length(FFT_timelocked.cfg.channel)) ')'], Interpreter="none")
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Freq analysis_butterfly audodd_' params.condition '.jpg']))

% Plot topographical FFT
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = params.layout;
cfg.showlabels      = 'yes';

h = figure;
ft_multiplotER(cfg, FFT_timelocked);
title(['Evoked ' params.modality ' - FFT ' params.condition ' (n_{chs}=' num2str(length(FFT_timelocked.cfg.channel)) ')'], Interpreter="none")
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Freq tag_FFT_' params.condition '.jpg']))

% Channel with highest peak
h = figure;
[~, pks] = max(FFT_timelocked.powspctrm(:,find(params.specificfreq == FFT_timelocked.freq)));
peak.values(peak.number_count,1) = pks;
peak.labels(peak.number_count,1) = {[params.modality '_Freq tag_FFT_' params.condition '_freq ' int2str(params.specificfreq) '_' params.pow_label]};
cfg.channel = pks;
ft_singleplotER(cfg, FFT_timelocked)
xlabel('frequency [Hz]')
ylabel(params.pow_label)
title(['Evoked ' params.modality ' - FFT ' params.condition ' (n_{chs}=' num2str(length(FFT_timelocked.cfg.channel)) ')'], Interpreter="none")
legend(FFT_timelocked.label(pks), Interpreter="none")
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Freq tag_FFT_max sensor ' params.condition '.jpg']))

% Channel with second highest peak
FFT_timelocked_ori = FFT_timelocked;
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = params.layout;
cfg.showlabels      = 'yes';
FFT_timelocked.powspctrm(pks,find(params.specificfreq == FFT_timelocked.freq)) = 0;
[~, pks2] = max(FFT_timelocked.powspctrm(:,find(params.specificfreq == FFT_timelocked.freq)));
h = figure;
cfg.channel = pks2;
ft_singleplotER(cfg, FFT_timelocked_ori)
xlabel('frequency [Hz]')
ylabel(params.pow_label)
title(['Evoked ' params.modality ' - FFT ' params.condition ' (n_{chs}=' num2str(length(FFT_timelocked_ori.cfg.channel)) ')'], Interpreter="none")
legend(FFT_timelocked.label(pks2), Interpreter="none")
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Freq tag_FFT_max sensor ' params.condition ' 2.jpg']))

% Plot activation in topoplot
h = figure;
cfg = [];
cfg.baselinetype    = 'relative';
cfg.baseline        = [-params.pre 0];
cfg.layout = params.layout;
cfg.xlim = [params.pretimwin params.posttimwin];
ft_topoplotTFR(cfg, FFT_timelocked)
title(['Evoked ' params.modality ' - FFT ' params.condition ' (n_{chs}=' num2str(length(FFT_timelocked.cfg.channel)) ')'], Interpreter="none")
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Freq tag_FFT_Topoplot_' params.condition '.jpg']))


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
title(['Evoked ' params.modality ' - TFR ' params.condition ' (n_{chs}=' num2str(length(TFRhann_multi.cfg.channel)) ')'], Interpreter="none")
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Freq tag_TFR_' params.condition '.jpg']))

% Channel with highest peak
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = params.layout;
cfg.showlabels      = 'yes';
cfg.baselinetype    = 'relative';
cfg.baseline        = [-params.pre 0];

h = figure;
cfg.channel = pks;
ft_singleplotTFR(cfg, TFRhann_multi)
xlabel('time [s]')
ylabel('frequency [Hz]')
title(['Evoked ' params.modality ' - TFR ' params.condition],TFRhann_multi.label(pks), Interpreter="none")
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Freq tag_TFR_max sensor ' params.condition '.jpg']))

% % Channel with second highest peak
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = params.layout;
cfg.showlabels      = 'yes';
cfg.baselinetype    = 'relative';
cfg.baseline        = [-params.pre 0];

h = figure;
cfg.channel = pks2;
ft_singleplotTFR(cfg, TFRhann_multi)
xlabel('time [s]')
ylabel('frequency [Hz]')
title(['Evoked ' params.modality ' - TFR ' params.condition], TFRhann_multi.label(pks2), Interpreter="none")
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Freq tag_TFR_max sensor ' params.condition ' 2.jpg']))

% Topoplot showing which sensor was picked with a white star
h = figure;
cfg = [];
cfg.baselinetype    = 'relative';
cfg.baseline        = [-params.pre 0];
cfg.layout = params.layout;
cfg.highlight = 'on';
cfg.highlightchannel = pks;
cfg.highlightcolor = 'white';
cfg.highlightsymbol = 'pentagram';
cfg.highlightsize = 12;
cfg.xlim = [params.freqwin(1) params.freqwin(end)];
cfg.ylim = [params.freqylim];
ft_topoplotTFR(cfg, TFRhann_multi)
title(['Evoked ' params.modality ' - TFR ' params.condition],TFRhann_multi.label(pks), Interpreter="none")
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Freq tag_TFR_' params.condition ' with channel selection.jpg']))

% Second Topoplot showing which sensor was picked with a white star
h = figure;
cfg = [];
cfg.baselinetype    = 'relative';
cfg.baseline        = [-params.pre 0];
cfg.layout = params.layout;
cfg.highlight = 'on';
cfg.highlightchannel = pks2;
cfg.highlightcolor = 'white';
cfg.highlightsymbol = 'pentagram';
cfg.highlightsize = 12;
cfg.xlim = [params.freqwin(1) params.freqwin(end)];
cfg.ylim = [params.freqylim];
ft_topoplotTFR(cfg, TFRhann_multi)
title(['Evoked ' params.modality ' - TFR ' params.condition],TFRhann_multi.label(pks2), Interpreter="none")
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Freq tag_TFR_' params.condition ' with channel selection 2.jpg']))

% Topoplot done in small steps
s = (params.freqwin(end)-params.freqwin(1))/params.freqsteps;
fignum = s/s;
li = params.freqwin(1);
while li < params.freqwin(end)
    h = figure;
    cfg = [];
    cfg.baselinetype    = 'relative';
    cfg.baseline        = [-params.pre 0];
    cfg.layout = params.layout;
    cfg.xlim = [li li+params.freqsteps];
    cfg.ylim = [params.freqylim];
    ft_topoplotTFR(cfg, TFRhann_multi)
    title(['Evoked ' params.modality ' - TFR ' params.condition ' (n_{chs}=' num2str(length(TFRhann_multi.cfg.channel)) ')'], Interpreter="none")
    colorbar
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Freq tag_TFR_' params.condition '_' int2str(fignum) '.jpg']))
    li = li+params.freqsteps;
    fignum = fignum+1;
end
close all
end