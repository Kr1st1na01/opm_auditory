function [peak] = MMN(MMN_data, TFR_data, params, save_path, peak) % data = MMN_data for modality
% Här gör vi MMN

cfg = [];
cfg.channel = params.chs;
MMN_data = ft_selectdata(cfg, MMN_data);
TFR_data = ft_selectdata(cfg, TFR_data);

%% NG trigger
G = [];

GL = find(MMN_data.trialinfo==params.trigger_code(3));
GH = find(MMN_data.trialinfo==params.trigger_code(5));

G = [G; GL];
G = [G; GH];

params.trials = G;
params.condition = 'Go trigger';
timelocked_G = timelock(MMN_data, params, save_path);
params.freqylim = [41 45];
freqanalysis(TFR_data, params, save_path, peak);
plot_butterfly(timelocked_G, params, save_path)

% Max channel
params.pretimwin = 0.08;
params.posttimwin = 0.12;
[peak, tmp] = findpeaks(timelocked_G, params, peak);
timelocked_G.avg = timelocked_G.avg(tmp.i_peaktime,:);
peak.values(end+1,1) = tmp.amplitude;
peak.labels(end+1,1) = {[params.modality '_MMN data_M100_' params.condition]};
plot_butterfly(timelocked_G, params, save_path)

%% Saving trigger indeces
% High NG
% Low NG
NG = [];
pre_NG = [];
H_NG = find(MMN_data.trialinfo==params.trigger_code(4)); % Finds index of trigger 11, High No Go
L_NG = find(MMN_data.trialinfo==params.trigger_code(2)); % Finds index of trigger 3, Low No Go

for i = 1:length(H_NG)
    NG = [Oddball; H_NG(i)]; % Adds index
    pre_NG = [pre_NG; H_NG(i)-1]; % Adds that preceding index
end
for i = 1:length(L_NG)
    NG = [Oddball; LowNG(i)];
    pre_NG = [pre_NG; LowNG(i)-1];
end

%% Timelocking data
% pre No Go
params.trials = pre_NG;
params.condition = 'pre No Go trigger';
timelocked_pre = timelock(MMN_data, params, save_path);
params.freqylim = [41 45];
freqanalysis(TFR_data, params, save_path, peak);
plot_butterfly(timelocked_pre, params, save_path)

% No Go
params.trials = NG;
params.condition = 'No Go trigger';
timelocked = timelock(MMN_data, params, save_path);
params.freqylim = [41 45];
freqanalysis(TFR_data, params, save_path, peak);
plot_butterfly(timelocked, params, save_path)

%% MMN
params.trials = Oddball;
params.condition = 'pre NG vs NG';
timelockedMMN = timelocked;
timelockedMMN.avg = timelocked.avg - timelocked_pre.avg; % Here the oddball timelocked data average is changed to represent the difference (MMN) in response
plot_butterfly(timelockedMMN, params, save_path)

%% M100 plots
% pre No Go
params.pretimwin = 0.08;
params.posttimwin = 0.12;

params.condition = 'pre No Go trigger';
[peak, tmp] = findpeaks(timelocked_pre, params, peak);
timelocked_pre.avg = timelocked_pre.avg(tmp.i_peaktime,:);
peak.values(end+1,1) = tmp.amplitude;
peak.labels(end+1,1) = {[params.modality '_MMN data_M100_' params.condition]};
plot_butterfly(timelocked_pre, params, save_path)

% No Go
params.condition = 'No Go trigger';
[peak, tmp] = findpeaks(timelocked, params, peak);
timelocked.avg = timelocked.avg(tmp.i_peaktime,:);
peak.values(end+1,1) = tmp.amplitude;
peak.labels(end+1,1) = {[params.modality '_MMN data_M100_' params.condition]};
plot_butterfly(timelocked, params, save_path)

%% MMN and M300
params.pretimwin = 0.27;
params.posttimwin = 0.33;
[peak, tmp] = findpeaks(timelockedMMN, params, peak);
params.condition = 'pre NG vs NG';

% Plot MMN with M300
h = figure;
plot(timelockedMMN.time*1e3,timelockedMMN.avg*params.amp_scaler);
hold on
xline(timelockedMMN.time(tmp.i_peaktime)*1000, '--') % Choose timepoint, x-axis
hold off
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelockedMMN.cfg.trials)) ')'])
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_MMN data_M300_' params.condition '.jpg']))

% Plotta max ch
timelockedMMN.avg = timelockedMMN.avg(tmp.i_peakch, :); % Choose peak ch, y-axis
params.condition = 'MMN';
plot_butterfly(timelockedMMN, params, save_path)
save(fullfile(save_path, [params.sub '_' params.modality '_MMN data_M300_' params.condition]), 'timelocked', '-v7.3'); 

h = figure;
plot(timelockedMMN.time*1e3,timelockedMMN.avg*params.amp_scaler);
hold on
xline(timelockedMMN.time(tmp.i_peaktime)*1000, '--') % Choose timepoint, x-axis
hold off
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelockedMMN.cfg.trials)) ')'])
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_MMN data_max ch_' params.condition '.jpg']))



