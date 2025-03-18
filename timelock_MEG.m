function [timelocked, peak] = timelock_MEG(MMN_data, TFR_data, params, save_path, peak)
% Here we timelock the Std, high and low triggers and find the sensor with
% the highest peak for each one.

% % Det borde vara rätt channel 
cfg = [];
cfg.channel = params.chs;
MMN_data = ft_selectdata(cfg, MMN_data);
TFR_data = ft_selectdata(cfg, TFR_data);

params.pretimwin = 0.08;
params.posttimwin = 0.12;

%% Timelock of Std
params.trials = find(MMN_data.trialinfo==params.trigger_code(1));
params.condition = 'Std';
timelocked = timelock(MMN_data, params, save_path);
plot_butterfly(timelocked, params, save_path)

% Sensor with highest peak
[peak, tmp] = findpeaks(timelocked, params, peak);

timelocked.avg = timelocked.avg(tmp.i_peakch,:);
peak.labels(end+1,1) = {[params.modality '_MMN data_M100_' params.condition]};
peak.values(end+1,1) = tmp.amplitude; % A(4:5,5:6) = [2 3; 4 5], this expands the matrix
params.freqylim = [37 41];
freqanalysis(TFR_data, params, save_path, peak);
plot_butterfly(timelocked, params, save_path)


%% Timelocked of Low
% Plot both low triggers
params.trials = find(ismember(MMN_data.trialinfo, params.trigger_code(2:3)));
params.condition = 'Low';
timelocked = timelock(MMN_data, params, save_path);
plot_butterfly(timelocked, params, save_path)

% Sensor with highest peak
[peak, tmp] = findpeaks(timelocked, params, peak);

timelocked.avg = timelocked.avg(tmp.i_peakch,:);
peak.labels(end+1,1) = {[params.modality '_MMN data_M100_' params.condition]};
peak.values(end+1,1) = tmp.amplitude;
params.freqylim = [41 45];
freqanalysis(TFR_data, params, save_path, peak);
plot_butterfly(timelocked, params, save_path)

%% Timelocked of High
% Plot both high triggers
params.trials = find(ismember(MMN_data.trialinfo, params.trigger_code(4:5)));
params.condition = 'High';
timelocked = timelock(MMN_data, params, save_path);
plot_butterfly(timelocked, params, save_path)

% Sensor with highest peak
[peak, tmp] = findpeaks(timelocked, params, peak);
timelocked.avg = timelocked.avg(tmp.i_peaktime,:);
peak.labels(end+1,1) = {[params.modality '_MMN data_M100_' params.condition]};
peak.values(end+1,1) = tmp.amplitude;
params.freqylim = [41 45];
freqanalysis(TFR_data, params, save_path, peak);
plot_butterfly(timelocked, params, save_path)

peak = MMN(MMN_data, TFR_data, params, save_path, peak); % The MMN is done on cropped data and TFR is for the frequency analysis

close all