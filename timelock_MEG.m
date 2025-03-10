function [timelocked] = timelock_MEG(MMN_data, TFR_data, params, save_path)
% Here we timelock the Std, high and low triggers and find the sensor with
% the highest peak for each one.

% % Det borde vara rätt channel 
cfg = [];
cfg.channel = params.chs;
MMN_data = ft_selectdata(cfg, MMN_data);
TFR_data = ft_selectdata(cfg, TFR_data);

%% Timelock of Std
params.trials = find(MMN_data.trialinfo==params.trigger_code(1));
params.condition = 'Std';
timelocked = timelock(MMN_data, params, save_path);
plot_butterfly(timelocked, params, save_path)

% Sensor with highest peak
[~, pks_i] = max(max(abs(timelocked.avg), [], 2)); % 2 indicates the max value of rows, max returns the column vector containing the max value of each row, therefore a second max
timelocked.avg = timelocked.avg(pks_i,:);
params.condition = 'Sensor with highest peak for Std trigger';
params.pk_val = freqanalysis(TFR_data, params, save_path);
plot_butterfly(timelocked, params, save_path)


%% Timelocked of Low
% Plot both low triggers
params.trials = find(ismember(MMN_data.trialinfo, params.trigger_code(2:3)));
params.condition = 'Low';
timelocked = timelock(MMN_data, params, save_path);
plot_butterfly(timelocked, params, save_path)

% Sensor with highest peak
[~, pks_i] = max(max(abs(timelocked.avg), [], 2));
timelocked.avg = timelocked.avg(pks_i,:);
params.condition = 'Sensor with highest peak for low trigger';
params.pk_val = freqanalysis(TFR_data, params, save_path);
plot_butterfly(timelocked, params, save_path)

%% Timelocked of High
% Plot both high triggers
params.trials = find(ismember(MMN_data.trialinfo, params.trigger_code(4:5)));
params.condition = 'High';
timelocked = timelock(MMN_data, params, save_path);
plot_butterfly(timelocked, params, save_path)

% Sensor with highest peak
[~, pks_i] = max(max(abs(timelocked.avg), [], 2));
timelocked.avg = timelocked.avg(pks_i,:);
params.condition = 'Sensor with highest peak for high trigger';
params.pk_val = freqanalysis(TFR_data, params, save_path);
plot_butterfly(timelocked, params, save_path)

%% NG and G for TFR
NG = [];
G = [];
NGL = find(TFR_data.trialinfo==params.trigger_code(2));
NGH = find(TFR_data.trialinfo==params.trigger_code(4)); % Finds index of trigger 3

GL = find(TFR_data.trialinfo==params.trigger_code(3));
GH = find(TFR_data.trialinfo==params.trigger_code(5)); % Finds index of trigger 3

if size(NGL) > size(NGH)
    for i = 1:length(NGL) % No Go index
        NG = [NG; NGL(i)]; % Adds index of trigger 3
        G = [G; GH(i)]; 
    end
else
    for i = 1:length(NGH) % No Go index
        NG = [NGH; NGH(i)]; % Adds that preceding index
        G = [G; GL(i)]; 
    end
end

params.trials = NG;
params.condition = 'No Go trigger';
timelocked_NG = timelock(TFR_data, params, save_path);
params.pk_val = freqanalysis(TFR_data, params, save_path);
plot_butterfly(timelocked, params, save_path)

params.trials = G;
params.condition = 'Go trigger';
timelocked_G = timelock(TFR_data, params, save_path);
params.pk_val = freqanalysis(TFR_data, params, save_path);
plot_butterfly(timelocked, params, save_path)

close all