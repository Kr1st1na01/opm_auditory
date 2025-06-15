function [timelocked, peak, ylimit] = evoked_analysis(data, params, save_path, ylimit)
% Here we timelock the Std, high and low triggers and find the sensor with
% the highest peak for each one.
peak = [];
peak.values = [];
peak.labels = {};
peak.number_count = 1;

% Finding peaks during this time intervall
params.pretimwin = 0.08;
params.posttimwin = 0.12;

%% Std tone
params.trials = find(data.trialinfo==params.trigger_code(1));
params.condition = 'Std';
timelocked = timelock(data, params, save_path);
ylimit = plot_butterfly(timelocked, params, save_path, ylimit);

% Sensor with highest peak at 100 ms
params.condition = 'M100_Std';
[tmp, peak] = find_peak(timelocked, params, peak); % Saves values
timelocked.avg = timelocked.avg(tmp.i_peakch,:);
plot_butterfly(timelocked, params, save_path, ylimit);

params.pretimwin = 0.27;
params.posttimwin = 0.33;

% Sensor with highest peak at 300 ms
params.condition = 'M300_Std';
[tmp, peak] = find_peak(timelocked, params, peak); % Saves values
timelocked.avg = timelocked.avg(tmp.i_peakch,:);
plot_butterfly(timelocked, params, save_path, ylimit);

%% Low tone
params.pretimwin = 0.08;
params.posttimwin = 0.12;

params.trials = find(ismember(data.trialinfo, params.trigger_code(2:3)));
params.condition = 'Low';
timelocked = timelock(data, params, save_path);
ylimit = plot_butterfly(timelocked, params, save_path, ylimit);

% Sensor with highest peak at 100 ms
params.condition = 'M100_Low';
[tmp, peak] = find_peak(timelocked, params, peak); % Saves values
timelocked.avg = timelocked.avg(tmp.i_peakch,:);
plot_butterfly(timelocked, params, save_path, ylimit)

%% High tone
params.trials = find(ismember(data.trialinfo, params.trigger_code(4:5)));
params.condition = 'High';
timelocked = timelock(data, params, save_path);
ylimit = plot_butterfly(timelocked, params, save_path, ylimit);

% Sensor with highest peak at 100 ms
params.condition = 'M100_High';
[tmp, peak] = find_peak(timelocked, params, peak);
timelocked.avg = timelocked.avg(tmp.i_peakch,:);
plot_butterfly(timelocked, params, save_path, ylimit)

%% Go trigger
G = [];
G = find(ismember(data.trialinfo, [params.trigger_code(3) params.trigger_code(5)])); % all Go triggers saved in one place

params.trials = G;
params.condition = 'Go';
timelocked_G = timelock(data, params, save_path);
ylimit = plot_butterfly(timelocked_G, params, save_path, ylimit);

timelocked_G100 = timelocked_G;
timelocked_G300 = timelocked_G;

% Sensor with highest peak at 100 ms
params.pretimwin = 0.08;
params.posttimwin = 0.12;
params.condition = 'M100_Go';
[tmp, peak] = find_peak(timelocked_G, params, peak);

timelocked_G100.avg = timelocked_G.avg(tmp.i_peakch,:);
plot_butterfly(timelocked_G100, params, save_path, ylimit);

% Sensor with highest peak at 300 ms
params.pretimwin = 0.27;
params.posttimwin = 0.33;
params.condition = 'M300_Go';
[tmp, peak] = find_peak(timelocked_G, params, peak);

timelocked_G300.avg = timelocked_G.avg(tmp.i_peakch,:);
plot_butterfly(timelocked_G300, params, save_path, ylimit);

%% Saving trigger indeces
pre_NG = [];
NG = find(ismember(data.trialinfo, [params.trigger_code(2) params.trigger_code(4)])); % all No Go triggers saved in one place

for i = 1:size(NG)
    pre_NG = [pre_NG; NG(i)-1]; % for every index we save the one before, so the standard tone is saved
end

%% Timelocking data
% pre No Go
params.trials = pre_NG;
params.condition = 'pre No Go';
timelocked_pre = timelock(data, params, save_path);
plot_butterfly(timelocked_pre, params, save_path, ylimit);

% No Go
params.trials = NG;
params.condition = 'No Go';
timelocked_NG = timelock(data, params, save_path);
ylimit = plot_butterfly(timelocked_NG, params, save_path, ylimit);

%% MMN
params.trials = NG;
params.condition = 'MMN';
timelockedMMN = timelocked_NG;
timelockedMMN.avg = timelocked_NG.avg - timelocked_pre.avg; % The averaged timelocked data of the NG and standard tones are subtracted to get the MMN
plot_butterfly(timelockedMMN, params, save_path, ylimit)

%% Plots
% pre No Go at 100 ms
params.pretimwin = 0.08;
params.posttimwin = 0.12;
params.condition = 'M100_pre No Go';
[tmp, peak] = find_peak(timelocked_pre, params, peak);
timelocked_pre.avg = timelocked_pre.avg(tmp.i_peakch,:);
plot_butterfly(timelocked_pre, params, save_path, ylimit)

timelocked_NG100 = timelocked_NG;
timelocked_NG300 = timelocked_NG;

% No Go at 100 ms
params.condition = 'M100_No Go';
[tmp, peak] = find_peak(timelocked_NG100, params, peak);
timelocked_NG100.avg = timelocked_NG100.avg(tmp.i_peakch,:);
plot_butterfly(timelocked_NG100, params, save_path, ylimit);

% No Go at 300 ms
params.pretimwin = 0.27;
params.posttimwin = 0.33;
params.condition = 'M300_No Go';
[tmp, peak] = find_peak(timelocked_NG300, params, peak);
timelocked_NG300.avg = timelocked_NG300.avg(tmp.i_peakch,:);
plot_butterfly(timelocked_NG300, params, save_path, ylimit);

%% MMN
params.pretimwin = 0.1; % for find_peak.m
params.posttimwin = 0.25;
params.condition = 'MMN peak';

[tmp, peak] = find_peak(timelockedMMN, params, peak);

% Plot butterfly MMN with xline at the found peak
h = figure;
plot(timelockedMMN.time*1e3,timelockedMMN.avg*params.amp_scaler);
hold on
xline(timelockedMMN.time(tmp.i_peaktime)*1000, '--') % Choose timepoint, x-axis
hold off
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelockedMMN.cfg.trials)) ')'], Interpreter="none")
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Evoked data_' params.condition '.jpg']))

% Plot max sensor with xline
h = figure;
plot(timelockedMMN.time*1e3,timelockedMMN.avg*params.amp_scaler);
hold on
xline(timelockedMMN.time(tmp.i_peaktime)*1000, '--') % Choose timepoint, x-axis
hold off
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelockedMMN.cfg.trials)) ')'], Interpreter="none")
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Evoked data_' params.condition '.jpg']))

close all