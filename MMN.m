function [] = MMN(data, TFR, params, save_path) % data = MMN_data for modality
% Här gör vi MMN

cfg = [];
cfg.channel = params.chs;
data = ft_selectdata(cfg, data);

%% Saving trigger indeces
% High NG
% Low NG
Oddball = [];
preOddball = [];
High = find(data.trialinfo==params.trigger_code(4));
Low = find(data.trialinfo==params.trigger_code(2)); % Finds index of trigger 3

for i = 1:length(High)
    Oddball = [Oddball; High(i)]; % Adds index of trigger 3
    preOddball = [preOddball; High(i)-1]; % Adds that preceding index
end
for i = 1:length(Low)
    Oddball = [Oddball; Low(i)]; % Adds index of trigger 3
    preOddball = [preOddball; Low(i)-1]; % Adds that preceding index
end

%% Timelocking data
params.trials = preOddball;
params.condition = 'Pre-oddball trigger for No Go';
timelocked_pre = timelock(data, params, save_path);
plot_butterfly(timelocked_pre, params, save_path)
[~, pks_i] = max(max(abs(timelocked_pre.avg), [], 2)); % 2 indicates the max value of rows, max returns the column vector containing the max value of each row, therefore a second max
%freqanalysis(TFR, params, pks_i, save_path)

params.trials = Oddball;
params.condition = 'Oddball trigger for No go';
timelocked = timelock(data, params, save_path);
plot_butterfly(timelocked, params, save_path)
[~, pks_i] = max(max(abs(timelocked_pre.avg), [], 2)); % 2 indicates the max value of rows, max returns the column vector containing the max value of each row, therefore a second max
%freqanalysis(TFR, params, pks_i, save_path)

%% MMN
% (Low No go + High No go) - (pre-Low No go + pre-High No go)
params.trials = Oddball;
params.condition = 'Trigger for No go vs pre-No go';
timelocked.avg = timelocked.avg - timelocked_pre.avg; % Here the oddball timelocked data average is changed to represent the difference (MMN) in response
plot_butterfly(timelocked, params, save_path)

%%
% Här vill jag hitta peaks för MMN. Plotta butterfly med linje om jag vill
% och sen peak sensor med linje
% Plotting

% Timelocked for MMN
tmp = [];
% Här får vi index av tiderna där det går från -0.1 s till 0.5 s
[~, interval_P300(1)] = min(abs(timelocked.time-0.027)); % find closest time sample to 270 ms
[~, interval_P300(2)] = min(abs(timelocked.time-0.033)); % find closest time sample to 330 ms

% Sensor med max peak i intervallet
[~, pks_i] = max(max(abs(timelocked.avg(:, interval_P300(1):interval_P300(2))), [], 2)); % Jag vill ha raden med max

% PLotta den
timelocked.avg_org = timelocked.avg;
timelocked.avg = timelocked.avg(pks_i, :);
params.condition = 'Sensor with the highest peak for MMN';
plot_butterfly(timelocked, params, save_path)
save(fullfile(save_path, [params.sub '_' params.modality '_' params.condition]), 'timelocked', '-v7.3'); 

% Hitta tiden för peak vid intervallet för den sensorn
[~, loc] = findpeaks(abs(timelocked.avg(:, interval_P300(1):interval_P300(2))),'SortStr','descend');

% Största värdet i datan, ger index för största värdet
if isempty(loc)
    loc = round((interval_P300(2)-interval_P300(1))/2);
    tmp.nopeak = true;
end

% % Peak sensor plot
% %[~, pks_i] = max(max(abs(timelocked.avg), [], 2)); % 2 indicates the max value of rows, max returns the column vector containing the max value of each row, therefore a second max
% [~, I] = max(timelocked.avg(:, loc+interval_P300(1)));
% timelocked.avg = timelocked.avg(I, :);
% params.condition = 'Sensor with the highest peak for MMN';
% plot_butterfly(timelocked, params, save_path)
% save(fullfile(save_path, [params.sub '_' params.modality '_' params.condition]), 'timelocked', '-v7.3'); 

% Plot
timelocked.avg = timelocked.avg_org;
h = figure;
plot(timelocked.time*1e3,timelocked.avg*params.amp_scaler);
hold on
xline(timelocked.time(loc+interval_P300(1))*10000)
hold off
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelocked.cfg.trials)) ')'])
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_MMN-' params.condition '.jpg']))




