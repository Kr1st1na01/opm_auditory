function [peak] = MMN(MMN_data, TFR_data, params, save_path, peak) % data = MMN_data for modality
% Här gör vi MMN

cfg = [];
cfg.channel = params.chs;
MMN_data = ft_selectdata(cfg, MMN_data);

%% Saving trigger indeces
% High NG
% Low NG
Oddball = [];
preOddball = [];
High = find(MMN_data.trialinfo==params.trigger_code(4));
Low = find(MMN_data.trialinfo==params.trigger_code(2)); % Finds index of trigger 3

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
timelocked_pre = timelock(MMN_data, params, save_path);
params.freqylim = [41 45];
freqanalysis(TFR_data, params, save_path, peak);
plot_butterfly(timelocked_pre, params, save_path)

params.trials = Oddball;
params.condition = 'Oddball trigger for No go';
timelocked = timelock(MMN_data, params, save_path);
params.freqylim = [41 45];
freqanalysis(TFR_data, params, save_path, peak);
plot_butterfly(timelocked, params, save_path)

%% MMN
% (Low No go + High No go) - (pre-Low No go + pre-High No go)
params.trials = Oddball;
params.condition = 'Trigger for No go vs pre-No go';
timelockedMMN = timelocked;
timelockedMMN.avg = timelocked.avg - timelocked_pre.avg; % Here the oddball timelocked data average is changed to represent the difference (MMN) in response
plot_butterfly(timelockedMMN, params, save_path)

%% MMN max ch butterfly plot
% Calculating the peak ch and value

tmp = [];

% The indexes of the times are saved
[~, interval_P300(1)] = min(abs(timelockedMMN.time-params.pretimwin)); % find closest time sample to 270 ms
[~, interval_P300(2)] = min(abs(timelockedMMN.time-params.posttimwin)); % find closest time sample to 330 ms

% Sensor med max peak i intervallet
[~, pks_time] = findpeaks(std(timelockedMMN.avg(:, interval_P300(1):interval_P300(2)), [], 1), 'SortStr','descend'); % Hitta tiden med bäst std. Använd tiden för att hitta relevanta peaks.
if isempty(pks_time)
    pks_time = round((interval_P300(2)-interval_P300(1))/2);
    tmp.nopeak = true;
end 

pks_time = interval_P300(1)+pks_time; % Correct index for peak time. 

tmp.time = timelockedMMN.avg(:,pks_time); % Save the values for all the sensors at the peak time

[tmp.maxval, maxch] = max(tmp.time, [], 1);
[tmp.minval, minch] = min(tmp.time, [], 1);

% Saving the biggest value and its index
if abs(tmp.maxval) > abs(tmp.minval)    
    tmp.peakch = timelockedMMN.label{maxch};
    tmp.amplitude = tmp.maxval;
    tmp.i_peakch = maxch;
else
    tmp.peakch = timelockedMMN.label{minch};
    tmp.amplitude = tmp.minval;
    tmp.i_peakch = minch;
end

%Save for statistical analysis
%peak.values(end+1,1) = tmp.amplitude;
%peak.labels(end+1,1) = {[params.modality '_' params.condition]};

% % Takes the choosen chs and finds the time with a peak
% [~, pks_time] = max(abs(timelockedMMN.avg(tmp.i_peakch, interval_P300(1):interval_P300(2))), [], 1); % Hitta tiden med bäst std. Använd tiden för att hitta relevanta peaks.
% pks_time = interval_P300(1) + pks_time;

% Plot
h = figure;
plot(timelockedMMN.time*1e3,timelockedMMN.avg*params.amp_scaler);
hold on
xline(timelockedMMN.time(pks_time)*10000, '--') % Choose timepoint, x-axis
hold off
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelocked.cfg.trials)) ')'])
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_MMN-' params.condition '.jpg']))

% Plotta max ch
timelockedMMN.avg = timelockedMMN.avg(tmp.i_peakch, :); % Choose peak ch, y-axis
params.condition = 'Sensor with the highest peak for MMN';
plot_butterfly(timelockedMMN, params, save_path)
save(fullfile(save_path, [params.sub '_' params.modality '_' params.condition]), 'timelocked', '-v7.3'); 


