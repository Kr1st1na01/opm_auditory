function [] = MMN(data, params, save_path)
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

params.trials = Oddball;
params.condition = 'Oddball trigger for No go';
timelocked = timelock(data, params, save_path);
plot_butterfly(timelocked, params, save_path)

%% MMN
% (Low No go + High No go) - (pre-Low No go + pre-High No go)
params.trials = Oddball;
params.condition = 'Trigger for No go vs pre-No go';

timelocked.avg = timelocked.avg - timelocked_pre.avg; % Here the oddball timelocked data average is changed to represent the difference (MMN) in response
plot_butterfly(timelocked, params, save_path)
save(fullfile(save_path, [params.sub '_' params.modality '_' params.condition]), 'timelocked', '-v7.3'); 


end