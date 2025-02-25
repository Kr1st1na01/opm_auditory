function [] = MMN(data, params, save_path)
% Här gör vi MMN

%% Saving trigger indeces of 3 and 11 and each Std trigger before
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

%% Timelocking

params.trials = preOddball;
params.condition = 'Trigger before oddball for No Go';
timelocked = timelock(data, params, save_path);
plot_butterfly(timelocked, params, save_path)

params.trials = Oddball;
params.condition = 'Oddball trigger for No go';
timelocked = timelock(data, params, save_path);
plot_butterfly(timelocked, params, save_path)

%% MMN
% (Low No go + High No go) - (pre-Low No go + pre-High No go)
params.trials = Oddball;
params.condition = 'Trigger for No go vs pre-No go MMN';

load([save_path '/' params.sub '_' params.modality '_timelocked_Trigger before oddball for No Go']);
MMN_pre = timelocked_data; % the loaded file is saved in a variable
load([save_path '/' params.sub '_' params.modality '_timelocked_Oddball trigger for No go.mat']);
MMN = timelocked_data;

MMN.filter = params.filter;
MMN.avg = MMN.avg - MMN_pre.avg; % Here the oddball timelocked data average is changed to represent the difference (MMN) in response
plot_butterfly(MMN, params, save_path)
save(fullfile(save_path, [params.sub '_' params.modality '_' params.condition]), 'MMN', '-v7.3'); 


end