function [peak] = freqtag_analysis(data, params, save_path, peak)
% % Det borde vara rätt channel 

cfg = [];
cfg.channel = params.chs;
data = ft_selectdata(cfg, data);

params.freqsteps = 0.2;
params.freqwin = [0.2 0.8];
% Timewindow for Topoplot in freqanalysis
params.pretimwin = 0.08;
params.posttimwin = 0.12;

%% Std tone
params.trials = find(data.trialinfo==params.trigger_code(1));
params.condition = 'Std';
params.freqylim = [37 41];
params.specificfreq = 39;

peak = freqanalysis(data, params, save_path, peak);

%% Low tone
params.trials = find(ismember(data.trialinfo, params.trigger_code(2:3)));
params.condition = 'Low';
params.freqylim = [41 45]; % Frequency limits for both high and low
params.specificfreq = 43;

peak = freqanalysis(data, params, save_path, peak);

%% High tone
params.trials = find(ismember(data.trialinfo, params.trigger_code(4:5)));
params.condition = 'High';
params.freqylim = [41 45]; % Frequency limits for both high and low
params.specificfreq = 43;

peak = freqanalysis(data, params, save_path, peak);

%% Go trigger
G = [];

GL = find(data.trialinfo==params.trigger_code(3));
GH = find(data.trialinfo==params.trigger_code(5));

G = [G; GL];
G = [G; GH];

params.trials = G;
params.condition = 'Go';
params.freqylim = [41 45];
params.specificfreq = 43;

peak = freqanalysis(data, params, save_path, peak);

%% Saving trigger indeces
pre_NG = [];
NG = find(ismember(data.trialinfo, [params.trigger_code(2) params.trigger_code(4)])); % all No Go triggers saved in one place

for i = 1:size(NG)
    pre_NG = [pre_NG; NG(i)-1]; % for every index we save the one before
end


%% Timelocking data
% pre No Go
params.trials = pre_NG;
params.condition = 'Std pre No Go';
params.freqylim = [37 41];
params.specificfreq = 39;

peak = freqanalysis(data, params, save_path, peak);

% No Go
params.trials = NG;
params.condition = 'No Go';
params.freqylim = [41 45];
params.specificfreq = 43;

peak = freqanalysis(data, params, save_path, peak);

end