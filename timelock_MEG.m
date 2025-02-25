function [timelocked] = timelock_MEG(data, save_path, params)

timelocked = cell(3,1); % cell(rader, columner)
% M100 = cell(5,1);
% h = figure; 
% hold on
% leg = [];   

cfg = [];
cfg.channel = params.chs;
data = ft_selectdata(cfg, data);

%% Timelock of Std

params.trials = find(data.trialinfo==params.trigger_code(1));
params.condition = 'Std';
timelocked = timelock(data, params, save_path);
plot_butterfly(timelocked, params, save_path)

% Sensor with highest peak
load([save_path '/' params.sub '_' params.modality '_timelocked_Std.mat']);
[~, pks_i] = max(max(abs(timelocked_data.avg), [], 2));
params.trials = timelocked_data.avg(pks_i,:);
params.condition = 'Sensor with highest peak for Std trigger';
plot_butterfly(timelocked, params, save_path)


%% Timelocked of Low

params.trials = find(ismember(data.trialinfo, params.trigger_code(2:3)));
params.condition = 'Low';
timelocked = timelock(data, params, save_path);
plot_butterfly(timelocked, params, save_path)

% Sensor with highest peak
load([save_path '/' params.sub '_' params.modality '_timelocked_Low.mat']);
[~, pks_i] = max(max(abs(timelocked_data.avg), [], 2));
params.trials = pks_i; 
params.condition = 'Sensor with highest peak for low trigger';
plot_butterfly(timelocked, params, save_path)

%% Timelocked of High

params.trials = find(ismember(data.trialinfo, params.trigger_code(4:5)));
params.condition = 'High';
timelocked = timelock(data, params, save_path);
plot_butterfly(timelocked, params, save_path)

% Sensor with highest peak
load([save_path '/' params.sub '_' params.modality '_timelocked_High.mat']);
[~, pks_i] = max(max(abs(timelocked_data.avg), [], 2));
params.trials = pks_i;
params.condition = 'Sensor with highest peak for high trigger';
plot_butterfly(timelocked, params, save_path)



%%
% leg = [];
% load([save_path '/' params.sub '_' params.modality '_timelocked_StdTone.mat']);
% dat = timelocked_data; % timelock averages the data
% [~, interval_M100(1)] = min(abs(dat.time-0.08)); % find closest time sample
% [~, interval_M100(2)] = min(abs(dat.time-0.125)); % find closest time sample
% [~, interval_M100(3)] = min(abs(dat.time-0)); % find closest time sample
% tmp = [];
% [~, l_peak_latency] = max((mean(abs(dat.avg(:,(interval_M100(1):interval_M100(2)))),1), SortStr='descend'), [], 2);
% disp(l_peak_latency);
% l_peak_latency = interval_M100(1)-1+l_peak_latency(1); % adjust for interval and pick first (=strongest) peak
% tmp.peak_latency = dat.time(l_peak_latency);
% disp(tmp.peak_latency);
% [tmp.max_amplitude, i_maxch] = max(dat.avg(:,l_peak_latency));
% [tmp.min_amplitude, i_minch] = min(dat.avg(:,l_peak_latency));
% tmp.max_channel = dat.label{i_maxch};
% tmp.min_channel = dat.label{i_minch};
% if abs(tmp.max_amplitude) > abs(tmp.min_amplitude)
%         tmp.peak_channel = tmp.max_channel;
%         tmp.peak_amplitude = abs(tmp.max_amplitude);
%         l_peakch = i_maxch;
%     else
%         tmp.peak_channel = tmp.min_channel;
%         tmp.peak_amplitude = abs(tmp.min_amplitude);
%         l_peakch = i_minch;
% end
% 
% %[~,i_peak_latency] = max(abs(dat.avg(i_maxch,interval_M100(1):interval_M100(2))));
% %tmp.peak_latency = dat.time(interval_M100(1)-1+i_peak_latency);
% tmp.prestim_std = std(dat.avg(l_peakch,1:interval_M100(3)));
% tmp.std_error = sqrt(dat.var(l_peakch,l_peak_latency));
% %for i_trl = find(data.trialinfo==params.trigger_code(i_trigger))'
% %   tmp.std_error = tmp.std_error + abs(tmp.max_amplitude - data.trial{i_trl}(i_maxch,interval_M100(1)-1+i_peak_latency));
% %end
% %tmp.std_error = tmp.std_error/length(find(data.trialinfo==params.trigger_code(i_trigger)));
% M100{or(params.trigger_code(2),params.trigger_code(3))} = tmp; % low
% 
% plot(dat.time*1e3, dat.avg(l_peakch,:)*params.amp_scaler)
% leg = [leg; [num2str(params.trigger_code(2),params.trigger_code(3)) ': ' strrep(tmp.peak_channel,'_','-')]];
% 
% hold off
% title([params.modality ' - Max channel'])
% ylabel(params.amp_label)
% xlabel('time [ms]')
% legend(leg)
% saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_evoked_maxchannels.jpg']))
% 
% save(fullfile(save_path, [params.sub '_' params.modality '_timelocked']), 'timelocked', '-v7.3'); 
% save(fullfile(save_path, [params.sub '_' params.modality '_M100']), 'M100', '-v7.3'); 
% 
% % Plot max channel with variation and peak time for low trigger
% dat = timelocked{or(params.trigger_code(2),params.trigger_code(3))};
% h = figure;
% hold on
% plot(data.time{data.trialinfo==find(or(params.trigger_code(2),params.trigger_code(3)))}*1e3, data.trial{(or(params.trigger_code(2),params.trigger_code(3)))}(i_peakch,:)*params.amp_scaler,'Color',[211 211 211]/255)
% plot(dat.time*1e3, dat.avg(i_peakch,:)*params.amp_scaler,'Color',[0 0 0]/255)
% ylimits = ylim;
% latency = 1e3*M100{length(params.trigger_code(2),params.trigger_code(3))}.peak_latency;
% plot([latency latency],ylimits,'r--')
% hold off
% title(['Peak ' params.modality ' channel: ' M100{length(params.trigger_code(2),params.trigger_code(3))}.peak_channel])
% ylabel(params.amp_label)
% xlabel('time [ms]')
% saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_evoked_peakchannel_ph' num2str(length(params.trigger_code(2),params.trigger_code(3))) '.jpg']))

% % Topoplot
% cfg = [];
% cfg.xlim = [M100{params.trigger_code > 1 & params.trigger_code < 8}.peak_latency-0.01 M100{params.trigger_code > 1 & params.trigger_code < 8}.peak_latency+0.01];
% cfg.layout = params.layout;
% cfg.parameter = 'avg';
% h = figure;
% ft_topoplotER(cfg, timelocked{params.trigger_code > 1 & params.trigger_code < 8});
% axis on
% colorbar
% saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_M100_topo_low-' params.trigger_labels{params.trigger_code > 1 & params.trigger_code < 8} '.jpg']))