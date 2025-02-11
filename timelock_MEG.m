function [timelocked] = timelock_MEG(data, save_path, params)

timelocked = cell(3,1); % cell(rader, columner)
% M100 = cell(5,1);
% h = figure; 
% hold on
% leg = [];   

cfg = [];
cfg.channel = params.chs;
data = ft_selectdata(cfg, data);

%% Saving trigger indeces of 3 and 11 and each Std trigger before
LowNG = [];
preLowNG = [];
params.trials = find(data.trialinfo==params.trigger_code(2)); % Finds index of trigger 3
for i = 1:length(params.trials)
    preLowNG = [preLowNG; params.trials(i)-1]; % Adds that preceding index
    LowNG = [LowNG; params.trials(i)]; % Adds index of trigger 3
end

HighNG = [];
preHighNG = [];
params.trials = find(data.trialinfo==params.trigger_code(4));
for i = 1:length(params.trials)
    preHighNG = [preHighNG; params.trials(i)-1];
    HighNG = [HighNG; params.trials(i)];
end

%% Timelock all the params.trials so I get 4 different timelocked trials (timelock is done when it is plotted).
params.trials = find(data.trialinfo==params.trigger_code(1));
params.condition = 'Std';
timelocked = timelock(data, params, save_path);
plot_butterfly(data, timelocked, params, save_path)

params.trials = LowNG;
params.condition = 'Low No Go';
timelocked = timelock(data, params, save_path);
plot_butterfly(data, timelocked, params, save_path)

params.trials = preLowNG;
params.condition = 'pre Low No Go';
timelocked = timelock(data, params, save_path);
plot_butterfly(data, timelocked, params, save_path)

params.trials = HighNG;
params.condition = 'High No Go';
timelocked = timelock(data, params, save_path);
plot_butterfly(data, timelocked, params, save_path)

params.trials = preHighNG;
params.condition = 'pre High No Go';
timelocked = timelock(data, params, save_path);
plot_butterfly(data, timelocked, params, save_path)

%% MMN
%Low_MNN = load(params.sub '_' params.modality '_timelocked_pre Low No Go') - load(params.sub '_' params.modality '_timelocked_Low No Go'); % Load the low tone file and pre low and subtract them from each other
butterfly_plot(Low_MNN, params, save_path)

%High_MNN = load(params.sub '_' params.modality '_timelocked_High No Go'); % Load the pre High tone file
%load(params.sub '_' params.modality '_timelocked_pre High No Go'); % Load the pre High tone file
butterfly_plot(High_MNN, params, save_path)

%% Normal trigger
params.trials = find(data.trialinfo==params.trigger_code(1));
params.condition = 'stdTone';
plot_butterfly(data, params, save_path)

% dat = timelocked{params.trigger_code(1)}; % timelock averages the data
% [~, interval_M100(1)] = min(abs(dat.time-0.08)); % find closest time sample
% [~, interval_M100(2)] = min(abs(dat.time-0.125)); % find closest time sample
% [~, interval_M100(3)] = min(abs(dat.time-0)); % find closest time sample
% tmp = [];
% [~, i_peak_latency] = findpeaks(mean(abs(dat.avg(:,(interval_M100(1):interval_M100(2)))),1),'SortStr','descend');
% disp(i_peak_latency);
% i_peak_latency = interval_M100(1)-1+i_peak_latency(1); % adjust for interval and pick first (=strongest) peak
% tmp.peak_latency = dat.time(i_peak_latency);
% disp(tmp.peak_latency);
% [tmp.max_amplitude, i_maxch] = max(dat.avg(:,i_peak_latency));
% [tmp.min_amplitude, i_minch] = min(dat.avg(:,i_peak_latency));
% tmp.max_channel = dat.label{i_maxch};
% tmp.min_channel = dat.label{i_minch};
% if abs(tmp.max_amplitude) > abs(tmp.min_amplitude) % compare the amplitudes of the min and max to see the largest
%         tmp.peak_channel = tmp.max_channel;
%         tmp.peak_amplitude = abs(tmp.max_amplitude);
%         n_peakch = i_maxch;
%     else
%         tmp.peak_channel = tmp.min_channel;
%         tmp.peak_amplitude = abs(tmp.min_amplitude);
%         n_peakch = i_minch;
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
% % Plot max channel with variation and peak time for high trigger
% dat = timelocked{params.trigger_code(1)};
% h = figure;
% hold on
% plot(data.time{data.trialinfo==find(params.trigger_code(1))}*1e3, data.trial{params.trigger_code(1)}(i_peakch,:)*params.amp_scaler,'Color',[211 211 211]/255)
% plot(dat.time*1e3, dat.avg(i_peakch,:)*params.amp_scaler,'Color',[0 0 0]/255)
% ylimits = ylim;
% latency = 1e3*M100{length(params.trigger_code(1))}.peak_latency;
% plot([latency latency],ylimits,'r--')
% hold off
% title(['Peak ' params.modality ' channel: ' M100{length(params.trigger_code(1))}.peak_channel])
% ylabel(params.amp_label)
% xlabel('time [ms]')
% saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_evoked_peakchannel_ph' num2str(length(params.trigger_code(1))) '.jpg']))

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

%% High trigger
% Butterfly plot
params.trials = find(ismember(data.trialinfo, params.trigger_code(4:5))); %Första är vid nr 17 (trigger 13)
params.condition = 'HighTone';
plot_butterfly(data, params, save_path)

% [~, interval_M100(1)] = min(abs(dat.time-0.08)); % find closest time sample
% [~, interval_M100(2)] = min(abs(dat.time-0.125)); % find closest time sample
% [~, interval_M100(3)] = min(abs(dat.time-0)); % find closest time sample
% tmp = [];
% [~, i_peak_latency] = findpeaks(mean(abs(dat.avg(:,(interval_M100(1):interval_M100(2)))),1),'SortStr','descend');
% disp(i_peak_latency);
% i_peak_latency = interval_M100(1)-1+i_peak_latency(1); % adjust for interval and pick first (=strongest) peak
% tmp.peak_latency = dat.time(i_peak_latency);
% disp(tmp.peak_latency);
% [tmp.max_amplitude, i_maxch] = max(dat.avg(:,i_peak_latency));
% [tmp.min_amplitude, i_minch] = min(dat.avg(:,i_peak_latency));
% tmp.max_channel = dat.label{i_maxch};
% tmp.min_channel = dat.label{i_minch};
% if abs(tmp.max_amplitude) > abs(tmp.min_amplitude)
%         tmp.peak_channel = tmp.max_channel;
%         tmp.peak_amplitude = abs(tmp.max_amplitude);
%         h_peakch = i_maxch;
%     else
%         tmp.peak_channel = tmp.min_channel;
%         tmp.peak_amplitude = abs(tmp.min_amplitude);
%         h_peakch = i_minch;
%end
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
% % Plot max channel with variation and peak time for high trigger
% dat = timelocked{or(params.trigger_code(4),params.trigger_code(5))};
% h = figure;
% hold on
% plot(data.time{data.trialinfo==find(or(params.trigger_code(4),params.trigger_code(5)))}*1e3, data.trial{(or(params.trigger_code(4),params.trigger_code(5)))}(i_peakch,:)*params.amp_scaler,'Color',[211 211 211]/255)
% plot(dat.time*1e3, dat.avg(i_peakch,:)*params.amp_scaler,'Color',[0 0 0]/255)
% ylimits = ylim;
% latency = 1e3*M100{length(params.trigger_code(4),params.trigger_code(5))}.peak_latency;
% plot([latency latency],ylimits,'r--')
% hold off
% title(['Peak ' params.modality ' channel: ' M100{length(params.trigger_code(4),params.trigger_code(5))}.peak_channel])
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


%% Low trigger
% Butterfly plot
params.trials = find(ismember(data.trialinfo, params.trigger_code(2:3)));
params.condition = 'LowTone';
plot_butterfly(data, params, save_path)


% [~, interval_M100(1)] = min(abs(dat.time-0.08)); % find closest time sample
% [~, interval_M100(2)] = min(abs(dat.time-0.125)); % find closest time sample
% [~, interval_M100(3)] = min(abs(dat.time-0)); % find closest time sample
% tmp = [];
% [~, l_peak_latency] = findpeaks(mean(abs(dat.avg(:,(interval_M100(1):interval_M100(2)))),1),'SortStr','descend');
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