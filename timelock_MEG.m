function [timelocked] = timelock_MEG(data, save_path, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
timelocked = cell(3,1); % cell(rader, columner)
M100 = cell(5,1);
h = figure; 
hold on
leg = [];

cfg = [];
cfg.channel = params.chs;
data = ft_selectdata(cfg, data);

%% Normal trigger
cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'prestim';
cfg.trials = find(data.trialinfo==params.trigger_code(1)); % data.trialinfo contains all data of the triggers so params.trigger_code(1) will put true (1) and false (0) on all the positions where the trigger code 1 appears, thus the function find will give every true (1) an index
timelocked{params.trigger_code == 1} = ft_timelockanalysis(cfg, data);
dat = timelocked{params.trigger_code == 1}; % timelock averages the data

% Butterfly plot
chs = find(contains(timelocked{params.trigger_code == 1}.label,ft_channelselection(params.chs,timelocked{params.trigger_code > 1 & params.trigger_code < 8}.label)));
h = figure;
plot(timelocked{params.trigger_code(1)}.time*1e3,timelocked{params.trigger_code(1)}.avg(chs,:)*params.amp_scaler)
hold on
ylimits = ylim;
%latency = 1e3*M100{params.trigger_code > 1 & params.trigger_code < 8}.peak_latency;
%plot([latency latency],ylimits,'k--')
hold off
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - trigger ' params.trigger_labels{params.trigger_code(1)} ' (n_{trls}=' num2str(length(timelocked{params.trigger_code(1)}.cfg.trials)) ')'])
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_butterfly_ph-' params.trigger_labels{params.trigger_code(1)} '.jpg']))

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


%% High trigger
cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'prestim';
cfg.trials = find(or(data.trialinfo==params.trigger_code(11),data.trialinfo==params.trigger_code(13)));
timelocked{params.trigger_code > 8} = ft_timelockedanalysis(cfg, data);
dat = timelocked{params.trigger_code > 8}; % timelock averages the data

% Butterfly plot
chs = find(contains(timelocked{params.trigger_code > 8}.label,ft_channelselection(params.chs,timelocked{params.trigger_code > 8}.label)));
h = figure;
plot(timelocked{params.trigger_code > 8}.time*1e3,timelocked{params.trigger_code > 8}.avg(chs,:)*params.amp_scaler)
hold on
ylimits = ylim;
%latency = 1e3*M100{params.trigger_code > 1 & params.trigger_code < 8}.peak_latency;
%plot([latency latency],ylimits,'k--')
hold off
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - trigger ' params.trigger_labels{params.trigger_code > 8} ' (n_{trls}=' num2str(length(timelocked{params.trigger_code > 8}.cfg.trials)) ')'])
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_butterfly_ph-' params.trigger_labels{params.trigger_code > 8} '.jpg']))


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
% end


%% Low trigger
cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'prestim';
cfg.trials = find(or(data.trialinfo==params.trigger_code(2),data.trialinfo==params.trigger_code(3)));
timelocked{params.trigger_code > 1 & params.trigger_code < 8} = ft_timelockedanalysis(cfg, data);
dat = timelocked{params.trigger_code > 1 & params.trigger_code < 8}; % timelock averages the data

% Butterfly plot
chs = find(contains(timelocked{params.trigger_code > 1 & params.trigger_code < 8}.label,ft_channelselection(params.chs,timelocked{params.trigger_code > 1 & params.trigger_code < 8}.label)));
h = figure;
plot(timelocked{params.trigger_code > 1 & params.trigger_code < 8}.time*1e3,timelocked{params.trigger_code > 1 & params.trigger_code < 8}.avg(chs,:)*params.amp_scaler)
hold on
ylimits = ylim;
%latency = 1e3*M100{params.trigger_code > 1 & params.trigger_code < 8}.peak_latency;
%plot([latency latency],ylimits,'k--')
hold off
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - trigger ' params.trigger_labels{params.trigger_code > 1 & params.trigger_code < 8} ' (n_{trls}=' num2str(length(timelocked{params.trigger_code > 1 & params.trigger_code < 8}.cfg.trials)) ')'])
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_butterfly_ph-' params.trigger_labels{params.trigger_code > 1 & params.trigger_code < 8} '.jpg']))

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
    

    %[~,i_peak_latency] = max(abs(dat.avg(i_maxch,interval_M100(1):interval_M100(2))));
    %tmp.peak_latency = dat.time(interval_M100(1)-1+i_peak_latency);
    tmp.prestim_std = std(dat.avg(l_peakch,1:interval_M100(3)));
    tmp.std_error = sqrt(dat.var(l_peakch,l_peak_latency));
    %for i_trl = find(data.trialinfo==params.trigger_code(i_trigger))'
    %    tmp.std_error = tmp.std_error + abs(tmp.max_amplitude - data.trial{i_trl}(i_maxch,interval_M100(1)-1+i_peak_latency));
    %end
    %tmp.std_error = tmp.std_error/length(find(data.trialinfo==params.trigger_code(i_trigger)));
    M100{params.trigger_code > 1 & params.trigger_code < 8} = tmp; % low
    
    plot(dat.time*1e3, dat.avg(l_peakch,:)*params.amp_scaler)
    leg = [leg; [num2str(params.trigger_code > 1 & params.trigger_code < 8) ': ' strrep(tmp.peak_channel,'_','-')]];


hold off
title([params.modality ' - Max channel'])
ylabel(params.amp_label)
xlabel('time [ms]')
legend(leg)
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_evoked_maxchannels.jpg']))

save(fullfile(save_path, [params.sub '_' params.modality '_timelocked']), 'timelocked', '-v7.3'); 
save(fullfile(save_path, [params.sub '_' params.modality '_M100']), 'M100', '-v7.3'); 

%% Plot max channel with variation and peak time
dat = timelocked{params.trigger_code > 1 & params.trigger_code < 8};
h = figure;
hold on
for i_trl = find(data.trialinfo == (params.trigger_code > 1 & params.trigger_code < 8))'
    plot(data.time{i_trl}*1e3, data.trial{i_trl}(l_peakch,:)*params.amp_scaler,'Color',[211 211 211]/255)
end
plot(dat.time*1e3, dat.avg(l_peakch,:)*params.amp_scaler,'Color',[0 0 0]/255)
ylimits = ylim;
latency = 1e3*M100{params.trigger_code > 1 & params.trigger_code < 8}.peak_latency;
plot([latency latency],ylimits,'r--')
hold off
title(['Peak ' params.modality ' channel: ' M100{params.trigger_code > 1 & params.trigger_code < 8}.peak_channel])
ylabel(params.amp_label)
xlabel('time [ms]')
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_evoked_peakchannel_l' num2str(params.trigger_code > 1 & params.trigger_code < 8) '.jpg']))

%% Butterfly & topoplot
chs = find(contains(timelocked{params.trigger_code > 1 & params.trigger_code < 8}.label,ft_channelselection(params.chs,timelocked{params.trigger_code > 1 & params.trigger_code < 8}.label)));
h = figure;
plot(timelocked{params.trigger_code > 1 & params.trigger_code < 8}.time*1e3,timelocked{params.trigger_code > 1 & params.trigger_code < 8}.avg(chs,:)*params.amp_scaler)
hold on
ylimits = ylim;
%latency = 1e3*M100{params.trigger_code > 1 & params.trigger_code < 8}.peak_latency;
%plot([latency latency],ylimits,'k--')
hold off
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - trigger ' params.trigger_labels{params.trigger_code > 1 & params.trigger_code < 8} ' (n_{trls}=' num2str(length(timelocked{params.trigger_code > 1 & params.trigger_code < 8}.cfg.trials)) ')'])
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_butterfly_ph-' params.trigger_labels{params.trigger_code > 1 & params.trigger_code < 8} '.jpg']))

cfg = [];
cfg.xlim = [M100{params.trigger_code > 1 & params.trigger_code < 8}.peak_latency-0.01 M100{params.trigger_code > 1 & params.trigger_code < 8}.peak_latency+0.01];
cfg.layout = params.layout;
cfg.parameter = 'avg';
h = figure;
ft_topoplotER(cfg, timelocked{params.trigger_code > 1 & params.trigger_code < 8});
axis on
colorbar
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_M100_topo_low-' params.trigger_labels{params.trigger_code > 1 & params.trigger_code < 8} '.jpg']))
