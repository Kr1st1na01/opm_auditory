function [butterfly_plot] = plot_butterfly(data, params, save_path)

cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'prestim';
cfg.trials = params.trials; % data.trialinfo contains all data of the triggers so params.trigger_code(1) will put true (1) and false (0) on all the positions where the trigger code 1 appears, thus the function find will give every true (1) an index
timelocked = ft_timelockanalysis(cfg, data);
save(fullfile(save_path, [params.sub '_' params.modality '_timelocked_' params.condition]), 'timelocked', '-v7.3'); 

h = figure;
butterfly_plot = plot(timelocked.time*1e3,timelocked.avg*params.amp_scaler);
hold on
%ylimits = ylim;
%latency = 1e3*M100{code}.peak_latency;
%plot([latency latency],ylimits,'k--')
hold off
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - trigger ' params.condition ' (n_{trls}=' num2str(length(timelocked.cfg.trials)) ')'])
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_butterfly_audodd-' params.condition '.jpg']))
