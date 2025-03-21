function [butterfly_plot] = plot_butterfly(timelocked, params, save_path)

h = figure;
butterfly_plot = plot(timelocked.time*1e3,timelocked.avg*params.amp_scaler);
xlabel('t [msec]')
ylabel(params.amp_label)
title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelocked.cfg.trials)) ')'], Interpreter="none")
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_MMN data_butterfly audodd_' params.condition '.jpg']))
