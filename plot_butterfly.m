function [ylimit] = plot_butterfly(timelocked, params, save_path, ylimit)

if contains('opm', params.modality)
    h = figure;
    plot(timelocked.time*1e3,timelocked.avg*params.amp_scaler);
    y = max(abs(ylim));
    ylimit(end+1,:) = {params.condition, y};
    if contains(params.condition, ylimit(:,1))
        name = split(replace(params.condition, '_', ' '));
        idx = find(strcmp(ylimit(:,1), name{end}));
        ylim([-ylimit{idx,2} ylimit{idx,2}])
    end
    xlabel('t [msec]')
    ylabel(params.amp_label)
    title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelocked.cfg.trials)) ')'], Interpreter="none")
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Evoked data_butterfly audodd_' params.condition '.jpg']))
elseif contains('squidmag', params.modality)
    h = figure;
    plot(timelocked.time*1e3,timelocked.avg*params.amp_scaler);
    if contains(params.condition, ylimit(:,1))
        name = split(replace(params.condition, '_', ' '));
        idx = find(strcmp(ylimit(:,1), name{end}));
        ylim([-ylimit{idx,2} ylimit{idx,2}])
    end
    xlabel('t [msec]')
    ylabel(params.amp_label)
    title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelocked.cfg.trials)) ')'], Interpreter="none")
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Evoked data_butterfly audodd_' params.condition '.jpg']))
else
    h = figure;
    plot(timelocked.time*1e3,timelocked.avg*params.amp_scaler);
    xlabel('t [msec]')
    ylabel(params.amp_label)
    title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelocked.cfg.trials)) ')'], Interpreter="none")
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Evoked data_butterfly audodd_' params.condition '.jpg']))
end