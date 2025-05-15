function [ylimit] = plot_butterfly(timelocked, params, save_path, ylimit)

if contains('opm', params.modality) % Save the ylim for OPMs
    h = figure;
    ax = axes(h);
    plot(ax, timelocked.time*1e3,timelocked.avg*params.amp_scaler);
    % Saving ylim
    y = max(abs(get(ax, "YLim")));
    ylimit(end+1,:) = {params.condition, y};
    if any(strcmp(ylimit{:,1}, replace(params.condition, '_', ' ')))
        idx = strcmp(ylimit(:,1), params.condition);
        ylim(ax, [-ylimit{idx,2}, ylimit{idx,2}])
    end
    xlabel('t [msec]')
    ylabel(params.amp_label)
    title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelocked.cfg.trials)) ')'], Interpreter="none")
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Evoked data_butterfly audodd_' params.condition '.jpg']))

elseif contains('squidmag', params.modality) % Plot the correct ylim for squidmags
    h = figure;
    plot(timelocked.time*1e3,timelocked.avg*params.amp_scaler);
    for i = 1:length(ylimit)
        if strcmp(ylimit{i}, params.condition)
%             ysquid = get(gca, 'Ylim');
%             Amax = max(ylimit{i,2}(2), ysquid(2));
%             Amin = min(ylimit{i,2}(1), ysquid(1));            
%             ylim([Amin,Amax]);
            ylim = ([-ylimit{i,2} ylimit{i,2}]);
        end
    end
    xlabel('t [msec]')
    ylabel(params.amp_label)
    title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelocked.cfg.trials)) ')'], Interpreter="none")
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Evoked data_butterfly audodd_' params.condition '.jpg']))

else % all other plots are plotted as normal
    h = figure;
    plot(timelocked.time*1e3,timelocked.avg*params.amp_scaler);
    xlabel('t [msec]')
    ylabel(params.amp_label)
    title(['Evoked ' params.modality ' - ' params.condition ' (n_{trls}=' num2str(length(timelocked.cfg.trials)) ')'], Interpreter="none")
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_Evoked data_butterfly audodd_' params.condition '.jpg']))
end