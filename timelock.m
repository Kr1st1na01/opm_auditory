function [timelocked_data] = timelock(data, params, save_path) % Here the data is timelocked

cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'prestim';
cfg.trials = params.trials;
tmp = ft_timelockanalysis(cfg, data);
cfg.channel = params.chs;
timelocked_data = ft_selectdata(cfg, tmp);
save(fullfile(save_path, [params.sub '_' params.modality '_timelocked_' params.condition]), 'timelocked_data', '-v7.3'); 

end