function [timelocked_data] = timelock(data, params, save_path) % Here the data is timelocked

cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'prestim';
cfg.trials = params.trials;
timelocked_data = ft_timelockanalysis(cfg, data);
save(fullfile(save_path, [params.sub '_' params.modality '_timelocked_' params.condition]), 'timelocked_data', '-v7.3'); 

end