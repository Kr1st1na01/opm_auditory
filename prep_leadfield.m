function [leadfield] = prep_leadfield(headmodels,sourcemodel, timelocked_data, name, save_path)
% Prepare leadfields

headmodel = headmodels.headmodel_meg;
sourcemodel.mom = surface_normals(sourcemodel.pos, sourcemodel.tri, 'vertex')';
sourcemodel.unit = 'cm';

%% MEG MAG
    cfg = [];
    cfg.grad             = timelocked_data.grad; % sensor positions
    cfg.senstype         = 'meg';            % sensor type
    cfg.sourcemodel      = sourcemodel;           % source points
    cfg.headmodel        = headmodel;          % volume conduction model
    leadfield = ft_prepare_leadfield(cfg,timelocked_data);

save(fullfile(save_path,'source analysis', [name '_leadfield']), "leadfield"); disp('done');

end