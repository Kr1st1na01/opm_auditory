function fit_mne(save_path, timelocked_data, headmodel, sourcemodel, sourcemodel_inflated, leadfield, params)


%% MNE invserse
MNE = [];

    if isfield(params,'use_cov_all') 
        if params.use_cov_all
            timelocked_data.cov = timelocked_data.cov_all;
        end
    end

    cfg = [];
    cfg.method              = params.method;
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.senstype            = 'meg';
    cfg.sourcemodel         = leadfield;
    cfg.keepfilter          = 'yes';
    tmp = ft_sourceanalysis(cfg, timelocked_data);
    tmp.tri = sourcemodel.tri;
    
    MNE = FullAreaHalfMax(tmp, sourcemodel, params, save_path);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = squidgrad_mne_MMN.peak_latency; %%%%%%%
    tmp.pos = sourcemodel_inflated.pos;
    tmp.tri = sourcemodel_inflated.tri;
    h = figure;
    ft_sourceplot(cfg, tmp)
    lighting gouraud
    material dull
    title(['SQUID-GRAD (FAHM=' num2str(MNE.fahm,3) 'cm^2; t=' num2str(round(squidgrad_mne_MMN.peak_latency*1e3)) 'ms)']) %%%%%%%%%%%
    saveas(h, fullfile(save_path,'source analysis', [params.sub '_squidgrad_mne_' params.trigger_labels '.jpg']))
    close all

    clear tmp



%     %% Overlaps (see how much the triangles, mesh, overlaps between
%     modalities
%     i_vertices = opm_mne_M60{i_phalange}.halfmax_distribution & squidmag_mne_M60{i_phalange}.halfmax_distribution;
%     [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
%     triangles = sourcemodel.tri(triangles,:);
%     opm_mne_M60{i_phalange}.overlap_squidmag = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
%     squidmag_mne_M60{i_phalange}.overlap_opm = opm_mne_M60{i_phalange}.overlap_squidmag;
% 
%     i_vertices = opm_mne_M60{i_phalange}.halfmax_distribution & squidgrad_mne_M60{i_phalange}.halfmax_distribution;
%     [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
%     triangles = sourcemodel.tri(triangles,:);
%     opm_mne_M60{i_phalange}.overlap_squidgrad = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
%     squidgrad_mne_M60{i_phalange}.overlap_opm = opm_mne_M60{i_phalange}.overlap_squidgrad;
% 
%     i_vertices = squidmag_mne_M60{i_phalange}.halfmax_distribution & squidgrad_mne_M60{i_phalange}.halfmax_distribution;
%     [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
%     triangles = sourcemodel.tri(triangles,:);
%     squidmag_mne_M60{i_phalange}.overlap_squidgrad = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
%     squidgrad_mne_M60{i_phalange}.overlap_squidmag = squidmag_mne_M60{i_phalange}.overlap_squidgrad;


save(fullfile(save_path, 'squidmag_mne_M60'), 'squidmag_mne_M60'); 


end