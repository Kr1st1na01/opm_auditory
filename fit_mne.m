function fit_mne(save_path, timelocked_data, headmodel, sourcemodel, sourcemodel_inflated, leadfield, params, peak)


%% MNE invserse
MNE = [];

    if isfield(params,'use_cov_all') 
        if params.use_cov_all
            timelocked_data.cov = timelocked_data.cov_all;
        end
    end

    cfg = [];
    cfg.method              = 'mne';
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

    h = figure;
    plot(tmp.time, tmp.avg.pow)
    title(params.modality)

    h = figure;
    plot(timelocked_data.time, timelocked_data.avg)
    title([params.modality '_timelocked data'], Interpreter="none")

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = peak; 
    tmp.pos = sourcemodel_inflated.pos;
    tmp.tri = sourcemodel_inflated.tri;
    h = figure;
    ft_sourceplot(cfg, tmp)
    lighting gouraud
    material dull
    view(-90,0)
    title([params.modality ' (FAHM=' num2str(MNE.fahm,3) 'cm^2; t=' num2str(round(peak*1e3)) 'ms)']) 
    saveas(h, fullfile(save_path,'source analysis', [params.sub '_' params.modality '_mne_.jpg']))
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


%save(fullfile(save_path, 'squidmag_mne_M60'), 'squidmag_mne_M60'); 


end