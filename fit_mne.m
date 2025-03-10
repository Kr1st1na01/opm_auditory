function [megmag_mne, megplanar_mne, opm_mne, megeeg_mne, opmeeg_mne] = fit_mne(save_path,megmag_timelocked,megplanar_timelocked,megeeg_timelocked,opm_timelocked, opmeeg_timelocked,headmodels,sourcemodel,latency,params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Prepare leadfields
cfg = [];
cfg.grad             = megmag_timelocked{1}.grad;              % sensor positions
cfg.channel          = 'megmag';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodels.headmodel_meg;          % volume conduction model
leadfield_megmag = ft_prepare_leadfield(cfg,megmag_timelocked{1});

cfg = [];
cfg.grad             = megplanar_timelocked{1}.grad;              % sensor positions
cfg.channel          = 'megplanar';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodels.headmodel_meg;          % volume conduction model
leadfield_megplanar = ft_prepare_leadfield(cfg,megplanar_timelocked{1});

cfg = [];
cfg.grad             = opm_timelocked{1}.grad;              % sensor positions
cfg.channel          = '*bz';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodels.headmodel_meg;          % volume conduction model
leadfield_opm = ft_prepare_leadfield(cfg,opm_timelocked{1});

if ~isempty(headmodels.headmodel_eeg)
    cfg = [];
    cfg.elec             = megeeg_timelocked{1}.elec;              % sensor positions
    cfg.channel          = 'eeg';                  % the used channels
    cfg.senstype         = 'eeg';            % sensor type
    cfg.grid.pos         = sourcemodel.pos;           % source points
    cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
    cfg.headmodel        = headmodels.headmodel_eeg;          % volume conduction model
    leadfield_megeeg = ft_prepare_leadfield(cfg,megeeg_timelocked{1});
    
    cfg = [];
    cfg.elec             = opmeeg_timelocked{1}.elec;              % sensor positions
    cfg.channel          = 'eeg';                  % the used channels
    cfg.senstype         = 'eeg';            % sensor type
    cfg.grid.pos         = sourcemodel.pos;           % source points
    cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
    cfg.headmodel        = headmodels.headmodel_eeg;          % volume conduction model
    leadfield_opmeeg = ft_prepare_leadfield(cfg,opmeeg_timelocked{1});
end

%% MNE invserse
% MEG-MAG
megmag_mne = [];
megmag_mne.avg = cell(5,1);
for i_phalange = 1:5
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodels.headmodel_meg;    % supply the headmodel
    %cfg.sourcemodel         = sourcemodel;
    %cfg.sourcemodel.leadfield = leadfield_megmag.leadfield;
    cfg.sourcemodel         = leadfield_megmag;
    cfg.senstype            = 'meg';            % sensor type
    cfg.channel             = 'megmag';         % which channels to use
    tmp = ft_sourceanalysis(cfg, megmag_timelocked{i_phalange});
    megmag_mne.avg{i_phalange} = [];
    megmag_mne.avg{i_phalange}.pow = tmp.avg.pow;
    megmag_mne.avg{i_phalange}.mom = tmp.avg.mom;
    megmag_mne.avg{i_phalange}.fahm = FullAreaHalfMax(megmag_mne,latency{i_phalange}.squidmag,i_phalange);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = latency{i_phalange}.squidmag;
    h = figure;
    ft_sourceplot(cfg, tmp)
    title(['SQUID-MAG (FAHM=' num2str(megmag_mne.avg{i_phalange}.fahm,3) ')'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_squidmag_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all
end
megmag_mne.time = tmp.time;
megmag_mne.cfg = tmp.cfg;
megmag_mne.method = tmp.method;
megmag_mne.pos = tmp.pos;
megmag_mne.tri = sourcemodel.tri;

save(fullfile(save_path, 'megmag_mne'), 'megmag_mne'); 
clear tmp megmag_mne leadfield_megmag

% MEG-GRAD
megplanar_mne = [];
megplanar_mne.avg = cell(5,1);
for i_phalange = 1:5
    cfg.channel             = 'meggrad';            % which channels to use
    %cfg.sourcemodel.leadfield = leadfield_megplanar.leadfield;
    cfg.sourcemodel         = leadfield_megplanar;
    tmp = ft_sourceanalysis(cfg, megplanar_timelocked{i_phalange});
    megplanar_mne.avg{i_phalange} = [];
    megplanar_mne.avg{i_phalange}.pow = tmp.avg.pow;
    megplanar_mne.avg{i_phalange}.mom = tmp.avg.mom;
    megplanar_mne.avg{i_phalange}.fahm = FullAreaHalfMax(megplanar_mne,latency{i_phalange}.squidgrad,i_phalange);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = latency{i_phalange}.squidgrad;
    h = figure;
    ft_sourceplot(cfg, tmp)
    title(['SQUID-GRAD (FAHM=' num2str(megplanar_mne.avg{i_phalange}.fahm,3) ')'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_squidgrad_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all
end
megplanar_mne.time = tmp.time;
megplanar_mne.cfg = tmp.cfg;
megplanar_mne.method = tmp.method;
megplanar_mne.pos = tmp.pos;
megplanar_mne.tri = sourcemodel.tri;

save(fullfile(save_path, 'megplanar_mne'), 'megplanar_mne'); 
clear tmp megplanar_mne leadfield_megplanar

% MEG-EEG
megeeg_mne = [];
megeeg_mne.avg = cell(5,1);
for i_phalange = 1:5
    megeeg_mne.avg{i_phalange} = [];
    if ~isempty(headmodels.headmodel_eeg)
        cfg.headmodel           = headmodels.headmodel_eeg;    % supply the headmodel
        cfg.senstype            = 'eeg';            % sensor type
        cfg.channel             = 'eeg';         % which channels to use
        %cfg.sourcemodel.leaddield = leadfield_megeeg.leadfield;
        cfg.sourcemodel         = leadfield_megeeg;
        tmp = ft_sourceanalysis(cfg, megeeg_timelocked{i_phalange});
        megeeg_mne.avg{i_phalange}.pow = tmp.avg.pow;
        megeeg_mne.avg{i_phalange}.mom = tmp.avg.mom;
        megeeg_mne.avg{i_phalange}.fahm = FullAreaHalfMax(megeeg_mne,latency{i_phalange}.squideeg,i_phalange);
        megeeg_mne.time = tmp.time;
        megeeg_mne.cfg = tmp.cfg;
        megeeg_mne.method = tmp.method;
        megeeg_mne.pos = tmp.pos;
        megeeg_mne.tri = sourcemodel.tri;

        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = latency{i_phalange}.squideeg;
        h = figure;
        ft_sourceplot(cfg, tmp)
        title(['SQUID-EEG (FAHM=' num2str(megeeg_mne.avg{i_phalange}.fahm,3) ')'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_megeeg_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
        close all
    end
end

save(fullfile(save_path, 'megeeg_mne'), 'megeeg_mne'); 
clear tmp megeeg_mne leadfield_megeeg

% OPM
opm_mne = [];
opm_mne.avg = cell(5,1);
for i_phalange = 1:5
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodels.headmodel_meg;    % supply the headmodel
    %cfg.sourcemodel         = sourcemodel;
    %cfg.sourcemodel.leadfield = leadfield_opm.leadfield;
    cfg.sourcemodel         = leadfield_opm;
    cfg.senstype            = 'meg';            % sensor type
    cfg.channel             = '*bz';         % which channels to use
    tmp = ft_sourceanalysis(cfg, opm_timelocked{i_phalange});
    opm_mne.avg{i_phalange} = [];
    opm_mne.avg{i_phalange}.pow = tmp.avg.pow;
    opm_mne.avg{i_phalange}.mom = tmp.avg.mom;
    opm_mne.avg{i_phalange}.fahm = FullAreaHalfMax(opm_mne,latency{i_phalange}.opm,i_phalange);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = latency{i_phalange}.opm;
    h = figure;
    ft_sourceplot(cfg, tmp)
    title(['OPM (FAHM=' num2str(opm_mne.avg{i_phalange}.fahm,3) ')'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_opm_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all
end
opm_mne.time = tmp.time;
opm_mne.cfg = tmp.cfg;
opm_mne.method = tmp.method;
opm_mne.pos = tmp.pos;
opm_mne.tri = sourcemodel.tri;

save(fullfile(save_path, 'opm_mne'), 'opm_mne'); 
clear tmp opm_mne leadfield_opm

% OPM-EEG
opmeeg_mne = [];
opmeeg_mne.avg = cell(5,1);
for i_phalange = 1:5
    opmeeg_mne.avg{i_phalange} = [];
    if ~isempty(headmodels.headmodel_eeg)
        cfg.headmodel           = headmodels.headmodel_eeg;    % supply the headmodel
        cfg.senstype            = 'eeg';            % sensor type
        cfg.channel             = 'eeg';         % which channels to use
        %cfg.sourcemodel.leadfield = leadfield_opmeeg.leadfield;
        cfg.sourcemodel         = leadfield_opmeeg;
        tmp = ft_sourceanalysis(cfg, opmeeg_timelocked{i_phalange});
        opmeeg_mne.avg{i_phalange}.pow = tmp.avg.pow;
        opmeeg_mne.avg{i_phalange}.mom = tmp.avg.mom;
        opmeeg_mne.avg{i_phalange}.fahm = FullAreaHalfMax(opmeeg_mne,latency{i_phalange}.opmeeg,i_phalange);
        opmeeg_mne.time = tmp.time;
        opmeeg_mne.cfg = tmp.cfg;
        opmeeg_mne.method = tmp.method;
        opmeeg_mne.pos = tmp.pos;
        opmeeg_mne.tri = sourcemodel.tri;
        
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = latency{i_phalange}.opmeeg;
        h = figure;
        ft_sourceplot(cfg, tmp)
        title(['OPM-EEG (FAHM=' num2str(opmeeg_mne.avg{i_phalange}.fahm,3) ')'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_opmeeg_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
        close all
    end
end

save(fullfile(save_path, 'opmeeg_mne'), 'opmeeg_mne'); 
clear tmp opmeeg_mne leadfield_opmeeg

end