function [squidmag_dipole, squidgrad_dipole, opm_dipole] = fit_dipoles(save_path,squidmag_timelocked,squidgrad_timelocked,opm_timelocked,headmodels,mri,MMN_squidmag, MMN_squidgrad, MMN_opm, params, peak)



colors = [[0 0.4470 0.7410]; % blue
    [0.8500 0.3250 0.0980]; % red
    [0.9290 0.6940 0.1250]; % yellow
    [0.4940 0.1840 0.5560]; % purple
    [0.4660 0.6740 0.1880]; % green
    [0.6350 0.0780 0.1840]]; % light blue

cfg              = [];
cfg.resolution   = 1;
cfg.tight        = 'yes';
cfg.inwardshift  = 0;
cfg.headmodel    = headmodels.headmodel_meg;
sourcemodel    = ft_prepare_sourcemodel(cfg);

%% Fit dipoles
    % MEG
    cfg = [];
    cfg.gridsearch      = 'yes';           
    cfg.numdipoles      = 1;                
    cfg.sourcemodel     = sourcemodel;           
    cfg.headmodel       = headmodels.headmodel_meg;    
    cfg.senstype        = 'meg';            
    cfg.channel         = 'megmag';         
    cfg.nonlinear       = 'yes';           
    cfg.latency         = MMN_squidmag.time(1,401) + [-0.01 0.01];
    cfg.dipfit.checkinside = 'yes';
    %cfg.dipfit.noisecov = meg_timelocked{i_phalange}.cov;
    squidmag_dipole = ft_dipolefitting(cfg, squidmag_timelocked);
    
    cfg.latency         = MMN.squidgrad + [-0.01 0.01];   
    cfg.channel         = 'megplanar';           
    squidgrad_dipole = ft_dipolefitting(cfg, squidgrad_timelocked);

    if ~isempty(headmodels.headmodel_eeg)
        cfg.headmodel       = headmodels.headmodel_eeg;    
        cfg.senstype        = 'eeg';            
        cfg.channel         = 'eeg';   
        cfg.latency         = latency{i_phalange}.squideeg + [-0.01 0.01];   
        squideeg_dipole = ft_dipolefitting(cfg, squideeg_timelocked);
    else 
        squideeg_dipole = [];
    end

    % OPM
    cfg = [];
    cfg.gridsearch      = 'yes';           
    cfg.numdipoles      = 1;                
    cfg.sourcemodel     = sourcemodel;            
    cfg.headmodel       = headmodels.headmodel_meg;    
    cfg.senstype        = 'meg';            
    cfg.channel         = '*bz';        
    cfg.nonlinear       = 'yes';            
    cfg.latency         = latency.opm + [-0.01 0.01];   
    cfg.dipfit.checkinside = 'yes';
    %cfg.dipfit.noisecov = opm_timelocked{i_phalange}.cov;
    opm_dipole = ft_dipolefitting(cfg, opm_timelocked);

    if ~isempty(headmodels.headmodel_eeg)
        cfg.headmodel       = headmodels.headmodel_eeg;  
        cfg.senstype        = 'eeg';           
        cfg.channel         = 'eeg';        
        cfg.latency         = latency{i_phalange}.opmeeg + [-0.01 0.01];   
        opmeeg_dipole = ft_dipolefitting(cfg, opmeeg_timelocked);
    else 
        opmeeg_dipole = [];
    end

    % Plot OPM vs SQUID
    pos_mag = squidmag_dipole.dip.pos;
    [~,idx] = max(vecnorm(squidmag_dipole.dip.mom,2,1));
    ori_mag = squidmag_dipole.dip.mom(:,idx);

    pos_gad = squidgrad_dipole.dip.pos;
    [~,idx] = max(vecnorm(squidgrad_dipole.dip.mom,2,1));
    ori_grad = squidgrad_dipole.dip.mom(:,idx);
    
    pos_opm = opm_dipole.dip.pos;
    %pos_opm = opm_trans.transformPointsInverse(pos_opm);
    [~,idx] = max(vecnorm(opm_dipole.dip.mom,2,1));
    ori_opm = -opm_dipole.dip.mom(:,idx);

    if ~isempty(headmodels.headmodel_eeg)
        pos_eeg = squideeg_dipole.dip.pos;
        [~,idx] = max(vecnorm(squideeg_dipole.dip.mom,2,1));
        ori_eeg = squideeg_dipole.dip.mom(:,idx);
    
        pos_opmeeg = opmeeg_dipole.dip.pos;
        %pos_opmeeg = opm_trans.transformPointsInverse(pos_opmeeg);
        [~,idx] = max(vecnorm(opmeeg_dipole.dip.mom,2,1));
        ori_opmeeg = -opmeeg_dipole.dip.mom(:,idx);

        h = figure;
        ft_plot_dipole(pos_eeg,ori_eeg,'color',colors(1,:))
        hold on;
        ft_plot_dipole(pos_opmeeg,ori_opmeeg,'color',colors(2,:))
        ft_plot_headmodel(headmodels.headmodel_meg,'EdgeAlpha',0,'FaceAlpha',0.3,'FaceColor',[229 194 152]/256,'unit','cm') 
        hold off
        title([params.trigger_labels ' (SQEEG-OPMEEG = ' num2str(norm((pos_eeg-pos_opmeeg))*10,'%.1f') 'mm)'])
        legend('SQUIDEEG','OPMEEG','brain')
        saveas(h, fullfile(save_path, 'figs', [params.sub '_dipfit_SQUIDEEGvOPMEEG_ph-' params.trigger_labels '.jpg']))
        close
    end

    h = figure;
    ft_plot_dipole(pos_mag,ori_mag,'color',colors(1,:))
    hold on;
    ft_plot_dipole(pos_opm,ori_opm,'color',colors(2,:))
    ft_plot_dipole(pos_gad,ori_grad,'color',colors(3,:))
    ft_plot_headmodel(headmodels.headmodel_meg,'EdgeAlpha',0,'FaceAlpha',0.3,'FaceColor',[229 194 152]/256,'unit','cm') 
    hold off
    title([params.phalange_labels{i_phalange} ' (SQMAG-OPM = ' num2str(norm(pos_mag-pos_opm)*10,'%.1f') 'mm / SQGRAD-OPM = ' num2str(norm(pos_gad-pos_opm)*10,'%.1f') 'mm)'])
    legend('SQUIDMAG','OPM','SQUIDPLANAR','brain')
    saveas(h, fullfile(save_path, 'figs', [params.sub '_dipfit_SQUIDvOPM_ph-' params.trigger_lab '.jpg']))
    close

close all

%% Plot phalanges jointly
% SQUID
params.modality = 'squidmag';
pos_mag = zeros(5,3);
ori_mag = zeros(5,3);
for i = 1:5
    pos_mag(i,:) = squidmag_dipole{i}.dip.pos;
    [~,idx] = max(vecnorm(squidmag_dipole{i}.dip.mom,2,1));
    ori_mag(i,:) = squidmag_dipole{i}.dip.mom(:,idx);
end
mean_pos = mean(pos_mag,1);
tmp = pos_mag;
tmp(:,1) = mean_pos(1)-1;
%100
h=figure;
hold on
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'resolution', 0.1)
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_mag(i,:),'color',colors(i,:))
end
hold off
view(-90,0)
title('SQUID-MAG')
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_dipfit-100.jpg']))
%010
h=figure;
hold on
tmp = pos_mag;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'resolution', 0.1)
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_mag(i,:),'color',colors(i,:))
end
hold off
view(0,0)
title('SQUID-MAG')
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_dipfit-010.jpg']))
%001
h=figure;
hold on
tmp = pos_mag;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'resolution', 0.1)       
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_mag(i,:),'color',colors(i,:))
end
hold off
title('SQUID-MAG')
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_dipfit-001.jpg']))
close all

% OPM
params.modality = 'opm';
pos_opm = zeros(5,3);
ori_opm = zeros(5,3);
for i = 1:5
    pos_opm(i,:) = opm_dipole{i}.dip.pos;
    [~,idx] = max(vecnorm(opm_dipole{i}.dip.mom,2,1));
    ori_opm(i,:) = -opm_dipole{i}.dip.mom(:,idx);
end
mean_pos = mean(pos_opm,1);
tmp = pos_opm;
tmp(:,1) = mean_pos(1)-1;
%100
h=figure;
hold on
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'resolution', 0.1)
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_opm(i,:),'color',colors(i,:))
end
hold off
view(-90,0)
title('OPM')
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_dipfit-100.jpg']))
%010
h=figure;
hold on
tmp = pos_opm;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'resolution', 0.1)
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_opm(i,:),'color',colors(i,:))
end
hold off
view(0,0)
title('OPM')
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_dipfit-010.jpg']))
%001
h=figure;
hold on
tmp = pos_opm;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'resolution', 0.1)       
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_opm(i,:),'color',colors(i,:))
end
hold off
title('OPM')
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_dipfit-001.jpg']))
close all

% SQUID-GRAD
params.modality = 'squidgrad';
pos_grad = zeros(5,3);
ori_grad = zeros(5,3);
for i = 1:5
    pos_grad(i,:) = squidgrad_dipole{i}.dip.pos;
    [~,idx] = max(vecnorm(squidgrad_dipole{i}.dip.mom,2,1));
    ori_grad(i,:) = squidgrad_dipole{i}.dip.mom(:,idx);
end
mean_pos = mean(pos_grad,1);
tmp = pos_grad;
tmp(:,1) = mean_pos(1)-1;
%100
h=figure;
hold on
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'resolution', 0.1)
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_grad(i,:),'color',colors(i,:))
end
hold off
view(-90,0)
title('SQUID-MAG')
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_dipfit-100.jpg']))
%010
h=figure;
hold on
tmp = pos_grad;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'resolution', 0.1)
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_grad(i,:),'color',colors(i,:))
end
hold off
view(0,0)
title('SQUID-MAG')
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_dipfit-010.jpg']))
%001
h=figure;
hold on
tmp = pos_mag;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'resolution', 0.1)       
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_grad(i,:),'color',colors(i,:))
end
hold off
title('SQUID-MAG')
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_dipfit-001.jpg']))
close all
%% Save 
save(fullfile(save_path, 'dipoles'), 'squidmag_dipole', 'squidgrad_dipole', 'opm_dipole', 'squideeg_dipole', 'opmeeg_dipole'); disp('done');
