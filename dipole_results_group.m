function dipole_results_group(base_save_path, subs, params)

n_subs = max(subs);
n_au = 3;
dist_sqmag_opm = nan(n_subs,1);
dist_sqgrad_opm = nan(n_subs,1);
dist_sqmag_sqgrad = nan(n_subs,1);
spread_opm = nan(n_subs,1);
spread_squidmag = nan(n_subs,1);
spread_squidgrad = nan(n_subs,1);
for i_sub = subs
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);
    clear squidmag_dipole squidgrad_dipole opm_dipole
    load(fullfile(save_path, 'source analysis', 'dipoles')); 
    dipole_squidmag{i_sub} = squidmag_dipole;
    dipole_squidgrad{i_sub} = squidgrad_dipole;
    dipole_opm{i_sub} = opm_dipole;
  
    % Metrics: 
    % - distance between dipoles for same phalange different systems
    % - over phalanges: average distance from mean location within distance
    pos_squidmag = zeros(n_subs,3);
    pos_squidgrad = zeros(n_subs,3);
    pos_opm = zeros(n_subs,3);
    
    pos_opm(i_sub,:) = dipole_opm{i_sub}.dip.pos(1,:);

    for s = 1:2
        pos_squidmag(i_sub,:) = dipole_squidmag{i_sub}.dip.pos(s,:);
        pos_squidgrad(i_sub,:) = dipole_squidgrad{i_sub}.dip.pos(s,:);
        tmp.mag_opm(s) = 1e1*norm(pos_squidmag-pos_opm);
        tmp.grad_opm(s) = 1e1*norm(pos_squidgrad-pos_opm);
    end
    if abs(tmp.mag_opm(1)) > abs(tmp.mag_opm(2))
        pos_squidmag(i_sub,:) = dipole_squidmag{i_sub}.dip.pos(2,:);
    else
        pos_squidmag(i_sub,:) = dipole_squidmag{i_sub}.dip.pos(1,:);
    end
    if abs(tmp.grad_opm(1)) > abs(tmp.grad_opm(2))
        pos_squidgrad(i_sub,:) = dipole_squidgrad{i_sub}.dip.pos(2,:);
    else
        pos_squidgrad(i_sub,:) = dipole_squidgrad{i_sub}.dip.pos(1,:);
    end

    dist_sqmag_opm(i_sub) = 1e1*norm(pos_squidmag-pos_opm);
    dist_sqgrad_opm(i_sub) = 1e1*norm(pos_squidgrad-pos_opm);
    dist_sqmag_sqgrad(i_sub) = 1e1*norm(pos_squidmag-pos_squidgrad);
    
    spread_opm(i_sub,:) = mean(1e1*vecnorm(pos_opm-mean(pos_opm,1),2,2))'; % mean distance from center
    spread_squidmag(i_sub,:) = mean(1e1*vecnorm(pos_squidmag-mean(pos_squidmag,1),2,2))'; % mean distance from center 
    spread_squidgrad(i_sub,:) = mean(1e1*vecnorm(pos_squidgrad-mean(pos_squidgrad,1),2,2))'; % mean distance from center
end

%% Plot distances
h = figure('DefaultAxesFontSize',16);
bar(1,mean(dist_sqmag_opm,1,'omitnan'));
hold on
er = errorbar(1,mean(dist_sqmag_opm,1,'omitnan'), mean(dist_sqmag_opm,1,'omitnan')-min(dist_sqmag_opm,[],1,'omitnan'), mean(dist_sqmag_opm,1,'omitnan')-max(dist_sqmag_opm,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQUID-MAG to OPM (mean = ' num2str(mean(mean(dist_sqmag_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
%xticklabels('test', 'test2', 'test3')
saveas(h, fullfile(base_save_path, 'statistics', 'dipole_squidmag_to_opm_dist.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1,mean(dist_sqgrad_opm,1,'omitnan'));
hold on
er = errorbar(1,mean(dist_sqgrad_opm,1,'omitnan'), mean(dist_sqgrad_opm,1,'omitnan')-min(dist_sqgrad_opm,[],1,'omitnan'), mean(dist_sqgrad_opm,1,'omitnan')-max(dist_sqgrad_opm,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQUID-GRAD to OPM (mean = ' num2str(mean(mean(dist_sqgrad_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
%xticklabels('test', 'test2', 'test3')
saveas(h, fullfile(base_save_path, 'statistics', 'dipole_squidgrad_to_opm_dist.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1,mean(dist_sqmag_sqgrad,1,'omitnan'));
hold on
er = errorbar(1,mean(dist_sqmag_sqgrad,1,'omitnan'), mean(dist_sqmag_sqgrad,1,'omitnan')-min(dist_sqmag_sqgrad,[],1,'omitnan'), mean(dist_sqmag_sqgrad,1,'omitnan')-max(dist_sqmag_sqgrad,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQUID-MAG to SQUID-GRAD (mean = ' num2str(mean(mean(dist_sqmag_sqgrad,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
%xticklabels('test', 'test2', 'test3')
saveas(h, fullfile(base_save_path, 'statistics', 'dipole_squidmag_to_squidgrad_dist.jpg'))
close

% for i_ph = 1:5
%     h = figure('DefaultAxesFontSize',16);
%     plot(subs,dist_sqmag_opm(subs,i_ph),'+-');
%     hold on
%     plot(subs,dist_sqgrad_opm(subs,i_ph),'x-');
%     plot(subs,dist_sqmag_sqgrad(subs,i_ph),'*-');
%     hold off
%     title([params.phalange_labels{i_ph} ' - dipole distances over subjects'])
%     ylabel('Distance [mm]')
%     xlabel('Subjects')
%     legend(['SQMAG-OPM   '; 'SQGRAD-OPM  '; 'SQMAG-SQGRAD'])
%     saveas(h, fullfile(base_save_path, 'figs', ['dipole_dist_vs_sub-' params.phalange_labels{i_ph} '.jpg']))
%     close
% end

%% Plot spread squidmag vs opm
data1 = spread_squidmag;
data2 = spread_opm;
data3 = spread_squidgrad;
mean1 = mean(data1,1,'omitnan');
mean2 = mean(data2,1,'omitnan');
mean3 = mean(data3,1,'omitnan');
min1 = min(data1,[],1,'omitnan');
min2 = min(data2,[],1,'omitnan');
min3 = min(data3,[],1,'omitnan');
max1 = max(data1,[],1,'omitnan');
max2 = max(data2,[],1,'omitnan');
max3 = max(data3,[],1,'omitnan');
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];
err3 = [mean3-min3; max3-mean3];

h = figure('DefaultAxesFontSize',16);
h.Position(3) = round(h.Position(3)*1.3);
bar(1,[mean1; mean2; mean3]','grouped');
hold on
errorbar(1-0.22,mean1(1),err1(1,1),err1(2,1),'k','linestyle','none');
errorbar(1,mean2(1),err2(1,1),err2(2,1),'k','linestyle','none');
errorbar(1+0.22,mean3(1),err3(1,1),err3(2,1),'k','linestyle','none');
p_values = zeros(1, 3);
[~, p_values(1)] = ttest(data1, data2);
[~, p_values(2)] = ttest(data2, data3);
[~, p_values(3)] = ttest(data1, data3);
sigstar({[1-0.22, 1]}, p_values(1));
sigstar({[1, 1+0.22]}, p_values(2));
sigstar({[1-0.22, 1+0.22]}, p_values(3));
hold off
title('Group level dipole spread')
ylabel('Dipoles spread [mm]')
legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
saveas(h, fullfile(base_save_path, 'statistics', 'dipole_spread.jpg'))

%% Plot spread squidgrad vs opm
data1 = spread_squidgrad;
data2 = spread_opm;
mean1 = mean(data1,1,'omitnan');
mean2 = mean(data2,1,'omitnan');
min1 = min(data1,[],1,'omitnan');
min2 = min(data2,[],1,'omitnan');
max1 = max(data1,[],1,'omitnan');
max2 = max(data2,[],1,'omitnan');
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];

h = figure('DefaultAxesFontSize',16);
bar(1,[mean1; mean2]','grouped');
hold on
k=1;
errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
[~, p_values] = ttest(data1(:, 1), data2(:, 1));

sigstar({[1, 1]}, p_values);
hold off
title('Group level M100 dipole spread')
ylabel('Dipoles spread [mm]')
%xlabel('Phalange')
legend({'squidgrad','opm'});
%xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'statistics', 'dipole_spread_squidgrad_opm.jpg'))

%% Plot spread squidgrad vs squidmag
data1 = spread_squidgrad;
data2 = spread_squidmag;
mean1 = mean(data1,1,'omitnan');
mean2 = mean(data2,1,'omitnan');
min1 = min(data1,[],1,'omitnan');
min2 = min(data2,[],1,'omitnan');
max1 = max(data1,[],1,'omitnan');
max2 = max(data2,[],1,'omitnan');
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];

h = figure('DefaultAxesFontSize',16);
bar(1,[mean1; mean2]','grouped');
hold on
k=1;
errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
[~, p_values] = ttest(data1(:, 1), data2(:, 1));

sigstar({[1, 1]}, p_values);
hold off
title('Group level M100 dipole spread')
ylabel('Dipoles spread [mm]')
%xlabel('Phalange')
legend({'squidgrad','squidmag'});
%xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'statistics', 'dipole_spread_squidgrad_squidmag.jpg'))

%% Boxplot
h = figure;
boxplot([dist_sqmag_opm, dist_sqgrad_opm, dist_sqmag_sqgrad], {'MAG vs OPM', 'GRAD vs OPM', 'MAG vs GRAD'})
title('Dipole fit distances - Std tone, M100')
ylabel('Distance [mm]')
saveas(h, fullfile(base_save_path, 'statistics', 'boxplot_comparing_dipole_distances.jpg'))

close all
end