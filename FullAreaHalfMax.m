function m60 = FullAreaHalfMax(sourcedistribution,sourcemodel,params,save_path)
% Calculates the full area at half max amplitude

[~,i1] = min(abs(sourcedistribution.time-0.175));
[~,i2] = min(abs(sourcedistribution.time-0.225));

dat = sourcedistribution.avg.pow(:,i1:i2);
[~,i_latency] = max(mean(abs(dat),1)); % max of mean across sources
i_latency = i1-1+i_latency;
peak_latency = sourcedistribution.time(i_latency);

% Half max level
[peak_pow, i_max] = max(sourcedistribution.avg.pow(:,i_latency));
half_max = peak_pow/2;
peak_loc = sourcemodel.pos(i_max,:);

% Find triangles that have at least one point with amplitude >= half max
i_halfmax_vertices = find(sourcedistribution.avg.pow(:,i_latency)>=half_max);
halfmax_distribution = sourcedistribution.avg.pow(:,i_latency)>=half_max;
[triangles,~] = find(ismember(sourcemodel.tri,i_halfmax_vertices)); 
triangles = sourcemodel.tri(triangles,:);

% Sum area of triangles and divide by 3 (since its a triangle per point).
fahm = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;  

m60 = []; 
m60.peak_latency = peak_latency;
m60.peak_loc = peak_loc;
m60.peak_pow = peak_pow;
m60.fahm = fahm;
m60.halfmax_distribution = halfmax_distribution;

h = figure;
plot(sourcedistribution.time*1e3,sum(abs(sourcedistribution.avg.pow)))
hold on
ylimits = ylim;
latency = 1e3*peak_latency;
plot([latency latency],ylimits,'k--')
hold off
xlabel('t [msec]')
ylabel('total power')
xlim([-params.pre params.post]*1e3);
title(['Summed source power ' params.modality ' - audodd Std'])
saveas(h, fullfile(save_path, 'source analysis', [params.sub '_' params.modality '_mne_sourcepow_audodd-Std.jpg']))

end