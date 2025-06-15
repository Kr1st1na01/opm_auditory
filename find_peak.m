function [tmp, peak] = find_peak(data, params, peak)
    tmp = [];
    
    % The indexes of the times are saved
    [~, interval_P300(1)] = min(abs(data.time-params.pretimwin)); % find closest time sample to 270 ms
    [~, interval_P300(2)] = min(abs(data.time-params.posttimwin)); % find closest time sample to 330 ms
    
    % Index for the time with sensor that have the overall highest peak during the specified times
    [~, tmp.i_peaktime] = findpeaks(std(data.avg(:, interval_P300(1):interval_P300(2)), [], 1), 'SortStr','descend');
    if isempty(tmp.i_peaktime)
        tmp.i_peaktime = round((interval_P300(2)-interval_P300(1))/2);
        tmp.nopeak = true;
    end 
    
    tmp.i_peaktime = interval_P300(1)+tmp.i_peaktime(1); % Correct index for peak time. 
    
    tmp.time = data.avg(:,tmp.i_peaktime); % Save the peak times for the indexes in tmp.i_peaktime
    
    [tmp.maxval, maxch] = max(tmp.time, [], 1); % The times (values) and indexes that are the highest
    [tmp.minval, minch] = min(tmp.time, [], 1); 
    
    % Saving the biggest value and its index
    if abs(tmp.maxval) > abs(tmp.minval)    
        tmp.peakch = data.label{maxch}; % The name of the peak ch
        tmp.amplitude = tmp.maxval; % The max amplitude
        tmp.i_peakch = maxch; % Index of the ch
        tmp.latency = data.time(:,tmp.i_peaktime);
    else
        tmp.peakch = data.label{minch};
        tmp.amplitude = tmp.minval;
        tmp.i_peakch = minch;
        tmp.latency = data.time(:,tmp.i_peaktime);
    end

    % Statistics (should make this better, not add at the end)
    peak.labels(peak.number_count,1) = {[params.modality '_Evoked data_' params.condition '_amplitude ' params.amp_label]}; % Save amplitude
    peak.values(peak.number_count,1) = tmp.amplitude;

    peak.labels(peak.number_count+1,1) = {[params.modality '_Evoked data_' params.condition '_latency [s]']}; % Save latency
    peak.values(peak.number_count+1,1) = tmp.latency;

    peak.number_count = peak.number_count+2;

end