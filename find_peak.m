function [tmp, peak] = find_peak(data, params, peak)
    tmp = [];
    
    % The indexes of the times are saved
    [~, interval_P300(1)] = min(abs(data.time-params.pretimwin)); % find closest time sample to 270 ms
    [~, interval_P300(2)] = min(abs(data.time-params.posttimwin)); % find closest time sample to 330 ms
    
    % Sensor med max peak i intervallet
    [~, tmp.i_peaktime] = findpeaks(std(data.avg(:, interval_P300(1):interval_P300(2)), [], 1), 'SortStr','descend'); % tmp.i_peaktime är index för tiden då peaken sker.
    if isempty(tmp.i_peaktime)
        tmp.i_peaktime = round((interval_P300(2)-interval_P300(1))/2);
        tmp.nopeak = true;
    end 
    
    tmp.i_peaktime = interval_P300(1)+tmp.i_peaktime(1); % Correct index for peak time. 
    
    tmp.time = data.avg(:,tmp.i_peaktime); % Save the values for all the sensors at the peak time
    
    [tmp.maxval, maxch] = max(tmp.time, [], 1);
    [tmp.minval, minch] = min(tmp.time, [], 1);
    
    % Saving the biggest value and its index
    if abs(tmp.maxval) > abs(tmp.minval)    
        tmp.peakch = data.label{maxch};
        tmp.amplitude = tmp.maxval;
        tmp.i_peakch = maxch;
    else
        tmp.peakch = data.label{minch};
        tmp.amplitude = tmp.minval;
        tmp.i_peakch = minch;
    end

    % Statistics
    peak.labels(end+1,1) = {[params.modality '_MMN data_' params.condition]};
    peak.values(end+1,1) = tmp.amplitude;

end