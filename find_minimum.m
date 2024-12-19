function [min_time, min_index, min_value] = find_smooth_minimum(signal, t, t_start, n)
    smoothing = 50; %smoothing factor

    % to get smooth signal: Moving-Average-Filter
    kernel = ones(smoothing, 1) / smoothing;
    smooth_signal = conv(double(signal), kernel, 'same');

    % starting index based on t_start
    start_index = find(t >= t_start, 1, 'first');

    % prepare subsignals to check inside them
    signal_sub = smooth_signal(start_index:end);
    t_sub = t(start_index:end);

    % find indices of the minima in the signal
    minima_indices_sub = find(signal_sub(1:end-1) < signal_sub(2:end));

    % determine nth minimum inside the subsignal
    min_index_sub = minima_indices_sub(n);

    % get the original index of the full signal
    min_index = start_index - 1 + min_index_sub;

    min_value = signal_sub(min_index_sub);
    min_time = t(min_index);
end

