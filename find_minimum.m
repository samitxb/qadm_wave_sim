function [min_time, min_index, min_value] = find_minimum(signal, t, t_start, n)

    % Find the starting index based on t_start
    start_index = find(t >= t_start, 1, 'first');

    % Slice the signal and time starting from t_start
    signal_sub = signal(start_index:end);
    t_sub = t(start_index:end);

    % Find all minima indices in the sliced signal
    minima_indices_sub = find(diff(sign(diff(signal_sub))) > 0);

    % Get the nth minimum index in the sliced signal
    min_index_sub = minima_indices_sub(n);

    % Convert back to the original index
    min_index = start_index - 1 + min_index_sub;

    % Get the minimum value
    min_value = signal(min_index);

    % Get the corresponding time
    min_time = t(min_index);
end

