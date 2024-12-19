function [signal_backshifted, shift] = backshift(signal1, signal2, t, t_start)

    [~, min_index1] = find_minimum(signal1, t, t_start, 1);


    [~, min_index2] = find_minimum(signal2, t, t_start, 1);

    shift = min_index2 - min_index1;

    signal_backshifted = circshift(signal2, -shift);
end
