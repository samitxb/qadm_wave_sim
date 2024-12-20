function [signal_backshifted, shift] = backshift(signal1, signal2, expand_sim_factor)

    [~, min_index1] = min(signal1((length(signal1)/2*expand_sim_factor):end));


    [~, min_index2] = min(signal2((length(signal2)/2*expand_sim_factor):end));

    shift = min_index2 - min_index1;

    signal_backshifted = circshift(signal2, -shift);
end
