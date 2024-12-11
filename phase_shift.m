function [phase_shift, delay_time] = phase_shift(signal1, signal2, fs)
    [xcorr_values, lags] = xcorr(signal1, signal2);

    [~, max_idx] = max(xcorr_values);
    lag = lags(max_idx);

    f_signal = fs / length(signal1);
    T_signal = 1 / f_signal;

    delay_time = lag / fs;

    phase_shift = mod((delay_time / T_signal) * 360, 360);
end

