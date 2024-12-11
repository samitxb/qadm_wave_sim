function quantized_signal = quantize(signal, res)
    max_value = 2^res-1;
    min_signal = min(signal);
    max_signal = max(signal);
    scale_factor = (max_value - 1) / (max_signal - min_signal);
    signal = (signal - min_signal) * scale_factor;
    quantized_signal = uint64(round(signal));
    quantized_signal = min(max(quantized_signal, 0), max_value - 1);
end


