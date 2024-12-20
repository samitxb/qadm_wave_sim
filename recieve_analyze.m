% ==============================================================================
% receive_analyze.m
% ------------------------------------------------------------------------------
%
% DESCRIPTION:
% This script analyzes a received signal by comparing the modulated signal
% (resampled carrier) and the filtered signal. It calculates the phase shift
% and visualizes the results for multiple frequencies using plots.
%
% PACKAGES: (pkg load)
% - control:  https://gnu-octave.github.io/packages/control/
% - signal:   https://gnu-octave.github.io/packages/signal/
%
% OUTPUTS:
% - Plots visualizing the comparison of signals, including resampled carrier,
%   filtered signal, and their difference.
% - Phase shift in degrees and radians.
%
% DATE: 2024-11-27
% GNU Octave, version 8.4.0
% ==============================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Visual explaination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visual explaination empty ^^; :'(.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recieve_sim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
difference_signal = resampled_carrier - filtered_signal;


figure ('Name','RECIEVE: modulated signal, filtered signal, and signaldifference',
        'NumberTitle','off',
        'Visible', vis_fig8);
% Moduliertes Signal
subplot(3, 1, 1);
plot(t_resampled, resampled_carrier);
title('resampled carrier');
xlabel('time in sec');
ylabel('Amplitude');
xlim([duration/2, duration]);

% Gefiltertes Signal
subplot(3, 1, 2);
plot(t_resampled, filtered_signal);
title('filtered signal');
xlabel('Time in sec');
ylabel('Amplitude');
xlim([duration/2, duration]);

% Unterschiedssignal
subplot(3, 1, 3);
plot(t_resampled, difference_signal);
title('Difference between modulated and filtered');
xlabel('Zeit (s)');
ylabel('Amplitude');
xlim([duration/2, duration]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all;
num_frequencies = length(frequencies);

figure ('Name', 'Receive: Analysis for all Frequencies',
        'NumberTitle', 'off',
        'Visible', vis_fig9);
grid on;

for i = 1:num_frequencies
    f = frequencies(i);

    rows = num_frequencies;
    cols = 4; % 4 Plots / Frequency

    % Plot 1: A * sin(phi) und A * cos(phi)
    subplot(rows, cols, (i - 1) * cols + 1);
    A_sin_phi = resampled_signal .* quantize(sin(2 * pi * f * t_resampled), dac_resolution);
    A_cos_phi = resampled_signal .* quantize(cos(2 * pi * f * t_resampled), dac_resolution);
    plot(t_resampled, A_sin_phi, t_resampled, A_cos_phi);
    legend('A * sin(phi)', 'A * cos(phi)', 'Location', 'northeast');
    title(['Frequency' num2str(f/1e3) ' kHz: Signal with Carrier']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([duration/2, duration]);

    % Plot 2: Resampled carrier

    [min_value, min_index] = min(resampled_carrier);
    min_time = t_resampled(min_index);

    subplot(rows, cols, (i - 1) * cols + 2);
    plot(t_resampled, resampled_carrier);
    hold on;
    line([min_time min_time],
         [min(resampled_carrier) max(resampled_carrier)],
         'Color', 'red',
         'LineStyle', '--'
    );
    plot(min_time, min_value, 'rx', 'MarkerSize', 10, 'LineWidth', 2);

    text(min_time, min_value, sprintf('Min: %.2f', min_value), ...
         'Color', 'blue', 'FontSize', 10);

    hold off;
    title('Resampled Carrier');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([duration/2, duration]);

    % Plot 3: Filtered signal
    A_sin_phi_filtered = quantize(filter(B, A, A_sin_phi), dac_resolution);
    A_cos_phi_filtered = quantize(filter(B, A, A_cos_phi), dac_resolution);

    [min_value, min_index] = min(A_sin_phi_filtered((length(A_sin_phi_filtered)/2*expand_sim_factor):end));
    min_time = t_resampled(min_index + (length(A_sin_phi_filtered)/2*expand_sim_factor));

    subplot(rows, cols, (i - 1) * cols + 3);
    plot(t_resampled, A_sin_phi_filtered, t_resampled, A_cos_phi_filtered);
    hold on;
    line([min_time min_time],
         [min(A_sin_phi_filtered) max(A_sin_phi_filtered)],
         'Color', 'red',
         'LineStyle', '--'
    );

    plot(min_time, min_value, 'rx', 'MarkerSize', 10, 'LineWidth', 2);

    text(min_time, min_value, sprintf('Min: %.2f', min_value), ...
         'Color', 'blue', 'FontSize', 10);
    hold off;
    title('Filtered Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([duration/2, duration]);

    % Plot 4: Difference signal
    difference_signal = double(resampled_carrier) - double(A_sin_phi_filtered);
    subplot(rows, cols, (i - 1) * cols + 4);
    plot(t_resampled, difference_signal);
    title('Difference Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([duration/2, duration]);
    ylim([-32e3, 32e3]);
end

figure ('Name', 'Receive: Analysis for all Frequencies',
        'NumberTitle', 'off',
        'Visible', vis_fig10);

for i = 1:num_frequencies
    f = frequencies(i);

    rows = num_frequencies;
    cols = 4; % 4 Plots / Frequency

    % Plot 1: A * sin(phi) und A * cos(phi)
    subplot(rows, cols, (i - 1) * cols + 1);
    A_sin_phi = resampled_signal .* quantize(sin(2 * pi * f * t_resampled), dac_resolution);
    A_cos_phi = resampled_signal .* quantize(cos(2 * pi * f * t_resampled), dac_resolution);
    plot(t_resampled, A_sin_phi, t_resampled, A_cos_phi);
    legend('A * sin(phi)', 'A * cos(phi)', 'Location', 'northeast');
    title(['Frequency' num2str(f/1e3) ' kHz: Signal with Carrier']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([duration/2, duration]);

    % Plot 2: sqrt(A_sin_phi² + A_cos_phi²)
    sinxcos = sqrt(A_sin_phi.*A_sin_phi .+ A_cos_phi.*A_cos_phi);

    A_sin_phi_filtered = quantize(filter(B, A, A_sin_phi), dac_resolution);
    A_cos_phi_filtered = quantize(filter(B, A, A_cos_phi), dac_resolution);

    sinxcos_filtered = sqrt(A_sin_phi_filtered.*A_sin_phi_filtered .+ ...
                            A_cos_phi_filtered.*A_cos_phi_filtered);

    subplot(rows, cols, (i - 1) * cols + 2);
    plot(t_resampled, sinxcos, t_resampled, sinxcos_filtered);

    title('sqrt(A_sin_phi² + A_cos_phi²) and same with filtered signals');
    legend ('sqrt(sin²+cos²)',
            'FILTERED SIGNALS: sqrt(sin²+cos²)',
            'Location', 'northeast');

    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([duration/2, duration]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure ('Name', 'Receive: Analysis for all Frequencies',
        'NumberTitle', 'off',
        'Visible', vis_fig11);
grid on;

sum_of_differences = zeros(size(t_resampled));
%sum_of_integrals_dif = 0;

for i = 1:num_frequencies
    f = frequencies(i);

    rows = num_frequencies;
    cols = 4; % 4 Plots / Frequency

    % Plot 1: A * sin(phi) und A * cos(phi)
    subplot(rows, cols, (i - 1) * cols + 1);
    A_sin_phi = resampled_signal .* quantize(sin(2 * pi * f * t_resampled), dac_resolution);
    A_cos_phi = resampled_signal .* quantize(cos(2 * pi * f * t_resampled), dac_resolution);
    plot(t_resampled, A_sin_phi, t_resampled, A_cos_phi);
    legend('A * sin(phi)', 'A * cos(phi)', 'Location', 'northeast');
    title(['Frequency' num2str(f/1e3) ' kHz: Signal with Carrier']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([duration/2, duration]);

    % Plot 2: Resampled carrier

    subplot(rows, cols, (i - 1) * cols + 2);
    plot(t_resampled, resampled_carrier);
    title('Resampled Carrier');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([duration/2, duration]);

    % Plot 3: Filtered signal
    A_sin_phi_filtered = quantize(filter(B, A, A_sin_phi), dac_resolution);
    A_cos_phi_filtered = quantize(filter(B, A, A_cos_phi), dac_resolution);

    A_sin_phi_filtered_backshifted = ...
    backshift(resampled_carrier, A_sin_phi_filtered, expand_sim_factor);

    A_cos_phi_filtered_backshifted = ...
    backshift(resampled_carrier, A_cos_phi_filtered, expand_sim_factor);


    subplot(rows, cols, (i - 1) * cols + 3);
    plot(t_resampled, A_sin_phi_filtered_backshifted, t_resampled, resampled_carrier);

    title('Filtered, Backshifted Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([duration/2, duration]);

    % Plot 4: Difference signal
    difference_signal = int64(resampled_carrier) - int64(A_sin_phi_filtered_backshifted);
    sum_of_differences = sum_of_differences .+ abs(difference_signal);

    %sum_of_integrals_dif = sum_of_integrals_dif + trapz(t_resampled, difference_signal);
    subplot(rows, cols, (i - 1) * cols + 4);
    plot(t_resampled, difference_signal);
    title('Difference Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([duration/2, duration]);
end
