% ==============================================================================
% ------------------------------------------------------------------------------
%
% DESCRIPTION:
% This script simulates the receive process in a signal processing pipeline. It
% performs the following key steps:
%
% 1. Signal Preparation:
%    - Modulates a `dac_signal` from `waveform_sim` by multiplying it with
%      a sine wave (`sin(wt)`) at its own frequency.
%    - Resamples the modulated signal using the factors `p` and `q` obtained
%      from the `rat()` function. The resampling adjusts the signal to a new
%      sample frequency using `resample()`.
%
% 2. Filtering:
%    - Applies an FIR low-pass filter to the resampled signal. The cutoff
%      frequency is adjusted based on the new sampling rate.
%    - The filter coefficients (`B`) are scaled and quantized before filtering.
%
% 3. Visualization:
%    - Plots the original, resampled, and filtered signals over time.
%    - Performs spectrum analysis on the resampled and filtered signals.
%    - Visualizes the filtered signals multiplied by sine and cosine waves
%      for each carrier frequency.
%
% PACKAGES: (pkg load)
% - control:  https://gnu-octave.github.io/packages/control/
% - signal:   https://gnu-octave.github.io/packages/signal/
%
% PARAMETERS:
% - This script uses configuration settings defined in `recieve_param.m`.
%
% OUTPUTS:
% - Time-domain plots of the modulated, resampled, and filtered signals.
% - Spectrum analysis of the resampled and filtered signals.
% - Per-frequency visualization of filtered signal components (`A * sin(phi)`
%   and `A * cos(phi)`).
%
% DATE: 2024-11-27
% GNU Octave, version 8.4.0
% ==============================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Visual explaination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%    [waveform_sim : dac_signal] --> X ----> OUTPUT
%                                    ^
%                                    |
%                             Störung: sin(wt)
%                                    ^
%                                    |
%              gaussian = A * exp(-((t - u).^2) / (2 * sigma^2));
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waveform_sim;

recieve_param;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nyquist
nyquist = f_sample_recieve/2;

%Gauß Bell Function depending of duration
A = 1;                 % Amplitude
u = duration / 2;     % midpoint
sigma = duration / 10; % Standardabweichung
gaussian = A * exp(-((t - u).^2) / (2 * sigma^2));


% generate the carrier signal
carrier = quantize(sin(2 * pi * carrier_frequency * t) .* gaussian, dac_resolution);



% multiplication with the carrier signal
modulated_signal = sum_of_sines .* carrier;

% resample of the modulated signal
[p, q] = rat(f_sample_recieve / f_sample); %resampling factors

% Resample of the signals modulated and carrier
resampled_carrier = quantize(resample(carrier, p, q), dac_resolution);


resampled_signal = resample(modulated_signal, p, q);


% Time axis for the new sampled signal
t_resampled = (0:length(resampled_signal)-1) / f_sample_recieve;

% FIR
% New W for the FIR bc of the new sampling rate
W = f_cutoff / nyquist;
B = fir1(N_fir, W, 'low'); % filtercoefficients in B
%Convert B from double to int
B = B * b_max;
B = quantize(B, b_resolution);

A = uint64(1); % FIR always 1 in the denominator (deutsch - "Nenner")

% filtered signal
filtered_signal = filter(B, A, resampled_signal);
filtered_signal = quantize(filtered_signal, dac_resolution);



figure('Name','RECIEVE: recieved signal before, after resampling and filtered',
       'NumberTitle','off');
% recieved signal (x sin(wt))
subplot(3, 1, 1);
plot(t, modulated_signal);
title('Signal (before resampling)');
xlabel('Time in sec');
ylabel('Amplitude');
xlim([0, duration]);

% resampled signal
subplot(3, 1, 2);
plot(t_resampled, resampled_signal);
title('Signal (after resampling)');
xlabel('Time in sec');
ylabel('Amplitude');
xlim([0, duration]);

% filtered signal
subplot(3, 1, 3);
plot(t_resampled, filtered_signal);
title('Filtered Signal (after resampling and filtering)');
xlabel('Time in sec');
ylabel('Amplitude');
xlim([0, duration]);

% spectrum analysis
figure('Name','RECIEVE: spectrum analysis of resampled signal and resampled filtered signal',
       'NumberTitle','off');
% spectrum of the resampled
N_resampled = length(resampled_signal);
freq_axis_resampled = linspace(0, nyquist, N_resampled/2+1);
Y_resampled = fft(resampled_signal, N_resampled);
Y_mag_resampled = abs(Y_resampled(1:N_resampled/2+1));

subplot(2, 1, 1);
stem(freq_axis_resampled, Y_mag_resampled);
title('Signal (after resampling)');
xlabel('Frequency in Hz');
ylabel('Amplitude');
xlim([0, nyquist]);

% spectrum of filtered signal
N_filtered = length(filtered_signal);
freq_axis_filtered = linspace(0, nyquist, N_filtered/2+1);
Y_filtered = fft(filtered_signal, N_filtered);
Y_mag_filtered = abs(Y_filtered(1:N_filtered/2+1));

subplot(2, 1, 2);
stem(freq_axis_filtered, Y_mag_filtered);
title('Filtered Signal');
xlabel('Frequency in Hz');
ylabel('Amplitude');
xlim([0, nyquist]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure ('Name','RECIEVE: sum multiplied by sin(wt) and cos(wt) of each frequency',
        'NumberTitle','off');
freq_count = 1;
for f = frequencies
    % Multiply with sin and cos carriers
    A_sin_phi = resampled_signal .* quantize(sin(2*pi*f*t_resampled), dac_resolution);
    A_cos_phi = resampled_signal .* quantize(cos(2*pi*f*t_resampled), dac_resolution);

    % Filter the scaled results
    A_sin_phi_filtered = filter(B, A, A_sin_phi);
    A_cos_phi_filtered = filter(B, A, A_cos_phi);

    A_sin_phi_filtered = quantize(A_sin_phi_filtered, dac_resolution);
    A_cos_phi_filtered = quantize(A_cos_phi_filtered, dac_resolution);

    % plotting every frequency
    subplot(length(frequencies), 1, freq_count);
    plot(t_resampled, A_sin_phi_filtered);
    hold on;
    plot(t_resampled, A_cos_phi_filtered);
    hold off;
    title([num2str(f/1e3) ' kHz']);
    legend('A * sin * phi', 'A * cos * phi');
    xlabel('time in sec');
    ylabel('amplitude');
    xlim([0, duration]);
    freq_count = freq_count + 1;
end


