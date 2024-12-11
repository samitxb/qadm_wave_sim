% ==============================================================================
% ------------------------------------------------------------------------------
% DESCRIPTION:
% This script simulates the signal processing of a sum of sines and performs
% the following steps:
%
% 1. Generate and Sum Sine Waves:
%    - Create individual sine waves for each frequency in `frequencies`.
%    - Sum the sines into a single signal (`sum_of_sines`).
%
% 2. Normalize and Quantize:
%    - Scale the summed signal to match the DAC resolution.
%    - Quantize the signal to imitate the output of the DAC in the QADM system.
%
% 3. Plot the Initial Results:
%    - Visualize the summed signal in the time domain.
%    - Perform a spectrum analysis (FFT) and plot the frequency domain.
%
% 4. FIR Filter Setup:
%    - Define the normalized cutoff frequency (`W`) based on Nyquist frequency.
%    - Calculate FIR filter coefficients (`B`) and quantize them.
%    - Set the denominator coefficients (`A`) to 1 for the FIR filter model:
%
%                      M
%                     SUM B(k+1) z^(-k)
%                     k=0
%            H(z) = ---------------------
%                        N
%                   1 + SUM A(k+1) z^(-k)
%                       k=1
%
% 5. Filter the Signal:
%    - Apply the FIR filter to the summed signal using the `filter()` function.
%
% 6. Plot Filtered Results:
%    - Visualize the filtered signal in the time domain.
%    - Perform and plot spectrum analysis of the filtered signal.
%
% 7. Multiply by Sine and Cosine:
%    - Multiply the summed signal by sine and cosine components of each
%      frequency.
%    - Quantize the results and filter them again.
%
% 8. Scale and Quantize:
%    - Scale the results to the DAC bit width specified in `waveform_param`.
%
% 9. Plot Final Results:
%    - Visualize the scaled and filtered signals for each frequency.
%    - Show individual sine and cosine components (`A*sin(phi)`, `A*cos(phi)`).
%
% 10. Optional: FIR Frequency Response:
%    - Plot the frequency response of the FIR filter using `freqz()`.
%
% PACKAGES: (pkg load)
% - control:  https://gnu-octave.github.io/packages/control/
% - signal:   https://gnu-octave.github.io/packages/signal/
%
% PARAMETERS:
% - The script relies on settings defined in `waveform_param.m`.
%
% OUTPUTS:
% - Plots for time-domain signals, frequency-domain spectra, and filtered results.
% - Visualization of the FIR filter frequency response (optional).
%
% DATE: 2024-11-27
% GNU Octave, version 8.4.0
% ==============================================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Visual explaination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               sin(wt)
%                         |------(*)-------(TPF) => A*sin(phi)
% sum_of_sines ---[DAC]---
%                         |------(*)-------(TPF) => A*cos(phi)
%                               cos(wt)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
waveform_param;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dac_max = 2^dac_resolution - 1;           % Maximum DAC value
b_max = 2^b_resolution - 1;               % Maximum B value - filter coefficient
nyquist = f_sample/2;

% time axis
t = 0:1/f_sample:duration-(1/f_sample);

% Generate individual sine waves and scale to DAC range
sum_of_sines = zeros(size(t));
for f = frequencies
    sin_wave = sin(2 * pi * f * t);
    sum_of_sines = sum_of_sines + sin_wave;
end

sum_of_sines = quantize(sum_of_sines, dac_resolution);

% plot of the summed sines
figure('Name','Sum of all sines and their FFT',
       'NumberTitle','off');
subplot(2, 1, 1);
plot(t, sum_of_sines);
title('sum of sines');
xlabel('time in sec');
ylabel('amplitude');
xlim([0, duration]);

% spectral analysis
N = length(sum_of_sines);
freq_axis = linspace(0, nyquist, N/2+1); % axis upto nyquist freq (f_sample/2)
Y = fft(sum_of_sines, N);
Y_mag = abs(Y(1:N/2+1));

subplot(2, 1, 2);
stem(freq_axis, Y_mag);
title('spectrum of summed sines');
xlabel('frequency f in Hz');
ylabel('amplitude');
xlim([0, f_max_plot]);

% FIR
W = f_cutoff / nyquist;
B = fir1(N_fir, W, 'low'); % filtercoefficients in B

%Convert B from double to int
B = quantize(B, b_resolution);


A = uint64(1); % FIR always 1 in the denominator (deutsch - "Nenner")

% answer of the filter
[H, W_filter] = freqz(B, A);

% filtered signal after the FIR
filtered_signal = filter(B, A, sum_of_sines);

figure('Name','Filtered sum of all sines and their FFT',
       'NumberTitle','off');
subplot(2, 1, 1);
plot(t, filtered_signal);
title('filtered signal');
xlabel('time in sec');
ylabel('amplitude');
xlim([0, duration]);

% spectral analysis
N = length(filtered_signal);
freq_axis = linspace(0, nyquist, N/2+1); % axis upto nyquist
Y = fft(filtered_signal, N);
Y_mag = abs(Y(1:N/2+1));

subplot(2, 1, 2);
stem(freq_axis, Y_mag);
title('spectrum of filtered signal of summed sines');
xlabel('frequency f in Hz');
ylabel('amplitude');
xlim([0, f_max_plot]);

% multiplication with sin(wt) und cos(wt) of each freq
figure ('Name','sum multiplied by sin(wt) and cos(wt) of each frequency',
        'NumberTitle','off');
freq_count = 1;
for f = frequencies
    % Multiply with sin and cos carriers
    A_sin_phi = sum_of_sines .* quantize(sin(2*pi*f*t), dac_resolution);
    A_cos_phi = sum_of_sines .* quantize(cos(2*pi*f*t), dac_resolution);

    % Filter the scaled results
    A_sin_phi_filtered = filter(B, A, A_sin_phi);
    A_cos_phi_filtered = filter(B, A, A_cos_phi);

    A_sin_phi_filtered = quantize(A_sin_phi_filtered, dac_resolution);
    A_cos_phi_filtered = quantize(A_cos_phi_filtered, dac_resolution);

    % plotting every frequency
    subplot(length(frequencies), 1, freq_count);
    plot(t, A_sin_phi_filtered);
    hold on;
    plot(t, A_cos_phi_filtered);
    hold off;
    title([num2str(f/1e3) ' kHz']);
    legend('A * sin * phi', 'A * cos * phi');
    xlabel('time in sec');
    ylabel('amplitude');
    xlim([0, duration]);
    freq_count = freq_count + 1;
end

% plot of the filter answer
figure('Name','Betrag und Phase mit Filterkoeffizienten B, A(=1).',
       'NumberTitle','off');
freqz(B, A, 4e6);


