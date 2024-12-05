% ============================================================================
% Simulation of recieving a sine like QADM project
% ----------------------------------------------------------------------------
% DESCRIPTION:
% This script performs the following steps:
% 1.  Sum the sines up
% 2.  Normalize and scale the analog sum into a digital version to immitate
%     the output from the dac of the QADM
% 3.  Ploting the result so far
% 4.  Preparing the Values of FIR Filter
%     - W: normalized Cutoff Value - fsample / fnyquist
%     - B: N+1 coefficients of the FIR
%     - A: Always 1 by the FIR Model
%             B(z)
%
%                    M
%                   SUM B(k+1) z^(-k)
%                   k=0
%          H(z) = ---------------------
%                      N
%                 1 + SUM A(k+1) z^(-k)
%                     k=1
%
% 5.  Filter by using the coefficients of step 4 and the filter()-function
% 6.  Ploting results again so far
% ==============================================================================
% 7.  Multiply the Sum of sines with each of its own frequencies
% 8.  Scale result to Bitwidth of DAC (declared in waveform_param)
% 9.  Filter the result as shown in the steps before
% 10. Ploting the final results
% ==============================================================================
% 11. Optional Plot: answer of the FIR using freqz()-function
%
% PACKAGES: (pkg load)
% - control:  https://gnu-octave.github.io/packages/control/
% - signal:   https://gnu-octave.github.io/packages/signal/
%
% PARAMETERS: waveform_param.m
%
% OUTPUTS:
%
% AUTHOR:   Sami Taieb,
%           sami.taieb@ic-design.de

% DATE: 2024-11-27
% GNU Octave, version 8.4.0
% ============================================================================

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

% Normalize and scale for DAC
sum_of_sines = sum_of_sines / max(abs(sum_of_sines)); % Normalize to [-1, 1]
sum_of_sines = (sum_of_sines + 1) * (dac_max / 2);    % Scale to [0, dac_max]
dac_signal = round(sum_of_sines);                    % Quantize
dac_signal = uint64(dac_signal);

% plot of the summed sines
figure('Name','Sum of all sines and their FFT',
       'NumberTitle','off');
subplot(2, 1, 1);
plot(t, dac_signal);
title('sum of sines');
xlabel('time in sec');
ylabel('amplitude');
xlim([0, duration]);

% spectral analysis
N = length(dac_signal);
freq_axis = linspace(0, nyquist, N/2+1); % axis upto nyquist freq (f_sample/2)
Y = fft(dac_signal, N);
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
B = B * b_max;
B = round(B);
B = uint64(B);


A = 1; % FIR always 1 in the denominator (deutsch - "Nenner")

% answer of the filter
[H, W_filter] = freqz(B, A);

% filtered signal after the FIR
filtered_signal = filter(B, A, dac_signal);

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
    A_sin_phi = sum_of_sines .* sin(2*pi*f*t);
    A_cos_phi = sum_of_sines .* cos(2*pi*f*t);

    % Scale to 16-bit DAC range
    A_sin_phi_dac = round((A_sin_phi + 1) * (dac_max / 2));
    A_cos_phi_dac = round((A_cos_phi + 1) * (dac_max / 2));

    % Filter the scaled results
    A_sin_phi_filtered = filter(B, A, A_sin_phi_dac);
    A_cos_phi_filtered = filter(B, A, A_cos_phi_dac);

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
freqz(B, A);


