% ============================================================================
% Simulation of recieving a sine like QADM project
% ----------------------------------------------------------------------------
%
% DESCRIPTION:
% This script performs the following steps:
% 1.  Signal Preparation
%     - dac_signal from waveform_sim multiplied by sin(wt) of his own frequency
%     - resample that signal by the factor p, q, which will be output from he
%       function rat() and the different sample frequencies.
%       To resample, the function resample() will  be used.
% 2.  Filter the signal with FIR
%     - adjusted the factor W on the new sample frequency (B, A stay the same
%       as in waveform_sim).
%     - using the function filter() with the coeeficent calculated before
% 3.  Ploting the results
%
% PACKAGES: (pkg load)
% - control:  https://gnu-octave.github.io/packages/control/
% - signal:   https://gnu-octave.github.io/packages/signal/
%
% PARAMETER: recieve_param.m
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
%
%    [waveform_sim : dac_signal] --> X ----> OUTPUT
%                                    ^
%                                    |
%                             Störung: sin(wt)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waveform_sim;

recieve_param;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nyquist
nyquist = f_sample_recieve/2;

% generate the carrier signal
carrier = sin(2 * pi * carrier_frequency * t);

% multiplication with the carrier signal
modulated_signal = sum_of_sines .* carrier;

% resample of the modulated signal
[p, q] = rat(f_sample_recieve / f_sample); %resampling factors

% Resample of the signal
resampled_signal = resample(modulated_signal, p, q);

% Time axis for the new sampled signal
t_resampled = (0:length(resampled_signal)-1) / f_sample_recieve;

% FIR
W = f_cutoff / nyquist;
B = fir1(N_fir, W, 'low'); % filtercoefficients in B
A = 1; % FIR always 1 in the denominator (deutsch - "Nenner")

% filtered signal
filtered_signal = filter(B, A, resampled_signal);



figure('Name','recieved signal before, after resampling and filtered',
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
figure('Name','spectrum analysis of resampled signal and resampled filtered signal',
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



