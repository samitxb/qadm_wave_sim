% - Frequencies: Frequencies of the sine waves to be summed (in Hz).
frequencies = [30e3, 60e3, 90e3, 120e3, 150e5];% Frequenies in Hz

% - Sampling Rate: Sampling frequency of the signal (in Hz).
f_sample = 40e6; % Sample-Frequency in Hz

% - xlimit in frequency domain
f_max_plot = 1.5e6;

% simulation time in sec
duration = 1e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR CONVERTION FROM DOUBLE TO INT
% ------------------------------------------
% done by resolution representing the bits avaliable on hardware and using the
% following methods:

% 1. max_value = 2^resolution-1

%DAC resolution in bits
dac_resolution = 16; % max 64

% See in code:
%--------------
%sum_of_sines = (sum_of_sines + 1) * (dac_max / 2);    % Scale to [0, dac_max]
%dac_signal = round(sum_of_sines);                    % Quantize
%dac_signal = uint64(dac_signal);

%Resolution of coefficients B in bits
b_resolution = 24; % max 64

%See in code
%-----------
%B = B * b_max;
%B = round(B);
%B = uint64(B);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%FIR Parameter
N_fir = 127; % order >> Tabs of B equals N_fir + 1
f_cutoff = 2e3; % Cut-Off-Frequency in Hz
