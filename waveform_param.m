% - Frequencies: Frequencies of the sine waves to be summed (in Hz).
frequencies = [50e3, 100e3, 250e3, 500e3, 1e6];% Frequenies in Hz

% - Sampling Rate: Sampling frequency of the signal (in Hz).
f_sample = 40e6; % Sample-Frequency in Hz

% - xlimit in frequency domain
f_max_plot = 1.5e6;

% simulation time in sec
duration = 1e-9;

%DAC resolution in bits
dac_resolution = 24;


%FIR Parameter
N_fir = 255; % order >> Tabs of B equals N_fir + 1
f_cutoff = 2e3; % Cut-Off-Frequency in Hz
