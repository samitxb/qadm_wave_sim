%störsignal frequenz
carrier_frequency = 1e3;

%abtastrate
f_sample_recieve = 4e6;

% New W for the FIR bc of the new sampling rate
N_fir = 255; % order >> Tabs of B equals N_fir + 1

