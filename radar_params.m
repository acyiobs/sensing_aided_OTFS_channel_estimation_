% DeepMIMO parameters utilized in this example
params.radiation_pattern = 0;
params.activate_array_rotation = 0;

% System parameters
params.num_ant_BS = [1, 1, 1];
params.num_ant_UE = [16, 1, 1];
params.ant_spacing_BS = 0.5;
params.ant_spacing_UE = 0.5;
params.scene_frequency = 30; %Hz
params.radar_channel_taps = 650; % Radar channel tap length to generate a time domain radar channel impulse response for one "chirp burst"
params.pulse_shaping = 2;
% params.rolloff_factor = 0.5; %not needed

%Chirp configuration
params.S = 30e12;
params.Fs = 30e6;
params.N_ADC = 512;
params.N_loop = 256;  % Number of chirp bursts (minimum of 1)
params.T_idle = 0e-6;
params.T_start = 0e-6;
params.T_excess = 0e-6;
params.duty_cycle = 1;
% params.F0 = 77e9 - params.S*params.T_start;
params.F0 = 28e9 - params.S*params.T_start;

%Derived configuration
params.T_active = params.N_ADC/params.Fs;
params.T_ramp = params.T_start + params.T_active + params.T_excess;
params.T_chirp = params.T_idle + params.T_ramp;
% params.T_gap = params.T_chirp - params.T_active;
params.T_gap = params.T_idle + params.T_start + params.T_excess;
params.BW_active = params.S*params.T_active;
params.BW_total = params.S*params.T_chirp;
params.fc = params.F0 + params.S*((params.T_active/2)+params.T_start);
params.f_step = params.S/params.Fs;
params.T_PRI = params.radar_channel_taps/params.Fs;  % Pulse repetition interval in seconds, for Doppler processing in FMCW radar
