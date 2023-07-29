function [CH_OTFS_DD, CH_OTFS_TD] = gen_channel(channel_param, sys_param)
% for converting the actual delay (sec) and deoppler frequency (Hz) into delay and doppler tap index --> refer to (2.11) and (2.14) in the book

%% read system params
M = sys_param.fft_size + sys_param.cp_size;
N = sys_param.num_symbol;
d = sys_param.antenna_interval;
num_tx = sys_param.num_tx;

num_path = size(channel_param.pathlosses, 2);

%% discrete-time baseband
z=exp(1j*2*pi/N/M);
delay_spread = channel_param.delay_taps(end);


%% time domain channel
gs=zeros(num_tx, delay_spread+1, N*M);

for p=1:num_path
    pathloss = channel_param.pathlosses(p);
    delay_tap = channel_param.delay_taps(p);
    doppler_tap = channel_param.doppler_taps(p);
    AoD = channel_param.AoDs(p);

    array_response = exp(-1j * 2 * pi * (0:num_tx-1)' * cos(AoD) * d); %% column vector
    tmp = pathloss * z.^( doppler_tap * ((0:N*M-1) - delay_tap) ) .* array_response;
    gs(:, delay_tap+1, :) = gs(:, delay_tap+1, :) + reshape(tmp, [size(tmp,1), 1, size(tmp,2)]);
end

CH_OTFS_TD = gs;

%% delay-doppler domain channel
H_dd_tmp = zeros(num_tx, M, N);
for p=1:num_path
    pathloss = channel_param.pathlosses(p);
    delay_tap = channel_param.delay_taps(p);
    doppler_tap = channel_param.doppler_taps(p);
    AoD = channel_param.AoDs(p);

    array_response = exp(-1j * 2 * pi * (0:num_tx-1)' * cos(AoD) * d);
    tmp = pathloss .* array_response;
    H_dd_tmp(:, delay_tap+1, doppler_tap+1) = H_dd_tmp(:, delay_tap+1, doppler_tap+1) + tmp;
end
CH_OTFS_DD = H_dd_tmp;

end
