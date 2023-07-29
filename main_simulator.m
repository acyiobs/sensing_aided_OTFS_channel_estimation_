function NMSE = main_simulator(sim_param, sys_param)
%% Simulation parameters
comm_foldername = sim_param.comm_foldername;  %' 6GHz';
radar_foldername = sim_param.radar_foldername;  %' 6GHz';
filename = sim_param.filename;  %' 546';
snr_db = sim_param.snr_db;
eta = sim_param.eta;   % 0.5;
guard_tao = sim_param.guard_tao;  % 14 -> The delay guard needs to be larger or equal to the maximum delay tap
support_size = sim_param.support_size;
% rng(1,'twister');

M = sys_param.fft_size + sys_param.cp_size;
M_ = sys_param.fft_size;
C = sys_param.cp_size;
N = sys_param.num_symbol;
num_tx = sys_param.num_tx;
z = exp(1j*2*pi/N/M);

%% Generate channel parameters
channel_param = gen_channel_param2(sys_param, comm_foldername, filename);
%% Generate channel
[CH_OTFS_DD2, CH_OTFS_TD] = gen_channel(channel_param, sys_param);

%% Generate data Symbol
x = sign(randn(sys_param.num_tx,sys_param.fft_size,sys_param.num_symbol));

%% Generate pilot Symbol
num_pilot_tao = eta * sys_param.fft_size;
num_pilot_tao = round(num_pilot_tao);
% num_pilot_tao = sys_param.num_tx * ceil(num_pilot_tao/sys_param.num_tx);
sys_param.num_pilot_tao = num_pilot_tao;
num_pilot_v = sys_param.num_symbol;
sys_param.num_pilot_v= num_pilot_v;

% p_Bernoul= sign(randn(sys_param.num_tx,num_pilot_tao,num_pilot_v));

p_Bernoul= sqrt(1/2) * (randn(sys_param.num_tx,num_pilot_tao,num_pilot_v) + 1j*randn(sys_param.num_tx,num_pilot_tao,num_pilot_v));

%% Generate transmit signal
x_p_input = x;  % DD domain data
x_p_input(:, 1:num_pilot_tao, 1:num_pilot_v) = p_Bernoul; % add pilot signal

x_p_input(:, num_pilot_tao+1:num_pilot_tao+guard_tao,:) = 0; % add guard symbol around pilot
x_p_input(:, end-guard_tao+1:end,:) = 0; % add guard symbol around pilot

x_p_input = cat(2, x_p_input(:,end-C+1:end,:), x_p_input);

Stx_symbol = ifft(x_p_input, [], 3) * sqrt(sys_param.num_symbol); % convert to delay-time domain
OTFS_signal = Stx_symbol;

%% Generate noise
Es_No = power(10, snr_db/10);
n = (randn(sys_param.fft_size,sys_param.num_symbol)+1i*randn(sys_param.fft_size,sys_param.num_symbol))*sqrt(1/2)/sqrt(Es_No);

%% Pass the Channel
Rx_signal_OTFS_Bernoul = pass_channel(OTFS_signal, CH_OTFS_TD, sys_param); % delay time domain receive signal
Rx_signal_OTFS_Bernoul = Rx_signal_OTFS_Bernoul(C+1:end,:);

%% Add noise
% Rx_signal_OTFS_Bernoul_no_noise = Rx_signal_OTFS_Bernoul;
Rx_signal_OTFS_Bernoul = Rx_signal_OTFS_Bernoul + n;
Rx_signal_OTFS_Bernoul_dd = fft(Rx_signal_OTFS_Bernoul, [], 2) / sqrt(N);

%% construct the pilot matrix -> this takes 0.17 seconds for (M, N) = (128, 4) debug2 -> this takes 10 seconds for (M, N) = (256, 14) [200 times faster than the initial version]
pilot = zeros(size(x_p_input));
pilot(:, C+1:C+num_pilot_tao, 1:num_pilot_v) = p_Bernoul;

Z = zeros(M*N, N*M);
P = zeros(M*N, num_tx, N*M);

N1 = reshape(0:N-1, N,1,1,1);
M2 = reshape(0:M-1, 1,M,1,1);
N3 = reshape(0:N-1, 1,1,N,1);
M4 = reshape(0:M-1, 1,1,1,M);

i = M*N1 + M2 + 1;
i_p = M*N3 + M4 + 1;

tmp = z .^ (mod(M2-M4, M) .* N3);
tmp = repmat(tmp, N, 1, 1, 1);
Z(i, i_p) = reshape(tmp, N*M, N*M);


tmp = pilot(:, mod(M2-M4, M)+1, mod(N1-N3, N)+1);
tmp = reshape(tmp, num_tx, M, M, N, N);
tmp = permute(tmp, [4 2 1 5 3]);
P(i, :, i_p) = reshape(tmp, [N*M, num_tx, N*M]);

%%
Z = reshape(Z, [size(Z, 1), 1, size(Z, 2)]);
W = reshape(Z.*P, M*N, num_tx*M*N);

%% extract the linear equation corresponding to the pilot
W = reshape(W, M, N, num_tx, M, N);
r_dd = reshape(Rx_signal_OTFS_Bernoul_dd, M_, N);
% remove CP
W = W(C+1:end, :, :, :, :);
% extract pilot
W = W(1:num_pilot_tao, 1:num_pilot_v, :, :, :);
r_dd = r_dd(1:num_pilot_tao, 1:num_pilot_v);
% reshape back
W = reshape(W, num_pilot_tao*num_pilot_v, num_tx*M*N);
r_dd = reshape(r_dd, num_pilot_tao*num_pilot_v, 1);


%% Channel estimation -- OMP with angle domain
beam_domain_pilot = reshape(W, size(W, 1), num_tx, []);
beam_domain_pilot = ifft(beam_domain_pilot, [], 2) * sqrt(num_tx);
beam_domain_pilot = reshape(beam_domain_pilot, size(W));


[h_dd_est_omp_, Omega2] = OMP2(r_dd, beam_domain_pilot, support_size, 0); %Omega = sort(Omega);
h_dd_est_omp_ = ifft(reshape(h_dd_est_omp_, num_tx, []), [], 1) * sqrt(num_tx);

H_DD_est_omp_ = reshape(h_dd_est_omp_, num_tx, M, N);
NMSE_OMP_angle_domain = sum(sum(sum((abs(CH_OTFS_DD2-H_DD_est_omp_)).^2)))/sum(sum(sum((abs(CH_OTFS_DD2)).^2)));

%% Channel estimation -- OMP with angle domain using radar sensing
% calculate support
bf_codebook = dftmtx(num_tx) / sqrt(num_tx);
supp__ = gen_supp_from_radar2(sys_param, bf_codebook, radar_foldername, filename, snr_db);


[h_dd_est_omp_radar, Omega3] = OMP2(r_dd, beam_domain_pilot,  2*size(supp__, 1), supp__); %min(2*size(supp__, 1), support_size)
h_dd_est_omp_radar = ifft(reshape(h_dd_est_omp_radar, num_tx, []), [], 1) * sqrt(num_tx);
H_DD_est_omp_radar = reshape(h_dd_est_omp_radar, num_tx, M, N);
NMSE_OMP_angle_domain_radar = sum(sum(sum((abs(CH_OTFS_DD2-H_DD_est_omp_radar)).^2)))/sum(sum(sum((abs(CH_OTFS_DD2)).^2)));


NMSE = [
    NMSE_OMP_angle_domain
    NMSE_OMP_angle_domain_radar
    ];
end
