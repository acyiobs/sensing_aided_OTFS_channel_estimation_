clear;%clc;
rng(2022,'twister');
tic
eta_all = 0.05:0.05:0.3;

%% simulation parameters
sim_param.comm_foldername = './radar_signal_generator/Raytracing_scenarios/6GHz_NLoS_fast/';
sim_param.radar_foldername = './radar_signal_generator/samples/28GHz_NLoS_fast/';
sim_param.filename = '0';
sim_param.guard_tao = 20;
sim_param.support_size = 64;

sim_param.snr_db = 10;
% sim_param.eta = 0.04;   % 0.5;

%% system parameters
sys_param.fft_size = 256;
sys_param.num_symbol = 14;
sys_param.cp_size = sys_param.fft_size * 0.25;

sys_param.num_tx = 32;
sys_param.antenna_interval = 0.5;
sys_param.delta_f=15e3;
sys_param.fc=6e9;
sys_param.c = 299792458;
sys_param.delay_resolution = 1/(sys_param.fft_size*sys_param.delta_f);
sys_param.T = 1 / sys_param.delta_f * (sys_param.fft_size + sys_param.cp_size) / sys_param.fft_size;
sys_param.Doppler_resolution = 1/(sys_param.num_symbol*sys_param.T);


list_of_files = dir([sim_param.radar_foldername, '*.mat']);
num_loop = size(list_of_files, 1);

NMSE_OMP_angle_domain_all = zeros(size(eta_all,2), num_loop);
NMSE_OMP_angle_domain_radar_all = zeros(size(eta_all,2), num_loop);

for eta_idx = 1:size(eta_all, 2)
    sim_param.eta = eta_all(eta_idx);
    for loop = 1:num_loop
        loop
        sim_param.filename = list_of_files(loop).name;
        NMSE = main_simulator(sim_param, sys_param);
        
        NMSE_OMP_angle_domain_all(eta_idx, loop) = NMSE(1);
        NMSE_OMP_angle_domain_radar_all(eta_idx, loop) = NMSE(2);

    end
end

NMSE_OMP_angle_domain_avg = mean(NMSE_OMP_angle_domain_all, 2)
NMSE_OMP_angle_domain_radar_avg = mean(NMSE_OMP_angle_domain_radar_all, 2)

toc

%%
figure;
semilogy(eta_all, NMSE_OMP_angle_domain_avg, 'r-s')
hold on
semilogy(eta_all, NMSE_OMP_angle_domain_radar_avg, 'k-s')

xlabel("eta")
ylabel("NMSE")
legend("OMP (angle domain)", "proposed")
grid on



