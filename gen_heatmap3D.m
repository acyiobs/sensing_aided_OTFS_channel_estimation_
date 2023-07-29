function [heatmap3D, x_angle_mat, R_mat, xx, yy] = gen_heatmap3D(param_file, RX_signal, angle_fft_size)
run(param_file)

%% Angle-Range x-y axes
% Radar KPI
Wavelength = physconst('LightSpeed')/params.fc;
Delta_r = physconst('LightSpeed')/(2*params.BW_active);
Max_r = Delta_r*(params.N_ADC-1);
Delta_v = Wavelength/(2*params.T_PRI*params.N_loop);
Max_v = [-1 ((params.N_loop-2)/params.N_loop)]*(Wavelength/(4*params.T_PRI));

% Range-Angle Axes
antDis = params.ant_spacing_UE;

cos_theta = fliplr(linspace((+0.5/antDis),-(0.5/antDis),(angle_fft_size+1)));
cos_theta = cos_theta(1:angle_fft_size);
sine_theta = sqrt(1-cos_theta.^2);

range_indices = 0:1:(params.N_ADC-1);
[R_mat, cos_theta_mat] = meshgrid((range_indices)*Delta_r,cos_theta);

x_angle_mat = acosd(cos_theta_mat);

%% Range-Doppler axes
[xx, yy] = meshgrid(Max_v(1):Delta_v:Max_v(2), Delta_r*(0:params.N_ADC-1));

%% 3D heatmap
range = fft(RX_signal, size(RX_signal, 1), 1);
%range = range - mean(range, 2);
range_doppler = fft(range, size(RX_signal, 2), 2);
range_doppler_angle = fft(range_doppler, angle_fft_size, 3);

zero_doppler_bins = sort(mod(0, angle_fft_size) + 1);

heatmap3D = fftshift(range_doppler_angle, 2);
heatmap3D = fftshift(heatmap3D, 3);
heatmap3D = abs(heatmap3D).^2;

end