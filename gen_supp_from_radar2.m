function supp = gen_supp_from_radar2(sys_param, bf_codebook, foldername, filename, snr_db)
M = sys_param.fft_size + sys_param.cp_size;
M_ = sys_param.fft_size;
C = sys_param.cp_size;
N = sys_param.num_symbol;
num_tx = sys_param.num_tx;
delay_resolution = sys_param.delay_resolution;
Doppler_resolution = sys_param.Doppler_resolution;
c = sys_param.c;
fc = sys_param.fc;


load([foldername, filename]);

%%%%%%%%%%% Add Noise %%%%%%%%%%%%%%%
% clustter removal
RX_signal = RX_signal - mean(RX_signal, 2);
% add noise
RX_signal = RX_signal / norm(RX_signal(:));
snr = power(10, snr_db/10);
noise = (randn(size(RX_signal))+1i*randn(size(RX_signal))) *sqrt(1/2) / sqrt(numel(RX_signal)) / sqrt(snr);
RX_signal = RX_signal + noise;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_file = 'radar_params.m';
[heatmap3D, x_angle_mat, R_mat, xx, yy] = gen_heatmap3D(param_file, RX_signal, 128);

static_doppler_bin = size(heatmap3D,2) / 2 + 1;

dynamic_heatmap3D = heatmap3D;
dynamic_heatmap3D(:, static_doppler_bin, :) = 0;

RA_map = squeeze(sum(dynamic_heatmap3D, 2));
RA_peak = peakdetection2d(RA_map, [5,5], [3,3], 10);
%%%%%%%%%%%%%%%%%
RA_peak_filt = colfilt(RA_peak, [5, 5], 'sliding', @max);
RA_peak_filt = (RA_peak==RA_peak_filt) & (RA_peak > (max(RA_peak, [], 2)/2));
%%%%%%%%%%%%%%%%%
RA_peak_idx = find(RA_peak_filt);
[i1, i2] = ind2sub(size(RA_peak), RA_peak_idx);
RA_peak_pwr = RA_peak(RA_peak_idx);

all_peak = [];
for p=1:size(RA_peak_idx, 1)
    Doppler_vec = squeeze(dynamic_heatmap3D(i1(p), :, i2(p)));
    doppler_peak = peakdetection2d(Doppler_vec, [0,10], [0,5], 10);

    doppler_peak_filt = colfilt(doppler_peak, [1, 5], 'sliding', @max);
    doppler_peak_filt = (doppler_peak==doppler_peak_filt) & (doppler_peak > (max(doppler_peak, [], 2)/2));

    doppler_peak_idx = find(doppler_peak_filt);
    doppler_peak_pwr = doppler_peak(doppler_peak_idx);
    for dp=1:size(doppler_peak_idx, 2)
        all_peak = [all_peak; [i1(p), doppler_peak_idx(dp), i2(p), RA_peak_pwr(p), doppler_peak_pwr(dp)]];
    end

end

all_peak = sortrows(all_peak, [4,5], "descend");


range = R_mat(1, all_peak(1:end,1));
Doppler_vel = xx(1, all_peak(1:end,2));
DoD_phi = x_angle_mat(all_peak(1:end,3),2);

range = range - min(range(:));

supp = [];
for p=1:size(range,2)
    AoD = DoD_phi(p) / 180 * pi;
    array_response = exp(-1j * 2 * pi * (0:num_tx-1) * cos(AoD) * sys_param.antenna_interval).';
    angle_response = bf_codebook * array_response;
    [angle_response_sort, angle_response_sort_idx] = sort(angle_response, 'descend');
    [~, max_angle_response_idx] = max(angle_response);
    angle_response_sort_idx = mod(max_angle_response_idx - 1 + (-5:5), num_tx) + 1;% (-5:5)
    for angle_idx = angle_response_sort_idx
        delay = range(p) / c / delay_resolution;

        delay_taps = round(delay);

        doppler = Doppler_vel(p) / c * fc / Doppler_resolution;

        doppler_taps = mod(round(doppler),N);
        for delay_tap=delay_taps
            for doppler_tap=doppler_taps
                idx = num_tx*M*doppler_taps + num_tx*delay_tap + angle_idx;
                supp = [supp, idx];
            end
        end
    end
end
supp = unique(supp','stable');
supp = sort(supp);


end