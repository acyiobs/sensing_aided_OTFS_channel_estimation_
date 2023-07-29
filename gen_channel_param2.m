function channel_param = gen_channel_param2(sys_param, foldername, filename)

load([foldername, filename]);
channel = channels{2}.paths;

fn = fieldnames(channel);
for k=1:numel(fn)
    tmp = channel.(fn{k});
    channel.(fn{k}) = tmp(1:3);
end


M = sys_param.fft_size + sys_param.cp_size;
M_ = sys_param.fft_size;
C = sys_param.cp_size;
N = sys_param.num_symbol;
num_tx = sys_param.num_tx;
delay_resolution = sys_param.delay_resolution;
Doppler_resolution = sys_param.Doppler_resolution;
c = sys_param.c;
fc = sys_param.fc;


channel_param.delay_taps = round((channel.ToA - min(channel.ToA))/ delay_resolution);
channel_param.doppler_taps = mod(round(channel.Doppler_vel / c * fc / Doppler_resolution), N);
channel_param.AoDs = channel.DoD_phi / 180 * pi;
channel_param.pathlosses = sqrt(1e-3*(10.^(0.1*(channel.power)))) .* exp(1j.*(channel.phase/180*pi));
channel_param.pathlosses = channel_param.pathlosses / sqrt(sum(abs(channel_param.pathlosses).^2));


[~,sortIdx] = sort(channel.ToA,'ascend');
channel_param.delay_taps = channel_param.delay_taps(sortIdx);
channel_param.doppler_taps = channel_param.doppler_taps(sortIdx);
channel_param.AoDs = channel_param.AoDs(sortIdx);
channel_param.pathlosses = channel_param.pathlosses(sortIdx);

end