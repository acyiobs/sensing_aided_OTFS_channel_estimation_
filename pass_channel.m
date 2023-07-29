function [RxSig]=pass_channel(Stx,CH_OTFS_TD,sys_param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [RxSig]=PassChannel(TxSig,CH_TD,Nfft,NumOFDMSyms)
%
% INPUTS:      TxSig: transmit signal
%              CH_TD: Time invariant multipath delay Channel
%              Nfft: FFT size for OFDM operation
%              NumOFDMSyms: number of OFDM symbols in each subframe(TTI)
%
% OUTPUT:      RxSig: Received signal after transmitting through the channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = sys_param.fft_size + sys_param.cp_size;
N = sys_param.num_symbol;
num_tx = sys_param.num_tx;

delay_spread = size(CH_OTFS_TD,2) - 1;

stx = reshape(Stx, num_tx, M*N); % convert from time-delay domain to time domain

%% Pass time domain channel
RxSig=zeros(num_tx, M*N);
for ell=0:delay_spread
    mask = ((0:N*M-1) >= ell);
    idx = max((1:N*M) - ell, 1);
    tmp = squeeze(CH_OTFS_TD(:, ell+1, :)) .* stx(:, idx) .*mask;
    RxSig(:, :) = RxSig(:, :) + tmp;
end

RxSig = squeeze(sum(RxSig, 1)); % summation over all transmit antennas
RxSig = reshape(RxSig, M, N); % time-delay domain signal
end
