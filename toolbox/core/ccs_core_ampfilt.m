function [alff,alff_nor,ts_filt] = ccs_core_ampfilt(ts,TR,f_lp,f_hp,filt_on)
%% CCS_CORE_AMPFILT Compute the amplitude and band-pass filtering 
% of low-frequency fluctuations (LFFs) by using stationary (classic FFT).
%
%   Detailed explanation:
%    INPUT:
%       ts -- original T*N time series (can be 1D or 2D)
%       TR -- sampling segment in time domain
%       f_lp -- low frequency point
%       f_hp -- high freqency point
%       filt_on - the on/off for filt operation
%    OUTPUT:
%       alff -- amplitude of low-frequency fluctuations
%       ts_filt -- timeseries filtered by [f_lp f_hp] band
%       alff_nor -- alff normalized by number of frequency samples
%       
% Credits:
%      Xi-Nian Zuo, PhD of Applied Mathematics
%      Beijing Normal University
%      Email: xinian.zuo@bnu.edu.cn
%      https://zuoxinian.github.io

%% setup for md fft
rate_samp = 1/TR;
if size(ts,1)==1
    ts = ts';
end
if nargin < 5
    filt_on = 'false';
end
nlp = numel(f_lp); nhp = numel(f_hp);
if nlp~=nhp
    disp('The frequency bands for band filtering must be paired.')
else
    num_samp = size(ts,1);
    if rem(num_samp,2)==0 
        freq_idx = 1:(num_samp/2 + 1);
    else
        freq_idx = 1:(num_samp+1)/2;
    end
    freq = rate_samp*(freq_idx-1)/num_samp; %normalized frequencies
    %ts_dm = ts - repmat(mean(ts),num_samp,1);
    ts_dft_full = fft(ts);
    ts_dft = ts_dft_full(freq_idx,:);
    ts_dft = ts_dft/num_samp;
    ts_dft(2:end-1,:) = 2*ts_dft(2:end-1,:);
    amp = abs(ts_dft);
    %estimate ALFF and filter the time series
    idx_lp = find(freq <= f_lp, 1, 'last') + 1;
    idx_hp = find(freq >= f_hp, 1, 'first') - 1;
    amp_filt = amp(idx_lp:idx_hp,:);
    alff = sum(amp_filt); alff_nor = mean(amp_filt);
    %band-pass filtering
    if strcmp(filt_on, 'true')
        ts_filt = bandpass(ts,[f_lp f_hp],rate_samp);
    else
        ts_filt = [];
    end
end

