%% test ALFF caculation using different methods
load('tmpts_amyg_lh.mat')
rep_time = 0.72; rate_samp = 1/rep_time;
num_samp = length(tmpts);
time = rep_time * (0:(num_samp-1));

%% plot the raw time series
subplot(411)
plot(time,tmpts,'b')
title('Raw Time Series from Left Amygdala in HCP')
xlabel('Time (sec)')
ylabel('BOLD Sginal (a.u.)')
axis tight; ax = gca; ax.TickDir = 'in'; ax.TickLength = [0.002 0.002];

%% plot the BOLD amplitude
tmpts = tmpts - mean(tmpts);
if rem(num_samp,2)==0 
  freq_idx = 1:(num_samp/2 + 1);
else
  freq_idx = 1:(num_samp+1)/2;
end
freq = rate_samp*(freq_idx-1)/num_samp; %normalized frequencies
ts_dft = fft(tmpts);
ts_dft = ts_dft(freq_idx);
ts_dft = ts_dft/num_samp;
ts_dft(2:end-1) = 2*ts_dft(2:end-1);
amp = abs(ts_dft);
subplot(412)
plot(freq,amp,'r')
title('Amplitude Estimation using FFT')
xlabel('Frequency (Hz)')
ylabel('BOLD Amplitude (a.u.)')
axis tight; ax = gca; ax.TickDir = 'in'; ax.TickLength = [0.002 0.002];

%% plot the BOLD amplitude using padding FFT
padding_length = 2*num_samp;
ts_dft = fft(tmpts,padding_length);
ts_dft = ts_dft(1:num_samp+1);
ts_dft = ts_dft/num_samp;
ts_dft(2:end-1) = 2*ts_dft(2:end-1);
amp = abs(ts_dft(1:2:end));
subplot(413)
plot(freq,amp)
xlabel('Frequency (Hz)')
ylabel('BOLD Amplitude (a.u.)')
axis tight

%% plot the BOLD amplitude using periodogram
[pts,wf1] = periodogram(tmpts,[],num_samp,'power');
amp_psd = sqrt(2*pts);
subplot(413)
plot(freq,amp_psd,'r')
title('Amplitude Estimation using Periodogram')
xlabel('Frequency (Hz)')
ylabel('BOLD Amplitude (a.u.)')
axis tight; ax = gca; ax.TickDir = 'in'; ax.TickLength = [0.002 0.002];

%% plot the BOLD amplitude using Welch method
[pts,wf2] = pwelch(tmpts,[],[],num_samp,'power');
amp_psd = sqrt(2*pts);
subplot(414)
plot(freq,amp_psd,'r')
title('Amplitude Estimation using Welch Method')
xlabel('Frequency (Hz)')
ylabel('BOLD Amplitude (a.u.)')
axis tight; ax = gca; ax.TickDir = 'in'; ax.TickLength = [0.002 0.002];

%% save the figure
set(gcf,'PaperPositionMode','auto')
print('demoALFF','-dpng','-r300')
