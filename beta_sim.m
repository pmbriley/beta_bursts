function [eeg,T] = beta_sim()
% basic function for simulating EEG data containing beta bursts
% version 1.0 - Paul M Briley (pmbriley@outlook.com)

srate = 500; % sample rate in Hz
totDur = 120; % total duration in seconds
noiseAmp = 0.5; % amplitude of background noise

events   = [10   20   35   50   63   72   80    90    105];  % times of beta bursts in secs
durs     = [0.2  0.6  0.5  0.4  0.6  0.2  0.1   0.4   0.3];  % duration of each beta burst in secs
freqs    = [13   20   25   20   40   15   8     15    20];   % peak frequency of each burst in Hz
% note that bursts 5 and 7 are NOT beta bursts, occurring at 40 Hz and 8
% Hz, respectively)

eeg = randn(1,srate*totDur)*noiseAmp; % create background noise
%eeg = (rand(1,srate*totDur)-0.5)*noiseAmp; % create background noise
T = 0:(1/srate):(totDur-(1/srate)); % create time vector for output

for i = 1:length(events)
    st = round((events(i) - durs(i)/2)*srate); % sample point corresponding to start of burst
    ed = round((events(i) + durs(i)/2)*srate); % end of burst
    t = 0:(1/srate):durs(i); % time vector for burst
    
    sig1 = sin(2*pi*t*freqs(i)); % burst
    sig2 = sin(2*pi*t*(0.5/durs(i))); % onset-offset ramp
    eeg(st:ed) = eeg(st:ed) + sig1.*sig2; % add to eeg signal
end
