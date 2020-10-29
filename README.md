# beta_bursts
Paul M Briley (pmbriley@outlook.com)  
  
beta_bursts.m - Matlab function for identifying beta-frequency bursts in a single EEG/MEG time course   
  
beta_bursts.m returns timings of beta bursts in sample points and in seconds, as well as spectral power and peak frequency of each burst  
it also returns burst duration and spectral width (see opt.propPwr; currently a test feature)  
it plots data time course, and time-frequency spectrograms, with beta bursts marked  
  
**bursts = beta_bursts(eeg,srate,showfigs,opt)**
  
based on work by Shin et al. (2017), eLife 6: e29086 (see also their beta burst identification code available at: https://github.com/hs13/BetaEvents)  

version 1.0 (24/6/2020) - PMB  
published version  

version 1.1 (07/10/2020) - PMB  
change to default value of opt.propPwr to improve identification of burst duration and spectral width  
  
version 1.3 (22/10/2020) - PMB  
outputs now returned as fields in structure 'bursts'  

** requires  

Matlab image processing toolbox  
mfeeg toolbox by Xiang Wu - http://sourceforge.net/p/mfeeg - for computing Morlet time-frequency spectograms (see doi: 10.1038/srep08346)  
Find_Peaks.m by Tony Fast - https://gist.github.com/tonyfast/d7f6212f86ee004a4d2b - for finding peaks in spectrograms using image dilatation method  
EEGLAB - uses eegplot to display time course  

** inputs  

eeg: row vector containing time course   
srate: sample rate in Hz  
showfigs: display time course and spectrograms (true/false)   
opt (options): structure with fields containing analysis parameters...  

opt.m: number of morlet cycles for time-frequency analysis  
opt.f0s: vector of frequencies for morlet analysis  
opt.nMeds: threshold for identifying beta events = median power for a frequency * opt.nMeds  
opt.propPwr: threshold for determining when a beta burst starts/ends (proportion of peak power)  
opt.filt2d: standard deviations for 2D gaussian filter applied to time-frequency spectrograms  
opt.peakFreqs: frequency window to identify peaks in time-frequency spectrogram  
opt.structElem: dimensions of structuring element for image dilatation used in peak identification procedure  
opt.eventGap: minimum gap between beta events in seconds  
opt.dispFreqs: frequency range used for plotting spectrograms and beta events (two elements)  
opt.dispBox: if true, encloses starts and ends, and lower and upper limits of spectral widths, of bursts on spectrograms  
opt.markDur: if true, marks starts and ends of bursts on the scrolling plot of eeg activity  

** outputs  
bursts: structure with fields containing burst properties...  
bursts.tp: locations of beta events in time points  
bursts.secs: locations of beta events in seconds  
bursts.freqs: peak frequency of each beta event in Hz  
bursts.pwr: spectral power of each beta event  
bursts.dur: duration of each beta event in ms (test feature)  
bursts.spec: spectral width of each beta event in Hz (test feature)  
bursts.thresh: threshold power values used at each frequency  
