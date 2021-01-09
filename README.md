# beta_bursts
**Version 1.8 (09/01/2021)**  
**Paul M Briley (brileypm@gmail.com)**  
  
**Citation: Briley PM, Liddle EB, Simmonite M, Jansen M, White TP et al. (2020). Regional Brain Correlates of Beta Bursts in Health and Psychosis: A Concurrent Electroencephalography and Functional Magnetic Resonance Imaging Study. Biological Psychiatry: Cognitive Neuroscience and Neuroimaging. Published online ahead of print. doi: 10.1016/j.bpsc.2020.10.018**  
  
Based on work by Shin et al. (2017), eLife 6: e29086 (see also their beta burst identification code available at: https://github.com/hs13/BetaEvents)  
  
beta_bursts.m - Matlab function for identifying beta-frequency bursts in a single EEG/MEG time course  
beta_bursts.m returns timings of beta bursts in sample points and in seconds, as well as spectral power and peak frequency of each burst  
also returns burst duration and spectral width, plots data time course, and time-frequency spectrograms, with beta bursts marked  
(peak picking element inspired by code by Tony Fast https://gist.github.com/tonyfast/d7f6212f86ee004a4d2b)  
    
**[bursts,tfrOut] = beta_bursts(eeg,srate,showfigs,opt,out)**  
  
**requires**  
Matlab image processing toolbox  
EEGLAB - uses eegplot to display time course (required if showsfigs = True)  
mfeeg toolbox by Xiang Wu et al. - http://sourceforge.net/p/mfeeg - for computing Morlet time-frequency spectrograms (if bursts.papf or bursts.bandsPhase needed, must also have modified function mf_tfcm2.m)  
  
**inputs**  
eeg: row vector containing time course  
srate: sample rate in Hz  
showfigs: display time course and spectrograms (true/false)  
opt (options): structure with fields containing analysis parameters...  
  
opt.m: number of morlet cycles for time-frequency analysis  
opt.f0s: vector of frequencies for morlet analysis  
opt.nMeds: threshold for identifying beta events = median power for a frequency * opt.nMeds  
opt.propPwr: threshold for determining when a beta burst starts/ends (proportion of peak power)  
opt.filt2d: standard deviations for 2D gaussian filter applied to time-frequency spectrograms (if empty then no filter applied)  
opt.peakFreqs: frequency window to identify peaks in time-frequency spectrogram  
opt.structElem: dimensions of structuring element for image dilatation used in peak identification procedure  
opt.eventGap: minimum gap between beta events in seconds  
opt.dispFreqs: frequency range used for plotting spectrograms and beta events (two elements)  
opt.dispBox: if true, encloses starts and ends, and lower and upper limits of spectral widths, of bursts on spectrograms (default: false)  
opt.markDur: if true, marks starts and ends of bursts on the scrolling plot of eeg activity (default: false)  
opt.bands: frequency bands for measuring power at the times of beta bursts (rows = bands, columns = edges of bands in Hz)  
opt.f0sForOutTFR: can provide a different frequency vector used to compute the optional output time-frequency spectrogram tfrOut.tfr  
  
out: a cell structure containing the output fields you want to compute (to speed up run time by excluding unwanted analyses), can include 'dur', 'spec', 'papf'  
  
**outputs**  
(note that only the first output - the bursts structure - is needed for most purposes)  
  
bursts: structure with fields containing burst properties...  
bursts.tp: locations of beta events in time points  
bursts.secs: locations of beta events in seconds  
bursts.freqs: peak frequency of each beta event in Hz  
bursts.pwr: spectral power of each beta event  
bursts.dur: duration of each beta event in ms  
bursts.spec: spectral width of each beta event in Hz  
bursts.papf: phase at peak frequency at time of burst (test feature, requires modified mfeeg function mf_tfcm2.m)  
bursts.thresh: threshold power values used at each frequency  
bursts.bandsPower: power in frequency bands specified in opt.bands at times of bursts  
bursts.bandsPhase: phase in frequency bands specified in opt.bands at times of bursts (uses midpoint of band) (test feature, requires modified mfeeg function mf_tfcm2.m)  
bursts.opt: returns the opt (options) structure for reference  
  
tfrOut.tfr: returns the full, filtered, time-frequency spectrogram (note this can be very large) - this can be of a different frequency resolution to that used in the peak finding algorithm by setting opt.f0sForOutTFR to be a different vector to opt.f0s  
tfrOut.f0s: vector of frequencies for tfrOut.tfr  
tfrOut.times: vector of time points (in seconds) for tfrOut.tfr  
  
note: it is possible to quickly get an opt structure containing the default parameter values - just run bursts = beta_bursts(nan,nan) then look in bursts.opt  
  
