# beta_bursts.m
**Version 2.32 (24/12/2021)**  
**Paul M Briley (brileypm@gmail.com)**  
  
**Citation: Briley PM, Liddle EB, Simmonite M, Jansen M, White TP et al. (2020). Regional Brain Correlates of Beta Bursts in Health and Psychosis: A Concurrent Electroencephalography and Functional Magnetic Resonance Imaging Study. Biological Psychiatry: Cognitive Neuroscience and Neuroimaging, 6, 1145-1156, doi: 10.1016/j.bpsc.2020.10.018**  
  
beta_bursts.m - Matlab function for identifying beta-frequency bursts in a single EEG/MEG time course  
beta_bursts.m returns timings of beta bursts in sample points and in seconds, as well as spectral power and peak frequency of each burst  
also returns burst duration and spectral width, plots data time course, and time-frequency spectrograms, with beta bursts marked  
can also extract power and phase in other frequency bands at the times of beta bursts  

inspired by Shin et al. (2017), eLife 6: e29086 (see also their beta burst identification code available at: https://github.com/hs13/BetaEvents)  
    
**usage**  
[bursts,tfrOut] = beta_bursts(eeg,srate,opt,showFigs)  
  
**requires**  
Matlab image processing toolbox  
mfeeg toolbox by Xiang Wu et al. - http://sourceforge.net/p/mfeeg - for computing Morlet time-frequency spectrograms  
EEGLAB - uses eegplot to display time course (required if showsfigs = True)  
  
**inputs**  
(only eeg and srate are required)  
eeg: row vector containing time course  
srate: sample rate in Hz  
opt (options): structure with fields containing analysis parameters  
showFigs: display time course and spectrograms (True/False)  
  
opt.m: number of morlet cycles for time-frequency analysis  
opt.f0s: vector of frequencies for morlet analysis  
opt.nMeds: threshold for identifying beta bursts = median power for a frequency * opt.nMeds  
opt.propPwr: threshold for determining when a beta burst starts/ends (proportion of peak power)  
opt.filt2d: standard deviations for 2D gaussian filter applied to time-frequency spectrograms (if empty then no filter applied)  
opt.peakFreqs: frequency window to identify peaks in time-frequency spectrogram  
opt.structElem: dimensions of structuring element for image dilatation used in peak identification procedure  
opt.burstGap: minimum gap between beta bursts in seconds  
opt.useAmp: use time-frequency amplitude instead of power  
opt.useHilbert: use Hilbert transform of bandpass-filtered data for computing power and phase in requested frequency bands (rather than extracting these from the Morlet time-frequency spectrograms)  
opt.dispFreqs: frequency range used for plotting spectrograms and beta events  
opt.dispBox: if True, encloses starts and ends, and lower and upper limits of spectral widths, of bursts on spectrograms  
opt.markDur: if True, marks starts and ends of bursts on the scrolling plot of eeg activity  
opt.bands: frequency bands for measuring power at the times of beta bursts (rows = bands, columns = edges of bands in Hz)  
opt.f0sForOutTFR: can provide a different frequency vector used to compute the optional output time-frequency spectrogram tfrOut.tfr  
opt.verbose: if False, suppresses most command line output (except that produced by functions from other toolboxes)  
   
**outputs**  
(note that only the first output - the bursts structure - is needed for most purposes)  
  
bursts: structure with fields containing burst properties...  
bursts.tp: locations of beta bursts in time points  
bursts.secs: locations of beta bursts in seconds  
bursts.freqs: peak frequency of each beta burst in Hz  
bursts.pwr: spectral power of each beta burst  
bursts.st: start of each burst in seconds  
bursts.ed: end of each burst in seconds  
bursts.stF: start frequency of each burst in Hz  
bursts.edF: end frequency of each burst in Hz  
bursts.dur: duration of each beta burst in ms  
bursts.spec: spectral width of each beta burst in Hz  
bursts.papf: phase at peak frequency at time of burst  
bursts.thresh: threshold power values used at each frequency  
bursts.bandsPower: power in frequency bands specified in opt.bands at times of bursts  
bursts.bandsPhase: phase in frequency bands specified in opt.bands at times of bursts  
bursts.opt: returns the opt (options) structure for reference  
bursts.myver: beta_bursts.m version number  
  
tfrOut.tfr: returns the full, filtered, time-frequency spectrogram (note this can be very large) - this can be of a different frequency resolution to that used in the peak finding algorithm by setting opt.f0sForOutTFR to be a different vector to opt.f0s  
tfrOut.f0s: vector of frequencies for tfrOut.tfr  
tfrOut.times: vector of time points (in seconds) for tfrOut.tfr  
  
note: it is possible to quickly get an opt structure containing the default parameter values - just run bursts = beta_bursts(nan,nan) then look in bursts.opt  
  
