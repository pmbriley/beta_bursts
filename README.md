# beta_bursts  
**Paul M Briley (pmbriley@outlook.com)**  
  
**Citation: Briley PM, Liddle EB, Simmonite M, Jansen M, White TP et al. (2020). Regional Brain Correlates of Beta Bursts in Health and Psychosis: A Concurrent Electroencephalography and Functional Magnetic Resonance Imaging Study. Biological Psychiatry: Cognitive Neuroscience and Neuroimaging. Published online ahead of print. doi: 10.1016/j.bpsc.2020.10.018**  
  
Based on work by Shin et al. (2017), eLife 6: e29086 (see also their beta burst identification code available at: https://github.com/hs13/BetaEvents)  
  
beta_bursts.m - Matlab function for identifying beta-frequency bursts in a single EEG/MEG time course  
beta_bursts.m returns timings of beta bursts in sample points and in seconds, as well as spectral power and peak frequency of each burst  
it also returns burst duration and spectral width (see opt.propPwr)  
it plots data time course, and time-frequency spectrograms, with beta bursts marked  
(peak picking element inspired by code by Tony Fast https://gist.github.com/tonyfast/d7f6212f86ee004a4d2b)  
  
**[bursts,opt] = beta_bursts(eeg,srate,showfigs,opt)**  
  
version 1.0 (24/6/2020) - Paul M Briley (PMB)  
published version  
  
version 1.1 (07/10/2020) - PMB  
change to default value of opt.propPwr to improve identification of burst duration and spectral width  
  
version 1.2 (08/10/2020) - PMB  
added calculation of power at time of beta bursts for frequency bands specified in opt.bands  
  
version 1.3 (22/10/2020) - PMB  
outputs now returned as fields in structure 'bursts'  
  
version 1.4 (29/10/2020) - PMB  
optionally outputs phase and amplitude in a different specified frequency band (or bands) at the times of beta bursts  
  
version 1.5 (05/11/2020) - PMB  
added optional argument 'out' - a cell structure containing the properties of beta bursts that you want to compute (to speed up run time by excluding unwanted analyses)  
  
version 1.6 (21/12/2020) - PMB  
added citation details, added optional argument checking and reduced dependencies  
  
requires  
Matlab image processing toolbox  
EEGLAB - uses eegplot to display time course (required if showsfigs = True)  
mfeeg toolbox by Xiang Wu et al. - http://sourceforge.net/p/mfeeg - for computing Morlet time-frequency spectrograms (if bursts.papf or bursts.bandsPhase needed, must also have modified function mf_tfcm2.m)  
  
inputs  
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
  
out: a cell structure containing the output fields you want to compute (to speed up run time by excluding unwanted analyses), can include 'dur', 'spec', 'papf'  
  
outputs  
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
opt: also returns the opt (options) structure for reference  
