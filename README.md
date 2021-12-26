# beta_bursts.m
**Version 3.0 (24/12/2021)**  
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
  
