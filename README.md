# beta_bursts.m
**Version 3.2 (29/01/2022)**  
**Paul M Briley (brileypm@gmail.com)**  
  
**Citation: Briley PM, Liddle EB, Simmonite M, Jansen M, White TP et al. (2021). Regional Brain Correlates of Beta Bursts in Health and Psychosis: A Concurrent Electroencephalography and Functional Magnetic Resonance Imaging Study. Biological Psychiatry: Cognitive Neuroscience and Neuroimaging, 6, 1145-1156, doi: 10.1016/j.bpsc.2020.10.018**  
  
beta_bursts.m - Matlab function for identifying beta-frequency bursts in a single EEG/MEG time course  
  
returns timings of beta bursts in sample points and seconds, as well as spectral power and peak frequency of each burst. also returns burst duration and spectral width, and can extract power and phase in other frequency bands at the times of beta bursts. can plot data time course, and time-frequency spectrograms, with beta bursts marked. inspired by Shin et al. (2017), eLife 6: e29086 (https://github.com/hs13/BetaEvents)  
    
**usage**  
bursts = beta_bursts(eeg,srate,opt,showFigs)  
  
**requires**  
Matlab image processing toolbox  
mfeeg toolbox by Xiang Wu et al. - http://sourceforge.net/p/mfeeg - for computing Morlet time-frequency spectrograms  
EEGLAB - https://sccn.ucsd.edu/eeglab/ - uses eegplot to display time course (required if showFigs = True)  
  
