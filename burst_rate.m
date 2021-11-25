function rate = burst_rate(bursts,events,win)
% Paul M Briley 18/04/2021 (brileypm@gmail.com)
% burst_rate - version 1.0
%
% Citation: Briley PM, Liddle EB, Simmonite M, Jansen M, White TP et al. (2020)
% Regional Brain Correlates of Beta Bursts in Health and Psychosis: A Concurrent Electroencephalography and Functional Magnetic Resonance Imaging Study
% Biological Psychiatry: Cognitive Neuroscience and Neuroimaging. Published online ahead of print. doi: 10.1016/j.bpsc.2020.10.018
%
% rate = burst_rate(bursts,events,win)
%
% Matlab function for calculating mean beta burst rate in a specific window
% relative to a set of event markers
%
% INPUTS
% bursts: vector of times of beta bursts in seconds
% events: vector of times of events in seconds
% win: time window relative to events in seconds (e.g., [0.5 1] would
% compute mean burst rate in the window +0.5 to +1 second after each event)
%
% OUTPUTS
% rate: mean rate in the specific window, in bursts per second
%

if nargin<3; error('rate = burst_rate(bursts,events,win)'); end

nWins = length(events);
burstsPerWin = nan(1,nWins);
for w = 1:nWins
    myWin = [events(w)+win(1) events(w)+win(2)];
    burstsPerWin(w) = sum((bursts>=myWin(1) & bursts<=myWin(2)));
end
rate = mean( burstsPerWin ./ (win(2) - win(1)) );