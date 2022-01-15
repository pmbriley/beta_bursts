function [slide, times] = slide_wins(bursts,events,params,showFig)
%
% Paul M Briley 18/04/2021 (brileypm@gmail.com)
% slide_wins - version 1.0
%
% CITATION
% Briley PM, Liddle EB, Simmonite M, Jansen M, White TP et al. (2020)
% Regional Brain Correlates of Beta Bursts in Health and Psychosis: A Concurrent Electroencephalography and Functional Magnetic Resonance Imaging Study
% Biological Psychiatry: Cognitive Neuroscience and Neuroimaging, 6, 1145-1156
% https://doi.org/10.1016/j.bpsc.2020.10.018
%
% REQUIRES
% burst_rate.m - function for calculating mean beta burst rate in a specific window relative to a set of event markers
%
% USAGE
% [slide, times] = slide_wins(bursts,events,params,showFig)
% Matlab function for calculating mean beta burst rate in sliding time windows relative to a set of event markers (see Fig. 3 in the above paper)
%
% INPUTS
% bursts: vector of times of beta bursts in seconds
% events: vector of times of events in seconds
% params (optional): sliding window parameters...
%       .slideOver: start and end of the time range to be analysed (secs)
%       .slideWindow: width of sliding window (secs)
%       .slideStep: sliding window will move in these steps (secs)
% showFig: if True then will display a plot of burst rate over time
%
% OUTPUTS
% slide: vector giving mean burst rate over time
% times: corresponding time points
%

if nargin<4; showFig = false; end
if nargin<3; params = ''; end
if nargin<2; error('[slide, times] = slide_wins(bursts,events,params)'); end
if ~exist('burst_rate.m','file'); error('requires burst_rate.m'); end
    
if ~isfield(params,'slideOver')
    params.slideOver = [-3 3]; % time range to be analysed (seconds)
    disp('using default for params.slideOver');
end
if ~isfield(params,'slideWindow')
    params.slideWindow = 0.5; % width of sliding window (seconds)
    disp('using default for params.slideWindow');
end
if ~isfield(params,'slideStep')
    params.slideStep = 0.1; % sliding window will move in the steps (seconds)
    disp('using default for params.slideStep');
end

times = params.slideOver(1):params.slideStep:params.slideOver(2);
slide = nan(1,length(times));
for w = 1:length(times)
    t = times(w);
    slide(w) = burst_rate(bursts,events,[t-(params.slideWindow/2) t+(params.slideWindow/2)]);
end

if showFig
    figure;
    plot(times,slide,'k-','linewidth',2); hold on;
    xlabel('Time relative to events (seconds)','FontSize',16);
    ylabel('Mean bursts per second','FontSize',16);
    set(gca,'FontSize',14);
    plot([0 0],get(gca,'YLim'),'k--','linewidth',1);
end
