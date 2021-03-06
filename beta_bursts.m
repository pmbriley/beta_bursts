function [bursts,tfrOut] = beta_bursts(eeg,srate,showfigs,opt)
% Paul M Briley 24/05/2021 (brileypm@gmail.com)
% beta_bursts - version 2.2
%
% Citation: Briley PM, Liddle EB, Simmonite M, Jansen M, White TP et al. (2020)
% Regional Brain Correlates of Beta Bursts in Health and Psychosis: A Concurrent Electroencephalography and Functional Magnetic Resonance Imaging Study
% Biological Psychiatry: Cognitive Neuroscience and Neuroimaging. Published online ahead of print. doi: 10.1016/j.bpsc.2020.10.018
%
% [bursts,tfrOut] = beta_bursts(eeg,srate,showfigs,opt)
% note that only the first output - the bursts structure - is needed for most purposes
%
% Matlab function for identifying beta-frequency bursts/events in single-channel electrophysiological data
% 
% returns timings of beta bursts in sample points and in seconds, as well as spectral power and peak frequency of each burst
% also returns burst duration and spectral width, plots data time course, and time-frequency spectrograms, with beta bursts marked
% (peak picking element inspired by code by Tony Fast https://gist.github.com/tonyfast/d7f6212f86ee004a4d2b)
%
% Based on work by Shin et al. (2017), eLife 6: e29086
% (see also their beta burst identification code available at: https://github.com/hs13/BetaEvents)
%
% version history is at end of file
%
% REQUIRES
% Matlab image processing toolbox
% EEGLAB - uses eegplot to display time course (required if showsfigs = True)
% mfeeg toolbox by Xiang Wu et al. - http://sourceforge.net/p/mfeeg - for computing Morlet time-frequency spectrograms
%
% INPUTS
% (only eeg and srate are required; if not specified, showfigs defaults to False and opt will use default parameters)
% eeg: row vector containing time course
% srate: sample rate in Hz
% showfigs: display time course and spectrograms (True/False)
% opt (options): structure with fields containing analysis parameters...
%
% opt.m: number of morlet cycles for time-frequency analysis
% opt.f0s: vector of frequencies for morlet analysis
% opt.nMeds: threshold for identifying beta events = median power for a frequency * opt.nMeds
% opt.propPwr: threshold for determining when a beta burst starts/ends (proportion of peak power) 
% opt.filt2d: standard deviations for 2D gaussian filter applied to time-frequency spectrograms (if empty then no filter applied)
% opt.peakFreqs: frequency window to identify peaks in time-frequency spectrogram
% opt.structElem: dimensions of structuring element for image dilatation used in peak identification procedure
% opt.eventGap: minimum gap between beta events in seconds
% opt.useAmp: use time-frequency amplitude instead of power (default: False)
% opt.useHilbert: use Hilbert transform of bandpass-filtered data for computing power and phase in requested frequency bands (rather than extracting these from the Morlet time-frequency spectrograms) (default: False)
% opt.dispFreqs: frequency range used for plotting spectrograms and beta events (two elements)
% opt.dispBox: if true, encloses starts and ends, and lower and upper limits of spectral widths, of bursts on spectrograms (default: False)
% opt.markDur: if true, marks starts and ends of bursts on the scrolling plot of eeg activity (default: False)
% opt.bands: frequency bands for measuring power at the times of beta bursts (rows = bands, columns = edges of bands in Hz)
% opt.f0sForOutTFR: can provide a different frequency vector used to compute the optional output time-frequency spectrogram tfrOut.tfr
% opt.verbose: if false, suppresses most command line output (except that produced by functions from other toolboxes) (default: True)
%
% OUTPUTS
% (note that only the first output - the bursts structure - is needed for most purposes)
% bursts: structure with fields containing burst properties...
% bursts.tp: locations of beta events in time points
% bursts.secs: locations of beta events in seconds
% bursts.freqs: peak frequency of each beta event in Hz
% bursts.pwr: spectral power of each beta event
% bursts.dur: duration of each beta event in ms
% bursts.spec: spectral width of each beta event in Hz
% bursts.papf: phase at peak frequency at time of burst
% bursts.thresh: threshold power values used at each frequency
% bursts.bandsPower: power in frequency bands specified in opt.bands at times of bursts
% bursts.bandsPhase: phase in frequency bands specified in opt.bands at times of bursts (uses midpoint of band)
% bursts.opt: returns the opt (options) structure for reference
% bursts.myver: beta_bursts.m version number
%
% tfrOut.tfr: returns the full, filtered, time-frequency spectrogram (note this can be very large) - this can be of a different frequency resolution to that used in the peak finding algorithm by setting opt.f0sForOutTFR to be a different vector to opt.f0s
% tfrOut.f0s: vector of frequencies for tfrOut.tfr
% tfrOut.times: vector of time points (in seconds) for tfrOut.tfr
%
% note: it is possible to quickly get an opt structure containing the default parameter values - just run bursts = beta_bursts(nan,nan) then look in bursts.opt

if nargin<2; error('bursts = beta_bursts(eeg,srate,showfigs,opt,out); only eeg and srate are required arguments; type help beta_bursts for more information'); end
if nargin<3; showfigs = false; end % default is to suppress figures
if nargin<4; opt = []; end % will use default values for all parameters (see below for default values)
if nargout<2; tfrOut.compute = false; else; tfrOut.compute = true; end % if second output not required, don't compute it

if ~ismatrix(eeg) || sum(size(eeg)>1)>1
    error('eeg should be a vector');
end

% prepare optional arguments
args = []; 
allOpts = {'m','f0s','nMeds','propPwr','filt2d','peakFreqs','structElem','eventGap','useAmp','useHilbert','dispFreqs','dispBox','markDur','bands','f0sForOutTFR','verbose'};
if ~isfield(opt,'m');                  opt.m = 5;                            else; args = [args 'm ']; end              % number of morlet cycles for time-frequency analysis
if ~isfield(opt,'f0s');                opt.f0s = 1:1:40;                     else; args = [args 'f0s ']; end            % vector of frequencies for morlet analysis
if ~isfield(opt,'nMeds');              opt.nMeds = 6;                        else; args = [args 'nMeds ']; end          % threshold for identifying beta events = median power for a frequency * opt.nMeds
if ~isfield(opt,'propPwr');            opt.propPwr = 0.5;                    else; args = [args 'propPwr ']; end        % threshold for determining when a beta burst starts/ends (proportion of peak power) 
if ~isfield(opt,'filt2d');             opt.filt2d = [1 3];                   else; args = [args 'filt2d ']; end         % standard deviations for 2D gaussian filter applied to time-frequency spectrograms (if empty then no filter applied)
if ~isfield(opt,'peakFreqs');          opt.peakFreqs = [13 30];              else; args = [args 'peakFreqs ']; end      % frequency window to identify peaks in time-frequency spectrogram
if ~isfield(opt,'structElem');         opt.structElem = [5 5];               else; args = [args 'structElem ']; end     % dimensions of structuring element for image dilatation used in peak identification procedure
if ~isfield(opt,'eventGap');           opt.eventGap = 0.2;                   else; args = [args 'eventGap ']; end       % minimum gap between beta events in seconds
if ~isfield(opt,'useAmp');             opt.useAmp = false;                   else; args = [args 'useAmp ']; end         % if True, will use time-frequency amplitude instead of power (default: false)
if ~isfield(opt,'useHilbert');         opt.useHilbert = false;               else; args = [args 'useHilbert ']; end     % if True, will use Hilbert transformed bandpass-filtered data for calculating phase and amplitude in opt.bands (instead of using the Morlet time-frequency spectrograms)
if ~isfield(opt,'dispFreqs');          opt.dispFreqs = [5 35];               else; args = [args 'dispFreqs ']; end      % frequency range used for plotting spectrograms and beta events (two elements)
if ~isfield(opt,'dispBox');            opt.dispBox = false;                  else; args = [args 'dispBox ']; end        % frequency range used for plotting spectrograms and beta events (two elements)
if ~isfield(opt,'markDur');            opt.markDur = false;                  else; args = [args 'markDur ']; end        % if true, marks starts and ends of bursts on the scrolling plot of eeg activity (default: false)
if ~isfield(opt,'bands');              opt.bands = [];                       else; args = [args 'bands ']; end          % frequency bands for measuring power at the times of beta bursts (rows = bands, columns = edges of bands in Hz)
if ~isfield(opt,'f0sForOutTFR');       opt.f0sForOutTFR = opt.f0s;           else; args = [args 'f0sForOutTFR ']; end   % can provide a different frequency vector used to compute the optional output time-frequency spectrogram 'tfrOut'
if ~isfield(opt,'verbose');            opt.verbose = true;                   else; args = [args 'verbose ']; end        % if false, suppresses most command line output (except that produced by other toolboxes)

bursts.myver = 2.2; % beta_bursts version number
bursts.opt = opt; % store the parameter values in the output bursts structure
if numel(eeg)==1 && isnan(eeg); return; end % this allows you to quickly grab the default parameter values by calling bursts = beta_bursts(nan,nan);

% check optional arguments (throws an error if invalid argument)
isSingInt(opt,{'m','nMeds'});
isVec(opt,{'f0s','f0sForOutTFR'});
isSingNum(opt,{'propPwr','eventGap'});
isTwoInt(opt,{'peakFreqs','structElem','dispFreqs'});
if ~isempty(opt.filt2d); isTwoInt(opt,{'filt2d'}); end
isSingLog(opt,{'dispBox','markDur','useAmp','useHilbert','verbose'});
if ~isempty(opt.bands); is2colMat(opt,{'bands'}); end
if numel(ones(opt.structElem))/2 == round(numel(ones(opt.structElem))/2)
    error('opt.structElem must specify dimensions of a matrix with an ODD number of elements');
end
if length(opt.f0s)==length(opt.f0sForOutTFR) && ~sum(opt.f0s~=opt.f0sForOutTFR)
    opt.f0sForOutTFR = nan;
end
verbose = opt.verbose;
unused = []; flds = fieldnames(opt); 
for a = 1:length(flds); if ~contains(allOpts,flds{a}); unused = [unused flds{a} ' ']; end; end

% check required files and introduce
if verbose; disp(' '); fprintf(1,'** beta_bursts v%s (PMB) **\n',num2str(bursts.myver)); disp('(see code for credits)'); disp(' '); end
if verbose; if isempty(args); disp('all arguments set to defaults (see bursts.opt for values)'); else; fprintf(1,'args accepted: %s\n',args); end; end
if ~isempty(unused); warning('args not recognised: %s',unused); end
if ~exist('mf_tfcm.m','file'); error('requires mfeeg toolbox'); end
if ~exist('imdilate.m','file'); error('requires Matlab image processing toolbox'); end
if (showfigs || (opt.useHilbert && ~isempty(opt.bands))) && ~exist('eegplot.m','file')
    if ~exist('eeglab.m','file'); error('requires EEGLAB to display figures, or compute power/phase in opt.bands at times of bursts when opt.useHilbert is True');
    else
        if verbose; disp(' '); disp('figures requested, or opt.useHilbert is True and opt.bands specified, loading EEGLAB... '); disp(' '); end
        eeglab; disp(' ');
    end
end
if opt.useAmp; tfrType = 'amplitude'; else; tfrType = 'power'; end
dataDur = length(eeg)/srate; dataMins = floor(dataDur/60); dataSecs = mod(dataDur,60);
if verbose; fprintf(1,'threshold: %.0f medians\nburst frequency range: %.0f to %.0f Hz\ninput data: %.0f mins, %.0f secs\n\n',opt.nMeds,opt.peakFreqs(1),opt.peakFreqs(2),dataMins,dataSecs); end

bbTimer = tic;

% compute time-frequency spectrograms using mfeeg toolbox
if verbose; disp('computing time-frequency spectrogram'); end
coeffs = mf_tfcm(eeg,opt.m,opt.f0s,srate,0,0,'coef'); % mfeeg toolbox function
[tfr,phs] = proc_coeffs(coeffs,tfrType);

if tfrOut.compute && ~sum(isnan(opt.f0sForOutTFR))
    if verbose; disp('(separate computation for tfrOut)'); end
    tfrOut.tfr = mf_tfcm(eeg,opt.m,opt.f0sForOutTFR,srate,0,0,'coef');
    tfrOut.tfr = proc_coeffs(tfrOut.tfr,tfrType);
end

% apply 2D gaussian filter
if ~isempty(opt.filt2d)
    if verbose; disp('applying 2D gaussian filter'); end
    tfr = imgaussfilt(tfr,opt.filt2d);
    if tfrOut.compute && ~sum(isnan(opt.f0sForOutTFR)); if verbose; disp('(separate computation for tfrOut)'); end; tfrOut.tfr = imgaussfilt(tfrOut.tfr,opt.filt2d); end
end
meds = median(tfr(:,srate:end),2); % calculate median across time points (ignoring first second) for each frequency

% find peaks in filtered time-frequency spectogram
if verbose; disp('finding peaks in filtered time-frequency spectogram'); end
fInds = find((opt.f0s>=opt.peakFreqs(1)) & (opt.f0s<=opt.peakFreqs(2))); % indices of frequencies for identifying time-frequency peaks
peaks = getPeaks(tfr,opt.structElem); % returns True at locations of peaks
[pksX,pksY] = ind2sub(size(tfr),find(peaks)); % pksX is frequency, pksY is time

% accept peaks that exceed threshold
if verbose; disp('accepting peaks exceeding threshold'); end
thresh = meds .* opt.nMeds;
accept = [];
for i = 1:length(pksY)
    if tfr(pksX(i),pksY(i)) >= thresh(pksX(i)) % if time-frequency power above threshold...
        if (pksX(i)>=fInds(1)) && (pksX(i)<=fInds(end)) % if frequency in range of interest...
            accept = [accept i];
        end
    end
end
pksX = pksX(accept);
pksY = pksY(accept);

% reject peaks that are too close together
if verbose; disp('rejecting peaks that are too close together'); end
keep = true(1,length(pksY));
for i = 1:length(pksY)
    this = pksY(i);
    inds = find(abs(pksY-this)<(opt.eventGap*srate)); % events within minimum gap of index event
    if numel(inds)>1
        pwr = tfr(:,pksY(inds));
        pwrMx = max(pwr,[],1); % maximum power at each event
        [~,toKeep] = max(pwrMx); % keep event with maximum power
        indsBool = false(1,length(inds));
        indsBool(toKeep) = true;
        keep(inds(~indsBool)) = false; % exclude the lower power events within the minimum gap
    end
end
keep(pksY<srate) = false; % exclude events in the first second of the time course
pksX = pksX(keep);
pksY = pksY(keep);

% create outputs
tp = pksY; % times of bursts in time points
secs = tp * (1/srate); % times of bursts in seconds
freqs = opt.f0s(pksX)'; % peak frequencies of bursts
papf = nan(1,length(tp)); % phase at peak frequency
pwr = nan(length(pksX),1); for i = 1:length(pwr); pwr(i) = tfr(pksX(i),pksY(i)); end % power of each burst
for i = 1:length(tp); papf(i) = phs(pksX(i),pksY(i)); end

% find beta event durations
if verbose; disp('finding event durations'); end
st = nan(length(pksX),1); ed = st;
for i = 1:length(st)
    if i==1; m1 = 1; else; m1 = pksY(i-1); end % time point of previous burst
    if i==length(st); p1 = size(tfr,2); else; p1 = pksY(i+1); end % time point of next burst
    prop = pwr(i) * opt.propPwr; % power threshold to determine start/end of burst
    prev = find(tfr(pksX(i),(pksY(i)-1):-1:m1) < prop,1); % time points below threshold before peak
    post = find(tfr(pksX(i),(pksY(i)+1):p1) < prop,1);  % time points below threshold after peak
    if ~isempty(prev); st(i) = pksY(i) - prev(1); end  % burst start sample point
    if ~isempty(post); ed(i) = pksY(i) + post(1); end  % burst end sample point
end
st = st / srate; % convert to seconds
ed = ed / srate;
dur = 1000 * (ed - st); % burst duration in ms

% find beta event spectral widths
if verbose; disp('finding spectral widths'); end
stF = nan(length(pksX),1); edF = stF;
for i = 1:length(stF)
    prop = pwr(i) * opt.propPwr; % power threshold to determine lower and upper frequency limits of burst
    prev = find(tfr((pksX(i)-1):-1:1,pksY(i)) < prop,1); % frequency indices below threshold before peak frequency
    post = find(tfr((pksX(i)+1):end,pksY(i)) < prop,1);  % frequency indices below threshold after peak frequency
    if ~isempty(prev); stF(i) = pksX(i) - prev(1); end % burst lower frequency limit
    if ~isempty(post); edF(i) = pksX(i) + post(1); end % burst upper frequency limit
end
xstF = nan(size(stF)); xedF = nan(size(edF));
xstF(~isnan(stF)) = opt.f0s(stF(~isnan(stF))); % convert to Hz
xedF(~isnan(edF)) = opt.f0s(edF(~isnan(edF)));
stF = xstF; edF = xedF;
spec = edF - stF; % burst spectral width in Hz

% find power in opt.bands at times of beta bursts
if isempty(opt.bands)
    bandsPower = []; bandsPhase = [];
else
    if verbose; disp('finding power and phase in requested frequency bands'); end
    nBands = size(opt.bands,1); % number of frequency bands
    nBursts = length(tp); % number of bursts
    bandsPower = nan(nBursts,nBands); bandsPhase = bandsPower;
    if opt.useHilbert
        for i = 1:nBands
            disp(' '); disp('*** using EEGLAB filter function ***');
            temp = eegfilt(eeg,srate,opt.bands(i,1),opt.bands(i,2)); % uses function from EEGLAB
            disp('************************************');
            temp = hilbert(temp); % hilbert transform of band-pass filtered data
            for ii = 1:nBursts
                bandsPower(ii,i) = abs(temp(pksY(ii)));
                bandsPhase(ii,i) = angle(temp(pksY(ii)));
            end
        end
    else % extract from already-computed Morlet time-frequency spectrograms
        for i = 1:nBands
            fInds = find((opt.f0s>=opt.bands(i,1)) & (opt.f0s<=opt.bands(i,2))); % indices of frequencies for selected band
            index = round(length(fInds)/2);
            for ii = 1:nBursts
                bandsPower(ii,i) = mean(tfr(fInds(index),pksY(ii))); % calculated from the already-derived time-frequency spectrogram
                bandsPhase(ii,i) = phs(fInds(index),pksY(ii)); % calculated from the phase information extracted alongside the time-frequency spectrogram
            end
        end
    end
end

elapsed = toc(bbTimer);
if verbose
    fprintf(1,'\nmean rate %.2f bursts per sec, mean duration %.0f ms\nmean spectral width %.0f Hz, mean peak frequency %.0f Hz\n',length(secs)/dataDur,nanmean(dur),nanmean(spec),nanmean(freqs));
    fprintf(1,'done (processing time: %.0f seconds)\n\n',elapsed);
end

if showfigs
    % create marker structure then display time course with eegplot (from EEGLAB toolbox)
    ind = 0;
    for i = 1:length(tp)
        ind = ind + 1; events(ind).type = 'beta'; events(ind).latency = tp(i);
        if opt.markDur; ind = ind + 1; events(ind).type = 'beta start'; events(ind).latency = st(i) * srate; end
        if opt.markDur; ind = ind + 1; events(ind).type = 'beta end';   events(ind).latency = ed(i) * srate; end        
    end
    eegplot(eeg,'srate',srate','events',events,'color','on');

    % display spectrogram for selected time windows
    done = false;
    while ~done
        disp(' '); disp('use the controls to scroll through the EEG time course, beta bursts are marked with vertical red lines'); disp(' ');
        disp('enter time window for spectrogram display in seconds (e.g., [2 3]), leave blank to exit): ');
        tWin = input('');
        if isempty(tWin)
            done = true;
        else
            if numel(tWin)==1; tWin = [tWin-0.5 tWin+0.5]; end % default time window is one-second long
            tInds = srate*[tWin(1) tWin(2)];
            fInds = [find(opt.f0s==opt.dispFreqs(1)) find(opt.f0s==opt.dispFreqs(2))];

            figure;
            imagesc(tWin,opt.dispFreqs,tfr(fInds(1):fInds(2),tInds(1):tInds(2)));
            hold on;
            set(gca,'YDir','normal');
            set(gca,'FontSize',14);
            xlabel('Time (ms)','FontSize',16);
            ylabel('Frequency (Hz)','FontSize',16);
            plot(secs,freqs,'k+','MarkerSize',10,'linewidth',2);
            if opt.dispBox % display a rectangle to mark the start and end, and lower and upper frequency limits, of each burst
                for i = 1:length(secs)
                    if ~isnan(st(i)) && ~isnan(ed(i)) && ~isnan(stF(i)) && ~isnan(edF(i))
                        rectangle('position',[st(i) stF(i) ed(i)-st(i) edF(i)-stF(i)],'edgecolor','r','linewidth',2);
                    end
                end
            end
        end
    end 
end

bursts.tp = tp; % times of bursts in time points
bursts.secs = secs; % times of bursts in seconds
bursts.freqs = freqs; % peak frequency of each burst in Hz
bursts.pwr = pwr; % spectral power of each burst
bursts.dur = dur; % duration of bursts in milliseconds
bursts.spec = spec; % spectral width of each burst in Hz
bursts.papf = papf'; % phase at peak frequency at time of burst (transposed to match the other outputs)
bursts.thresh = thresh'; % threshold power values used at each frequency (transposed to match opt.f0s)
bursts.bandsPower = bandsPower; % power in frequency bands specified in opt.bands at times of bursts
bursts.bandsPhase = bandsPhase; % phase in frequency bands specified in opt.bands at times of bursts (uses midpoint of band)
if tfrOut.compute
    tfrOut = rmfield(tfrOut,'compute');
    tfrOut.times = 0:(1/srate):(length(eeg)/srate); tfrOut.times = tfrOut.times(1:end-1); % time index for use alongside tfrOut.tfr
    if sum(isnan(opt.f0sForOutTFR)) % in this case, requested output time-frequency spectrogram has same frequency axis as that used for peak finding
        tfrOut.tfr = tfr; % filtered time-frequency spectrogram (note the frequency axis can be different to that used for peak finding by setting opt.f0sForOutTFR to a different frequency vector to opt.f0s)
        tfrOut.f0s = opt.f0s; 
    else % in this case, requested output time-frequency spectrogram has different frequency axis to that used for peak finding
        % tfrOut.tfr has already been computed
        tfrOut.f0s = opt.f0sForOutTFR;
    end
end 

function peaks = getPeaks(tfr,structElem) % get peaks in time-frequency spectrogram using image dilation method
% inspired by Tony Fast (https://gist.github.com/tonyfast/d7f6212f86ee004a4d2b)
f = ones(structElem);
f(ceil(numel(f)./2)) = 0;
peaks = tfr > imdilate(tfr,f);

function [tfr,phs] = proc_coeffs(coeffs,tfrType)
phs = angle(coeffs);
switch tfrType
    case 'amplitude'
        tfr = abs(coeffs);
    case 'power'
        tfr = abs(coeffs).^2;
    case 'db'
        tfr = 10*log10(abs(coeffs).^2);
end

function isSingInt(opt,flds) % check that specified opt fields contain single integers
for i = 1:length(flds)
    if ~isnumeric(opt.(flds{i})) || numel(opt.(flds{i}))~=1 || opt.(flds{i})~=round(opt.(flds{i}))
        error('opt.%s should be a single integer',flds{i});
    end
end

function isVec(opt,flds) % check that specified opt fields contain vectors
for i = 1:length(flds)
    if ~ismatrix(opt.(flds{i})) || sum(size(opt.(flds{i}))>1)>1
        error('opt.%s should be a vector',flds{i});
    end
end

function isSingNum(opt,flds) % check that specified opt fields contain single numbers
for i = 1:length(flds)
    if ~isnumeric(opt.(flds{i})) || numel(opt.(flds{i}))~=1
        error('opt.%s should be a single number',flds{i});
    end
end

function isTwoInt(opt,flds) % check that specified opt fields contain two integers
for i = 1:length(flds)
    if ~isnumeric(opt.(flds{i})) || numel(opt.(flds{i}))~=2 || sum(round(opt.(flds{i}))~=opt.(flds{i}))
        error('opt.%s should contain two integers',flds{i});
    end
end

function isSingLog(opt,flds) % check that specified opt fields contain single logical values
for i = 1:length(flds)
    if numel(opt.(flds{i}))~=1 || ~islogical(opt.(flds{i}))
        error('opt.%s should be a single logical value',flds{i});
    end
end

function is2colMat(opt,flds) % check that specified opt fields contains 2-column matrices
for i = 1:length(flds)
    if ~ismatrix(opt.(flds{i})) || size(opt.(flds{i}),2)~=2
        error('opt.%s should be an nx2 matrix',flds{i});
    end
end

% version 1.0 (24/6/2020) - PMB
% original published version
%
% version 1.1 (07/10/2020) - PMB
% change to default value of opt.propPwr to improve identification of burst duration and spectral width
%
% version 1.2 (08/10/2020) - PMB
% added calculation of power at time of beta bursts for frequency bands specified in opt.bands
%
% version 1.3 (22/10/2020) - PMB
% outputs now returned as fields in structure 'bursts'
%
% version 1.4 (29/10/2020) - PMB
% optionally outputs phase and amplitude in a different specified frequency
% band (or bands) at the times of beta bursts
%
% version 1.5 (05/11/2020) - PMB
% added optional argument 'out' - a cell structure containing the properties of beta
% bursts that you want to compute (to speed up run time by excluding
% unwanted analyses)
%
% version 1.6 (21/12/2020) - PMB
% added citation details, added optional argument checking and reduced dependencies
%
% version 1.7 (26/12/2020) - PMB
% optionally outputs the full, filtered, time-frequency spectrograms plus a time index
%
% version 1.8 (09/01/2021) - PMB
% simplication of output arguments to be the main bursts structure, then an
% optional time-frequency spectrogram
%
% version 1.9 (10/01/2021) - PMB
% automatic creation of modified mfeeg function mf_tfcm2.m if needed 
% (required to obtain phase information from time-frequency spectrograms)
%
% version 2.0 (25/04/2021) - PMB
% option to calculate amplitude and phase in opt.bands at times of beta
% bursts using Hilbert-transformed bandpass filtered data rather than the
% Morlet time-frequency spectrograms
%
% version 2.1 (27/04/2021) - PMB
% removed need for modified mfeeg toolbox function for computing phase
% information
%
% version 2.2 (24/05/2021) - PMB
% sped up calculation of burst duration; removed optional argument 'out' 
% as it no longer leads to performance improvement
%