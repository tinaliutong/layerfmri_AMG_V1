% simExptDesign.m
%
%      usage: simExptDesign()
%         by: eli merriam
%       date: 12/01/21
%    purpose: To simulate fMRI time series with different experimental designs
%             This code uses functiosn in the mrTools distribution (https://github.com/justingardner/mrTools)
%             This code may be run with the following commands:
%
%             d = simExptDesign('fixTRs=6');
%             d = simExptDesign('fixTRs=12');
%             d = simExptDesign('fixTRs=24');
%
function d = simExptDesign(varargin)

% check arguments
if ~any(nargin == [0:3])
    help simExptSdesign
    return
end

% evaluate the input arguments
getArgs(varargin, [], 'verbose=0');

if ieNotDefined('numCond'), numCond = 3; end
if ieNotDefined('numBlocks'), numBlocks = 10; end
if ieNotDefined('stimTRs'), stimTRs = 12; end
if ieNotDefined('fixTRs'), fixTRs = 6; end
if ieNotDefined('noiseLevel'), noiseLevel = 0.4; end
if ieNotDefined('respAmp'), respAmp = [1 1.2 1.4]; end
if ieNotDefined('runLength'), runLength = 400; end
if ieNotDefined('doPlots'), doPlots = 1; end

% create design with block-randomized trial order
design = [];
for iBlock = 1:numBlocks
    design = cat(1, design, randperm(numCond));
end
design = design(:);
% generate blocks of trials of stimTRs length
design = repmat(design, 1, stimTRs);
% insert fixation blocks between each stimulus block
design = cat(2, design, zeros(numBlocks*numCond, fixTRs))';
% make a trial order vector
design = design(:);

% get the cannonical HRF with standard parameters
params.x = 6; params.y = 16; params.z = 6; params.stimDur = 0.01; params.incDeriv=0;
% requires mrTools
hrf = hrfDiffGamma(1.5, params);

% create a systemtic time series with respAmp response amplitudes, 
% and a stimulus convolution matrix (scm) used for analysis
scm = zeros(length(design), numCond);
for iCond=1:numCond
    scm(design==iCond, iCond) = 1;
    tSeries(:,iCond,:) = conv(scm(:,iCond), hrf*respAmp(iCond));
end
tSeries = sum(tSeries,2);

% add random noise to the full time series
% this could be extended by adding additional trial-to-trial variability
% (i.e., amplitude noise), and also option for different noise
% distributions (e.g., Rician)
tSeries = tSeries + (randn(size(tSeries)) .* noiseLevel);
%tSeries = 100*(tSeries-mean(tSeries)/mean(tSeries));
tSeries = tSeries/100+1;

% truncate to be the length of the experiment
tSeries = tSeries(1:length(scm));

% ok, now that we have our simulated time series, recover amplitudes
% Two ways to do this. First, we could use a GLM, as shown here:
scm = conv2(scm, hrf);
scm = scm(1:length(tSeries), :);
% truncate
scm = scm(1:runLength,:);
tSeries = tSeries(1:runLength,:);
% and compute betas
betas = pinv(scm)*tSeries;


% Second, could do deconvolution to recover the hemodynamic response
% to do this, we need a stimulus convolution matrix with a 1 for each block onset
newscm = zeros(length(tSeries), numCond);
for iCond = 1:numCond
  temp = [];
  temp(design==iCond, 1) = 1;
  edges = getedges(temp', 0.1);
  newscm(edges.rising', iCond) = 1;
end

% then expand along the diagonals
sc = [];
for i=1:numCond
    sc = cat(2, sc, stimconv(newscm(:,i)',stimTRs+fixTRs));
end

% and get the estimated hemodynamic response using mrTools
tSeries = reshape(tSeries, 1, 1, 1, length(tSeries));
d = getr2timecourse(tSeries, numCond, stimTRs+fixTRs, sc(1:length(tSeries),:), 1.5);

% finally, plot the results!
if doPlots
  smartfig('simtseries', 'reuse');
  clf
  myerrorbar(d.time, d.ehdr(1,:), 'yError', d.ehdrste(1,:), 'MarkerSize=12', 'Color=[0 0 0]', 'MarkerEdgeColor=[1 1 1]');
  hold on
  myerrorbar(d.time, d.ehdr(3,:), 'yError', d.ehdrste(3,:), 'MarkerSize=12', 'Color=[1 0 0]', 'MarkerEdgeColor=[1 1 1]');
  xlabel('Time (s)');
  ylabel('fMRI response (a.u.)');
  xaxis([0 56]);
  yaxis([-2.5 3]);
  drawPublishAxis('xAxisMax', 56, 'yTick', [-2:1:3], 'yAxisMax', 3);
  legend off
end

%% helper functions


% getedges.m
%
%      usage: getedges(timeseries,cutoff)
%         by: justin gardner
%       date: 12/08/03
%       e.g.: getedges(timeseries,4.5)
%    purpose: returns edge fall and rise times
%
function edges = getedges(timeseries,cutoff,dilation)

% check command line arguments
if (nargin == 2)
  dilation = 1;
elseif (nargin ~=3)
  help getedges;
  return
end

% find when timecourse exceeds cutoff value (if cutoff is negative
% assume that means less than cutoff. If cutoff is positive assume
% we are looking for values larger than cutoff).
if (cutoff < 0)
  cutofftimes = timeseries < cutoff;
else
  cutofftimes = timeseries > cutoff;
end

% no events, give up
if (isempty(find(cutofftimes)))
  edges.cutofftimes = cutofftimes;
  edges.rising = [];
  edges.falling = [];
  edges.n = 0;
  edges.dilated = edges.cutofftimes;
  return
end

% find rising edges
rising = [0 find(diff(find(cutofftimes)) > 1)]+1;
% make sure that last one is not problematic
if (rising(length(rising)) > length(cutofftimes))
  rising(length(rising)) = risingedgretimes(length(rising))-1;
end

% find falling edges
falling = [find(diff(find(cutofftimes)) > 1) length(find(cutofftimes))];

% match length
if (length(rising) ~= length(falling))
  falling = falling(1:length(rising));
end

% get times where the signal is over cutoff
findcutofftimes = find(cutofftimes);

% pack return structure
edges.cutofftimes = cutofftimes;
edges.rising = findcutofftimes(rising);
edges.falling = findcutofftimes(falling);

% dilate edges 
dilatedrise = edges.rising-dilation;
dilatedrise(dilatedrise <= 0) = 1;
dilatedfall = edges.falling+dilation;
dilatedfall(dilatedfall > length(timeseries)) = length(timeseries);

% set dilated edges to true
edges.dilated = cutofftimes;
for i = 1:length(dilatedrise);
  edges.dilated(dilatedrise(i):dilatedfall(i)) = 1;
end

edges.n = length(edges.rising);

% get edge length
edges.completeN = min(length(edges.rising),length(edges.falling));
edges.len = edges.falling(1:edges.completeN) - edges.rising(1:edges.completeN);

% hrfDiffGamma.m
%
%      usage: hrfDiffGamma('params'), hrfDiffGamma(tr, params)
%         by: farshad moradi
%       date: 14/06/07
%    purpose: returns a canonical hrf
%
function [hrf] = hrfDiffGamma(tr, params)

if ~any(nargin == [1 2])
  help hrfDiffGamma
  return
end

if tr=='params'
    hrf = {...
        {'description', 'hrfDiffGamma', 'comment describing the hdr model'},...
        {'x', 6, 'hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
        {'y', 16, 'hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
        {'z', 6, 'hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
        {'stimDur', 0.01, 'duration of stimulation/event (seconds, min=0.01s). a boxcar function that is convolved with hrf'},...
        {'incDeriv',0,'type=checkbox','include derivative of the hrf in the model?'},...
    };
    return
end

tmax = max(params.y*3, 20);

if isfield(params, 'tmax')
    tmax = params.tmax;
end

shift = 0;
if isfield(params, 'shift')
    shift = params.shift;
end

dt = 0.01;

t = 0:dt:tmax;
HRF = gampdf(t, params.x, 1) - gampdf(t, params.y, 1)/params.z;
HRF = convn(HRF, ones(1, max(1, ceil(params.stimDur/dt))) );

if shift<0
    HRF = [zeros(1, ceil(-shift/dt)), HRF];
elseif shift>0
    HRF = HRF( ceil(shift/dt):end );
end
    
% subsample hrf
t = [0:length(HRF)-1]*dt;
h_intp = interp1(t, HRF, tr/2:tr:max(t));

% remove mean
h_intp = h_intp - mean(h_intp);
% normalize
h_intp = h_intp / norm(h_intp'); 


if params.incDeriv
    
    % take the derivative
    HRFD = [diff(HRF), 0];
    
    % subsample hrf derivative
    hd_intp = interp1(t, HRFD, tr/2:tr:max(t));
    
    % remove mean
    hd_intp = hd_intp - mean(hd_intp);
    % orthogonalize
    hd_intp = hd_intp - h_intp*(hd_intp/h_intp);
    % normalize
    hd_intp = hd_intp / norm(hd_intp');

    % return as column vectors, zero at time zero
    hrf = [h_intp; hd_intp]';
else
    % return as column vector, zero at time zero
    hrf = h_intp';
end
