% SAK_corr.m (02/07/2017)

close all
clc
clear all
%Params---------------------------------------------------------------
%numSessions= length(sessionFiles);
numCells = 955;         % total number of neurons
numPxl = 512;           % is this necessary??
%trialLength = maxTrialsize;     % max trial array length
timeFrame = 76;         % Each frame is 187.2 ms
direction = 2;          % R-2 and L-3

RAWINDIVTRACE = 1;
TRIALMEAN = 1;
POPMEAN = 1;
SESSMEAN = 1;

sessionFiles = dir('session_data/LD187_*.mat');
sessionNum = 33;
sessionFiles = sessionFiles(sessionNum);        % Session number of interest

load(['session_data/' sessionFiles.name])

% Signal correlation per session
sigCorr = corr(sp');
imagesc(sigCorr);colormap gray;
colorbar
title(sprintf('Signal Correlation: Day %d',sessionNum))
ylabel('Cell ID')
xlabel('Cell ID')

figure
idx = find(tril(ones(size(sigCorr))));
distrCorr = sigCorr(idx);

histogram(distrCorr)
sampleVar = var(distrCorr(distrCorr~=1));
