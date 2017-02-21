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

cellID=ParseCells(sessionNum);

sp = sp(cellID,:);
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

sample = distrCorr(distrCorr~=1);
sample_idx = find(distrCorr~=1);	% key back to the original list of upper triangle
figure
%sampleVar = nanvar(sample);
%sampleMean = nanmean(sample);
%H=ztest(sample,sampleMean,sampleVar);

% Plot the histograms
hold on
g=histogram(sample);	% using number of bins = sqrt(number of elements)
g(1).FaceAlpha=0.2;

% Shuffle time independently from cell-to-cell in the original data to generate null distribution (100 times)
[m,n]=size(sp);
clearvars null
for simNum = 1:100
clearvars spNull nullCorr null_temp
for c = 1:m
    spNull(c,:) = sp(c,randperm(n));
end
nullCorr = corr(spNull');
null_temp = nullCorr(find(tril(ones(size(nullCorr)))));
[nullValues, nullEdges]= histcounts(null_temp(null_temp~=1),g.BinEdges);
null(simNum,:) = nullValues;	% simulation number x histogram bins
end

%null = sampleMean+sqrt(sampleVar).*randn(length(sample),1);	% not sure if this is valid thing to do. The area under null distribution does not match that of sample distribution. 

% Compute values of 5% and 95% of the distribution of each bin

nullBounds=prctile(null,[5 95]);

xcoordEdges=repelem(g.BinEdges,2);
xcoordEdges=xcoordEdges(2:length(xcoordEdges)-1);
ycoordEdges=[repelem(nullBounds(1,:),2);repelem(nullBounds(2,:),2)];
h=fill([xcoordEdges, fliplr(xcoordEdges)],[ycoordEdges(2,:) , fliplr(ycoordEdges(1,:))],'k');    % Coloring the significant area between 5 and 95 percentile
h.FaceAlpha = 0.2; h.EdgeColor = 'none'; h.EdgeAlpha = 0;

legend('sample','null')

title(sprintf('Sig-Corr Histogram: Session %d',sessionNum))
%h=histogram(null);
%h(1).FaceAlpha= 0.2;
%g=histfit(sample,'BinEdges',h.BinEdges);		
%sampleZ = zscore(sample);
%outlier_idx1 = find(sampleZ>1.96);
%[sortZ,sort_idx] = sort(sampleZ(outlier_idx1));
%outlier_idx2 = idx(sample_idx(outlier_idx1));

%[sampleZ,sortZ_idx] = sort(sample);

%for jj = 1:length(g.Values)		% Indices in the histogram format in the sorting order (smallest to largest)
%    if jj == 1 
%    strIdx = 1;
%    endIdx = g.Values(jj);
%    elseif jj == length(g.Values)
%    strIdx = sumValues +1;
%    endIdx = length(g.Values);
%    else 
%    strIdx = sumValues+1;
%    endIdx = sumValues+g.Values(jj)-1;
%    end
%    
%    M{jj} = sortZ_idx(strIdx:endIdx);
%    sumValues = endIdx;
%end

%outlier_idx1 = find(sampleZ>1.96);
%[sortZ,sort_idx] = sort(sampleZ(outlier_idx1));
%outlier_idx2 = idx(sample_idx(outlier_idx1));
thresCorr = 0.035; 	% Just from eyeballing the histogram. Need to come up with a justifiable method for setting a threshold value
outlier_idx1 = find(sample>thresCorr);	
[sortV,sort_idx] = sort(sample(outlier_idx1));	% Sort by sig-corr strength
outlier_idx2 = idx(sample_idx(outlier_idx1));

NeuronPair= zeros(size(sigCorr));
NeuronPair(outlier_idx2)=1;		% Location map
NeuronPairStrength=zeros(size(sigCorr));
NeuronPairStrength(outlier_idx2(sort_idx))=sortV;
%NeuronPairStrength = NeuronPairStrength + diag(diag(ones(size(sigCorr))));	% add 1 for diagonal entries
figure
imagesc(1-NeuronPair);colormap gray; hold on; line([1 78],[1 78],'Color','b')
title(sprintf('Connectivity: Session %d, threshold %d',sessionNum,round(thresCorr,3)))
figure
imagesc(NeuronPairStrength);colorbar; hold on; line([1 78],[1 78],'Color','w')
title(sprintf('Connectivity strength: Session %d, threshold %d',sessionNum,round(thresCorr,3)))

