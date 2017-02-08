function outmat = SAK_load_data(sessionFiles,maxTrialsize);

numSessions= length(sessionFiles);
numCells = 955;		% total number of neurons
numPxl = 512;		% is this necessary??
trialLength = maxTrialsize;	% max trial array length
timeFrame = 76;		% Each frame is 187.2 ms
direction = 2;		% R-2 and L-3

outmat = nan(numSessions,direction,numCells,trialLength,timeFrame);

for i = 1:numSessions
    display(['Loading Session ' num2str(i)]);
    for j = 1:2 %map to 1-R 2-L
    clear rasterMean
    load(['session_data/' sessionFiles(i).name],'rasterActivity');
    a = rasterActivity{j+1};	% cell ID, trials, time
    dim = size(rasterActivity{j});	
    curTrialSize = dim(2);
    A = cat(2,a,zeros(numCells,abs(trialLength-curTrialSize),timeFrame));
    outmat(i,direction,:,:,:) = A;% cellID, trials, time
    %{
    % variable(sessionID, cellID, trialType)
    meanF= squeeze(mean(rasterMean(j+1,:,:),3)); %trial type is in 2,3
    %meanF = meanF.*cellsMatch; % only count the cells with matching ROIs
    expMask = ~isnan(meanF);
    %maxMean = squeeze(max(rasterMean(j+1,:,:),[],3)./meanF);

    F = squeeze(rasterMean(j+1,:,:));
    %F = F.*repmat(cellsMatch,trialLength,1)';

    sessionF(:,:,i,j) = F;
    sessionMeanF(:,i,j) = meanF;
    %sessionExpMask(:,i,j) = expMask;
    %sessionMeanAct(:,i,j) = maxMean;
    %figure; plot(trialAct,'o');
    %}
    end
end
end 
