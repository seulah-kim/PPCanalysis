 function sessionF = load_LDMS2_Data(sessionFiles)

numSessions = length(sessionFiles);
numCells = 955;
numPxl = 512;
trialLength = 76;

for i = 1:numSessions
    display(['Loading Session ' num2str(i)]);
    for j = 1:2 %map to 1-R 2-L
    clear rasterMean    
    load(['session_data/' sessionFiles(i).name]);         
    sessionF(:,:,i,j) = squeeze(rasterMean(j+1,:,:));
    
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
