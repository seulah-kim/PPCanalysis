% AK 20160923 
% Load and process Laura's 2P data for LD187

cd('~/Documents/MATLAB/PPC/LD187');
sessionFiles = dir('session_data/LD187_*.mat');
numSessions = length(sessionFiles);
numCells = 955;
numPxl = 512;
trialLength = 76;

load('labels_all.mat');
load('wei_masks.mat');
%% Parce Cell Matching Scores

% labels_all, 1 is good, 2 is ok, 3 is not match
match1 = labels_all == 1;
match2 = labels_all == 2;
matchCells = match1 | match2;

cells1 = sum(match1,1)>0; numCells1 = sum(cells1);
cells2 = sum(match2,1)>0; numCells2 = sum(cells2);
cellsMatch = sum(matchCells,1)>0; numMatchCells = sum(cellsMatch);

%%

%% Load Data
sessionMeanF = NaN(numCells,numSessions,2);
sessionExpMask = NaN(numCells,numSessions,2);
sessionMeanAct = NaN(numCells,numSessions,2);
sessionF = NaN(numCells,trialLength,numSessions,2);

for i = 1:numSessions
    tic
    for j = 1:2 %map to 1-R 2-L
    clear rasterMean
    load(['session_data/' sessionFiles(i).name]);         
   
    % variable(sessionID, cellID, trialType)
    meanF= squeeze(mean(rasterMean(j+1,:,:),3)); %trial type is in 2,3  
    meanF = meanF.*cellsMatch; % only count the cells with matching ROIs
    expMask = ~isnan(meanF);
    maxMean = squeeze(max(rasterMean(j+1,:,:),[],3)./meanF);
    
    F = squeeze(rasterMean(j+1,:,:));
    F = F.*repmat(cellsMatch,trialLength,1)';
    
    sessionF(:,:,i,j) = F;
    sessionMeanF(:,i,j) = meanF;
    sessionExpMask(:,i,j) = expMask;
    sessionMeanAct(:,i,j) = maxMean;
    %figure; plot(trialAct,'o');
    end
    toc
end

activeCells = sessionMeanF>0;

% Assumes if cell is R active, it will also show 
% up in L (at least won't be NaN)
numSessionActive = sum(activeCells(:,:,1),2);

% Calculate L R Activity
cellMeanAct = squeeze(nanmean(sessionMeanAct,2));
cellMeanActR = cellMeanAct(:,1);
cellMeanActL = cellMeanAct(:,2);
rlActFactor = cellMeanActR - cellMeanActL;


%% RL statistics plots

figure;plot(cellMeanActR, cellMeanActL,'o');
xlabel('Right Activity (\Delta F/F)');
ylabel('Left Activity (\Delta F/F)');
figure; histogram(rlActFactor);
xlabel('R L Activity Factor (\Delta(\Delta F / F))');
ylabel('Num of Cells');


%% Activity Plots

figure;
subplot(2,1,1); imagesc(sessionMeanAct(:,:,1));
subplot(2,1,2); imagesc(sessionMeanAct(:,:,2));
colormap jet;

figure; imagesc(activeCells(:,:,1));
xlabel('Session Number');
ylabel('Cells');
title('Active Cells by Session (R)');

figure; histogram(numSessionActive);
xlabel('Number of Active Sessions (R)');
ylabel('Number of Cells');

%% Extract ROIs

comX = nan(955,33,4);
comY = nan(955,33,4);
%sessionsActive = zeros(955,33);
dummyImgX = repmat(1:512,512,1);
dummyImgY = repmat((1:512)',1,512);

for cIdx = 1:numCells
    sessionsActive{cIdx} = find(activeCells(cIdx,:,1) | activeCells(cIdx,:,2));
    numActSessions(cIdx) = length(sessionsActive{cIdx});
    tic
    if numActSessions(cIdx) > 0
        for z = 1:4            
            for sIdx = 1:length(sessionsActive{cIdx})
                sessionId = sessionsActive{cIdx}(sIdx);
                cellMask = squeeze(roi_im_all{sessionId}(z,:,:)==cIdx);
                nPixCell = sum(cellMask(:));
                if nPixCell > 0
                    comX(cIdx,sessionId,z) = sum(dummyImgX(cellMask))/nPixCell;
                    comY(cIdx,sessionId,z) = sum(dummyImgY(cellMask))/nPixCell;
                end
            end
        end
    end
    toc
end
%% Plot Cell Maps vs Session
figure; hold on;
cmap = jet(numCells);
for cIdx = 1:numCells
    plot(comX(cIdx,:),comY(cIdx,:),'o',...
        'MarkerFaceColor',cmap(cIdx,:,:),...
        'Color',cmap(cIdx,:,:),...
        'MarkerSize',2);
end

axis square;
title('Cell Positions, all sessions, all planes');
%com(cellInx,zplane,sessionIdx)
%figure; imagesc(squeeze(cellMask(z,:,:)));

%% Averaged over Sessions and Planes
% Colors by number of active sessions
comXavg = squeeze(nanmean(nanmean(comX,2),3));
comYavg = squeeze(nanmean(nanmean(comY,2),3));

h = figure('Position', [100, 100, 1049, 895]);
hold on;
axis([0 512 0 512]);
xlabel('Pixels'); ylabel('Pixels');
box on;


% 2 color map for R L Activity factor
cmax = max(abs(rlActFactor));
%cmap = cool(200);
%rlColor = round(rlActFactor/cmax*99)+101;
dimPower = 0.7;
rldim = (abs(rlActFactor)/cmax).^dimPower;
rldim(isnan(rldim)) = 0;
sizeBoost = 3;

for cIdx = 1:numCells    
    if numSessionActive(cIdx) > 0 %&& ~isnan(rlColor(cellIdx))        
        if isnan(rlActFactor(cIdx))  
            rlColor = [.5 .5 .5];
        elseif rlActFactor(cIdx) < 0
            rlColor = [1 0 0] + [0 1 1]*(1 - rldim(cIdx));
        else
            rlColor = [0 0 1] + [1 1 0]*(1 - rldim(cIdx));
            %}
        end
        plot(comXavg(cIdx),comYavg(cIdx),'o',...
            'MarkerFaceColor',rlColor,...
            'Color','k',...%rlColor,...
            'MarkerSize',numSessionActive(cIdx)^.5*sizeBoost);
    end
end


% Manual Legend
numEx = 5;
yLoc = linspace(500,400,numEx);
xLoc = linspace(400,400,numEx);
rectangle('Position',[388 388 512 512])

numSessionActLeg = linspace(1,max(numSessionActive),numEx);
rlActFactorLeg = linspace(-cmax, cmax, numEx);
sizeLeg = numSessionActLeg.^.5*sizeBoost;
rldimLeg = (abs(rlActFactorLeg)/cmax).^dimPower;

for i = 1:numEx
    if rlActFactorLeg(i) < 0
        rlColor = [1 0 0] + [0 1 1]*(1 - rldimLeg(i));
    else
        rlColor = [0 0 1] + [1 1 0]*(1 - rldimLeg(i));
    end
    plot(xLoc(i),yLoc(i),'o',...
        'MarkerSize',sizeLeg(i),...
        'MarkerFaceColor',rlColor,...
        'Color', 'k');
    text(xLoc(i)+20,yLoc(i),[num2str(numSessionActLeg(i)) ' Sessions, '...
        num2str(rlActFactorLeg(i),2) ' Act']);
end

%% Do closeness in sequence correlate with anatomical closeness

% sessionF(cellID,t,sessionID,RL)
% want lists of distances, and matching list of correlations for each
% session

% Need to run ROI plot section first
allCells = 1:numCells;
cellMask = numActSessions > 0;
goodCells = allCells(cellMask);
numGoodCells = length(cellId);

pos = nan(2, numGoodCells);
actR = nan(trialLength,numGoodCells);
actL = nan(trialLength,numGoodCells);
tmR = nan(numGoodCells,1); %time of maximum
tmL = nan(numGoodCells,1);


%%


for cIdx = 1:numGoodCells
    tic
    cellId = goodCells(cIdx);
    for sIdx = 1:length(sessionsActive{cellId})
        sessionId = sessionsActive{cellId}(sIdx);
        xAvg = nanmean(comX(cellId,sessionId,:)); %avg across z planes
        yAvg = nanmean(comY(cellId,sessionId,:));
        pos(:,cIdx) = vertcat(xAvg, yAvg);
        
        actR(:,cIdx) = sessionF(cellId,:,sIdx,1);
        [peakR, tmR(cIdx)] = max(actR(:,cIdx));
        if tmR(cIdx) < 2, tmR(cIdx) = NaN; end %disqualify peaks at 1st frame
        actL(:,cIdx) = sessionF(cellId,:,sIdx,2);
        [peakL, tmL(cIdx)] = max(actL(:,cIdx));
        if tmL(cIdx) < 2, tmL(cIdx) = NaN; end
        clear dists;
        
        clear dists tdistsR tdistsL;
        dists = pdist(pos');
        tdistsR = pdist(tmR);
        tdistsL = pdist(tmL);
    end
    [RHO,PVAL] = corr(dists(~isnan(tdistsR))',tdistsR(~isnan(tdistsR))','type','Spearman');
    rho(sIdx) = RHO;
    toc
end
%%
dataBlock = squeeze(sessionF(goodCells,:,sIdx,1));
figure;imagesc(dataBlock); hold on;
plot(tmR,1:numGoodCells,'wo');

[RHO,PVAL] = corr(dists(~isnan(tdistsR))',tdistsR(~isnan(tdistsR))','type','Spearman')


figure; plot(dists,tdistsR,'b.');%,dists,tdistsL,'r.');

%%
for cIdx = 1:2%numCells
    sessionsActive{cIdx} = find(activeCells(cIdx,:,1) | activeCells(cIdx,:,2));
    numActSessions(cIdx) = length(sessionsActive{cIdx});
    tic
    if numActSessions(cIdx) > 0
        for z = 1:4            
            for sIdx = 1:length(sessionsActive{cIdx})
                cellMask = squeeze(roi_im_all{sessionsActive{cIdx}(sIdx)}(z,:,:)==cIdx);
                nPixCell = sum(cellMask(:));
                if nPixCell > 0
                    comX(cIdx,sIdx,z) = sum(dummyImgX(cellMask))/nPixCell;
                    comY(cIdx,sIdx,z) = sum(dummyImgY(cellMask))/nPixCell;
                end
            end
        end
    end
    toc
end
    
    