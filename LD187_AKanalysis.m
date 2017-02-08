% AK 20160923 
% Load and process Laura's 2P data for LD187

%% Load Data Files
%cd('...LDMS2/MATLAB/');
sessionFiles = dir('session_data/LD187_*.mat'); 
%dataF(cellID,t_pt,sessionID,rlID)
dataF = SAK_load_data(sessionFiles);

%% Parse Cell Matching Scores
load('labels_all.mat'); 
% labels_all, (1)same (2)probably same (3)not same cell 
% (4)filled same cell (5)filled different cell (6) hard to see

% Identify only 1's and 2's for further analysis
sessionMask = labels_all == 1 | labels_all == 2;
cellID = find(nanmean(sessionMask,1)>0); % cellID now has only good cells.
numCells = length(cellID);
dataF = dataF(cellID,:,:,:); % dataF now only has good cells


%% Mean Activity
sessionMeanF = squeeze(nanmean(dataF,2));
cellMeanF = squeeze(nanmean(sessionMeanF,2));
cellMaxF = squeeze(nanmax(sessionMeanF,[],2));
rlActFactor = cellMaxF(:,1)./cellMeanF(:,1)-cellMaxF(:,2)./cellMeanF(:,2);
activeF = sessionMeanF(:,:,1)>0;
%% RL statistics plots
figure;plot(cellMeanF(:,1),cellMeanF(:,2),'o');
xlabel('Right Activity (\Delta F/F)');
ylabel('Left Activity (\Delta F/F)');
figure; histogram(rlActFactor);
xlabel('R L Activity Factor (\Delta(\Delta F / F))');
ylabel('Num of Cells');

%% Activity Plots
figure; imagesc(sessionMeanF(:,:,1)>0);
xlabel('Session Number');
ylabel('Cells');
title('Active Cells by Session (R)');

figure; histogram(squeeze(sum(sessionMeanF(:,:,1)>0,2)));
xlabel('Number of Active Sessions (R)');
ylabel('Number of Cells');

%% Extract ROIs

cellLocs = extractROIcenters('wei_masks.mat',cellID,activeF);
%figure;plot(cellLocs(:,1),cellLocs(:,2))


%% Averaged over Sessions and Planes
% Size by number of active sessions
numSessionActive = squeeze(sum(activeF>0,2));

h = figure('Position', [100, 100, 1049, 895]); hold on;
axis([0 512 0 512]); xlabel('Pixels'); ylabel('Pixels');box on;

% Plot vasculature image in background
% vas = flipud(imread('vasculature.tif')); imagesc(vas); colormap gray;

% 2 color map for R L Activity factor
cmax = max(abs(rlActFactor));
dimPower = 1;
rldim = (abs(rlActFactor)/cmax).^dimPower; rldim(isnan(rldim)) = 0;
sizeBoost = 3;

for cIdx = 1:length(cellID)
        if rlActFactor(cIdx) < 0
            rlColor = [1 0 0] + [0 1 1]*(1 - rldim(cIdx));
        else
            rlColor = [0 0 1] + [1 1 0]*(1 - rldim(cIdx));
        end
        plot(cellLocs(cIdx,1),cellLocs(cIdx,2),'o',...
            'MarkerFaceColor',rlColor,...
            'Color','k',...%rlColor,...
            'MarkerSize',numSessionActive(cIdx)^.5*sizeBoost);
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

%% Plot sequences averaged over sessions

% For now, average all session together
trialMean = squeeze(nanmean(dataF,3));
zTrialMean = zscore(trialMean,0,2);
%figure; imagesc(zTrialMean(:,:,1));


% find peaks
[peaks, tPeaks] = max(zTrialMean,[],2); tPeaks = squeeze(tPeaks);
[temp,cellOrderR] = sort(tPeaks(:,1),1);
[temp,cellOrderL] = sort(tPeaks(:,2),1);

figure; colormap jet; % Just a figure to look at sequences averaged over all sessions
subplot(2,2,1);imagesc(zTrialMean(cellOrderR,:,1));title('R sorted by R');
subplot(2,2,2);imagesc(zTrialMean(cellOrderR,:,2));title('L sorted by R');
subplot(2,2,3);imagesc(zTrialMean(cellOrderL,:,1));title('R sorted by L');
subplot(2,2,4);imagesc(zTrialMean(cellOrderL,:,2));title('L sorted by L');

%% Are cells that cluster functionally close to each other?
% frames 1-26: 12 fr before trial onset to 12 fr after running
% frames 27-51: 12 fr before cue to 12 fr after cue
% frames 52-76: 12 frame before end to 12 fr after end
% correspond to values 1,2,3

for rl = 1:2
cellCategory(:,rl) = 1*(tPeaks(:,rl) < 27)+...
2*(tPeaks(:,rl) >= 27 & tPeaks(:,rl) <= 52)+...
3*(tPeaks(:,rl) >= 53);
end

figure; hold on;
ttl = {'R categorization'; 'L categorization'}; 
cmap = lines(3);

for rl = 1:2 
subplot(1,2,rl); 
scatter(cellLocs(cellID,1),cellLocs(cellID,2),[],cmap(cellCategory(:,rl),:,:),'filled');
title(ttl{rl}); axis square;
end        

figure; hold on;
for catIdx = 1:3
    for rl = 1:2
        ct = (catIdx-1)*2+rl;
        dist_within{catIdx,rl} = pdist(cellLocs(cellCategory(:,rl)==catIdx,:));
        dist_between{catIdx,rl} = pdist2(cellLocs(cellCategory(:,rl)==catIdx,:),cellLocs(cellCategory(:,1)~=catIdx,:));
        subplot(3,2,ct);hold on; 
        h1(ct) = histogram(dist_within{catIdx,rl});hold on; h1(ct).Normalization = 'pdf';
        h2(ct) = histogram(dist_between{catIdx,rl}(:)); h2(ct).Normalization = 'pdf';
    end
end

%% Are cells that overlap in their activity closer anatomically
% correlation vs anatomical distance
% Specifically interested in cells with functionally very similar responses
figure; hold on;
for rl = 1:2
    corrDistF{rl}=pdist(zTrialMean(:,:,rl),'correlation');
    anatDist = pdist(cellLocs(cellID,:),'euclidean');
    
end
    corr_anat_h = histogram2(corrDistF{1},anatDist,'DisplayStyle','tile');
    xlabel('Correlation Distance'); ylabel('Anatomical Distance (um)');
    corr_anat_h.NumBins = [40 40];
  
corrVals = 0:.1:2;clear anatDistCumCorr;
for rl = 1:2
    for i = 1:length(corrVals)
        anatDistCumCorr{rl}(i) = nanmean(anatDist(corrDistF{rl}<corrVals(i)));
        anatDistCumCorr75{rl}(i) = prctile(anatDist(corrDistF{rl}<corrVals(i)),75);
        anatDistCumCorr25{rl}(i) = prctile(anatDist(corrDistF{rl}<corrVals(i)),25);
    end
end

figure;
for rl = 1:2
    subplot(1,2,rl); plot(corrVals,anatDistCumCorr{rl});
    xlabel('Correlation Distance');ylabel('Anatomical Distance (um)');
    hold on; plot(corrVals,anatDistCumCorr75{rl},'k-');
    plot(corrVals,anatDistCumCorr25{rl},'k-');
end
%% Do closeness in sequence correlate with anatomical closeness?

% dataF(cellID,t,sessionID,RL)
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
        
        actR(:,cIdx) = dataF(cellId,:,sIdx,1);
        [peakR, tmR(cIdx)] = max(actR(:,cIdx));
        if tmR(cIdx) < 2, tmR(cIdx) = NaN; end %disqualify peaks at 1st frame
        actL(:,cIdx) = dataF(cellId,:,sIdx,2);
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
dataBlock = squeeze(dataF(goodCells,:,sIdx,1));
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
  

    
