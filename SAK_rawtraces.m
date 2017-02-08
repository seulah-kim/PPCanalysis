%% SAK_rawtraces.m
lalalal
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

%% Parse Cell Matching Scores------------------------------------------
load('labels_all.mat'); 
% labels_all, (1)same (2)probably same (3)not same cell 
% (4)filled same cell (5)filled different cell (6) hard to see

% Identify only 1's and 2's for further analysis
sessionMask = labels_all == 1 | labels_all == 2;
healthy_any = find(sum(sessionMask,1)>0); % These cells were healthy at least 1 day
healthy_10 = find(sum(sessionMask,1)>3);     % These cells were healthy at least 3 days (~10%)
healthy_50 = find(sum(sessionMask,1)>16);     % These cells were healthy at least 16 days (~50%)
healthy_all = find(sum(sessionMask,1)==33);     % These cells were healthy all 33 days (100%)

cellID = healthy_any;		% This set will be used for the rest of analysis 
numCells = length(cellID);


% Only look at Session 1 throughout the script--------------------------------

colorlabels = {'r','b','g','k','m'};

sessionFiles = dir('session_data/LD187_*.mat');
sessionNum = 5;
sessionFiles = sessionFiles(sessionNum);	% Session number of interest

if RAWINDIVTRACE == 1

for j = 1:2 %map to 1-R 2-L
    load(['session_data/' sessionFiles(1).name],'rasterActivity');
    F = rasterActivity{j+1};	% cell ID, trials, time
%   row = find(~isnan(F(:,1,1)));    % find rows that contain number
    F = F(cellID,:,:);		     % select only good cells based on mask (across sessions)
%    [row,col] = find(~isnan(iF(:,:,1)));    % find rows that contain number
    [row,col] = find(nanmean(F,3)>0);    % find rows&cols that contain number (active cells/trial comb)
    ActiveCells = nanmean(F,3)>0;    % template for active cells (should we impose threshold?)
    row = find(sum(ActiveCells,2)>0);
    if j ==1
       cell_idx{j} = sort(datasample(row,5,'Replace',false));		% 5 random cell IDs selected from [row]
    else
         if ~isempty(intersect(row,cell_idx{j-1}))
	  overlapCell = intersect(row,cell_idx{j-1});
	      if length(overlapCell)<5
 	  	addCell = sort(datasample(setdiff(row,overlapCell),5-length(overlapCell)));
                cell_idx{j} = [overlapCell,addCell];
              elseif length(overlapCell)==5
                cell_idx{j} = overlapCell;
	      end
         else
          cell_idx{j} = datasample(row,5);
         end
     end
% Raw single trial trace---------------------------------------------
    temp_cell = cell_idx{j};
    figure
    col_array = [];
        for n = 1:length(temp_cell)
		%col = find(~isnan(F(cellID(n),:,1)));	% find columns that contain number
    		col_idx = col(row==temp_cell(n));	% find columns that contain number
		col_idx = datasample(col_idx,1);	% pick 1 random trial from a given cell
		randTrial = squeeze(F(temp_cell(n),col_idx,:));	
		
		plot(randTrial,colorlabels{n})
		hold on
		col_array = [col_array,col_idx];
	end
    xL = get(gca,'XLim');
    yL = get(gca,'YLim');
    Label= {sprintf('cell %d trial%d',temp_cell(1),col_array(1)),
    	    sprintf('cell %d trial%d',temp_cell(2),col_array(2)),
    	    sprintf('cell %d trial%d',temp_cell(3),col_array(3)),
    	    sprintf('cell %d trial%d',temp_cell(4),col_array(4)),
    	    sprintf('cell %d trial%d',temp_cell(5),col_array(5))};
    legend(Label)
    if j==1
	Dir = 'R';
    else
	Dir = 'L';
    end
    title(sprintf('raw traces (Stimulus-%s)',Dir))
    xlabel('Time frame')
    plot([13 13],yL,'Color','k','LineStyle',':')	% trial onset/cue onset??
%    plot([14 14],yL,'Color','k','LineStyle',':')	  
    area([39-12 39+12],[max(yL) max(yL)],'FaceColor','y','FaceAlpha',0.2, ...
	'EdgeColor','none','EdgeAlpha',0)
    
    plot([39 39],yL,'Color','k','LineStyle',':')	% cue offset / delay period onset
    area([54-12 54+12],[max(yL) max(yL)],'FaceColor','b','FaceAlpha',0.2, ...
	'EdgeColor','none','EdgeAlpha',0)
    plot([54 54],yL,'Color','k','LineStyle',':')	% trial end / turn onset
    ylabel('Spike Activity (a.u)')
    saveas(gcf,sprintf('RawTrace_%s.png',Dir))
end
end	% End of loop for RawIndivTrace

commonCell = datasample(intersect(cell_idx{1},cell_idx{2}),1);
maxActiveTrial = find(sum(ActiveCells,1)==max(sum(ActiveCells,1)));		% which trial had the most number of active cells
maxActiveTrial = maxActiveTrial(1);

for jj = 1:2 %map to 1-R 2-L
    load(['session_data/' sessionFiles(1).name],'rasterActivity');
    F = rasterActivity{jj+1};	% cell ID, trials, time
    F = F(cellID,:,:);		     % select only good cells based on mask (across sessions)

% Average across trials given a stimulus. Trial variability per neuron
    TrialMean = nanmean(squeeze(F(commonCell,:,:)),1);
    TrialStd = nanstd(squeeze(F(commonCell,:,:)),1)./sqrt(sum(ActiveCells(commonCell,:)));	% standard error = standard dev/sqrt(N) 

% Average across cells given a stimulus.neuronal variability per trial
    PopulMean = nanmean(squeeze(F(:,maxActiveTrial,:)),1);
    PopulStd = nanstd(squeeze(F(:,maxActiveTrial,:)),1)./sqrt(sum(ActiveCells(:,maxActiveTrial)));	% standard error = standard dev/sqrt(N) 

% Average session activity, variability across cells and trials
    SessionMean = squeeze(nanmean(nanmean(F,1),2));
    SessionStd = squeeze(nanmean(nanmean(F,1),2))./sqrt(sum(sum(ActiveCells))); 
 
    if TRIALMEAN == 1
 	plot_timetrace(TrialMean,TrialStd,'trial',commonCell,sessionNum,jj)
    end
    
% Average across neurons given a trial. neuronal variability per trial-------------------------------
    if POPMEAN == 1
     	plot_timetrace(PopulMean,PopulStd,'cell',maxActiveTrial,sessionNum,jj)	% type is 'cell' because this is mean of many cells per trial

    if SESSMEAN == 1
    plot_timetrace(SessionMean,SessionStd,'session',sessionNum,sessionNum,jj)    


    end

end 

%% Extract Decision information from Behavioral Data

load(['session_data/' sessionFiles(1).name],'vData');


%correct{1,c1}(t) = (mod(c1,2)==1 && ...
%                        sum(VirmenCombined(2,trial_ind))<-sum(VirmenCombined(2,trial_ind)) ||...
%                        mod(c1,2)==0 && sum(VirmenCombined(2,trial_ind))>-sum(VirmenCombined(2,trial_ind)));
%	% c1 I think is the stimulus info, so if c1 is odd, L, if c1 is even, R. 
%correct{1,c1}(t) = (mod(c1,2)==1 && ...
%                        sum(VirmenCombined(2,trial_ind))>-sum(VirmenCombined(2,trial_ind)) ||...
%                        mod(c1,2)==0 && sum(VirmenCombined(2,trial_ind))<-sum(VirmenCombined(2,trial_ind)));
%
