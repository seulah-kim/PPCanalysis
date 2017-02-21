function plot_timetrace(Mean,Std,type,input,sessionNum,RorL)
% This function plots a timetrace for given combination of input variables. It will plot the mean and standard error mean (not sure if this is valid analysis method for deconvoluted data because its underlying distribution may not be Gaussian. Maybe plot median and 5 and 95 percentiles instead?)
% Description of inputs:
% 	type -	1) cell : cell average across trials
%		2) trial: trial average across cells
%		3) session: session average across trials and cells
% 	input - corresponding number associated with type
% 	RorL - 1 is R, 2 is L
 figure
    plot(squeeze(Mean),'k')	% mean trace
    hold on
    Mean = reshape(Mean,[1 76]); Std = reshape(Std,[1 76]);
    h =fill([1:76, fliplr(1:76)],[Mean+Std , fliplr(Mean-Std)],'k');	% Coloring stderr
    h.FaceAlpha = 0.2; h.EdgeColor = 'none'; h.EdgeAlpha = 0;
    
    xL = get(gca,'XLim');
    yL = get(gca,'YLim');
    if RorL==1
        Dir = 'R';
    else
        Dir = 'L';
    end
    
    if isequal(type,'session')	% define axes ranges based on plot type 
    axis([0 max(xL) [yL]]);
    else
    axis([0 max(xL) -0.002 max(yL)]);
    end

% Below assigns appropriate title and legend labels for plot type
    if isequal(type,'trial')
    Label= {sprintf('cell %d',input)};
    title(sprintf('Average across trial in Session %d (Stimulus-%s)',sessionNum, Dir))
 
    elseif isequal(type,'cell')
    Label= {sprintf('trial %d',input)};
    title(sprintf('Population average in Session %d (Stimulus-%s)',sessionNum,Dir))
   
    elseif  isequal(type,'session')
    Label= {sprintf('session %d',sessionNum)};
    title(sprintf('Average in Session %d (Stimulus-%s)',sessionNum, Dir))

    end

    legend(Label)
    xlabel('Time frame')
    plot(xL,[0 0],'Color','k','LineStyle',':')  % trial onset/cue onset??
    plot([13 13],yL,'Color','k','LineStyle',':')        % trial onset/cue onset??
%    plot([14 14],yL,'Color','k','LineStyle',':')         
%    area([39-12 39+12],[max(yL) max(yL)],'FaceColor','y','FaceAlpha',0.2, ...
%       'EdgeColor','none','EdgeAlpha',0)

    plot([39 39],yL,'Color','k','LineStyle',':')        % cue offset / delay period onset
%    area([54-12 54+12],[max(yL) max(yL)],'FaceColor','b','FaceAlpha',0.2, ...
%       'EdgeColor','none','EdgeAlpha',0)
    plot([54 54],yL,'Color','k','LineStyle',':')        % trial end / turn onset
    ylabel('Spike Activity (a.u)')
  
% Save current figure 
    if isequal(type,'trial')
    saveas(gcf,sprintf('AverageAcrossTrials_cell%d_%s.png',input,Dir))
 
    elseif isequal(type,'cell')
    saveas(gcf,sprintf('AverageAcrossCells_trial%d_%s.png',input,Dir))

    elseif isequal(type,'session')

    saveas(gcf,sprintf('AverageInSession%d_%s.png',input,Dir))
    end
end
