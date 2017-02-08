function spreadsheet = SAK_read_data(sessionFiles)
% Read data information and serve as key

%% Params-----------------------------------------------------------------
numSessions= length(sessionFiles);
numCells = 955;		% total numb of cells
numPxl = 512;           % ??  
trialLength = 76;       % dimension of time
direction = 2;		% R and L

spreadsheet = cell(numSessions,direction,3);
%outmat = cell(numSessions,direction);

for i = 1:numSessions
    display(['Loading Session ' num2str(i)]);
    for j = 2:3 %map to 2-R 3-L (cell 1 and 4 are empty)
    clear rasterMean
    load(['session_data/' sessionFiles(i).name],'rasterActivity');
    spreadsheet(i,j-1,:)= num2cell(size(rasterActivity{j}));
%    outmat(i,j-1) = mat2cell(rasterActivity{j});
    end
end

end

