function cellLocs = extractROIcenters(masksFileName,cellID,activeF)

load(masksFileName); % ROI masks identified manually by students

comX = nan(955,33,4);
comY = nan(955,33,4);

dummyImgX = repmat(1:512,512,1);
dummyImgY = repmat((1:512)',1,512);

for cIdx = 1:length(cellID)
    sessionsActive = find(activeF(cIdx,:));
    display(['Extracting ROI for cell ' num2str(cIdx)]);
        for z = 1:4            
            for sIdx = 1:length(sessionsActive)
                sessionId = sessionsActive(sIdx);
                %note, roi_im_all is uses original cell ids (955)
                cellMask = squeeze(roi_im_all{sessionId}(z,:,:)==cellID(cIdx));
                nPixCell = sum(cellMask(:));
                if nPixCell > 0
                    comX(cIdx,sessionId,z) = sum(dummyImgX(cellMask))/nPixCell;
                    comY(cIdx,sessionId,z) = sum(dummyImgY(cellMask))/nPixCell;
                end
            end
        end

end

comXavg = squeeze(nanmean(nanmean(comX,3),2));
comYavg = squeeze(nanmean(nanmean(comY,3),2));
cellLocs =  horzcat(comXavg,comYavg);