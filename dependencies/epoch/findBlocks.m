function [blockBoundariesCell, blockBoundariesCellTimes] = findBlocks(stimIn,isOn,behavioralData, addBuffer, adjBehavFlag, currBlock)
% args
% stimIn indices for stimulus
% isOn: 0 for offStimulus, 1 for onStimulus
% behavioralData behavioralData structure
% addBuffer adds a 2s buffer to each end of a block

stimDiff = diff(stimIn);
blockBoundaries = find(stimDiff ~= 1);


blockBoundariesCell = cell(length(blockBoundaries)+1,1); % +1 is experimental
blockBoundariesCellTimes = cell(length(blockBoundaries)+1,1);  % +1 is experimental

for k = 1:length(blockBoundariesCell)
    if ~isOn
        if k-1 == 0
            if (adjBehavFlag && currBlock == 1)
                blockBoundariesCell{k} = stimIn(1):stimIn(blockBoundaries(k));
            else
                blockBoundariesCell{k} = 1:stimIn(blockBoundaries(k));
            end
        elseif k > length(blockBoundaries) %else if is experimental
            blockBoundariesCell{k} = stimIn(blockBoundaries(k-1)+1):stimIn(end);
        else
            blockBoundariesCell{k} = stimIn(blockBoundaries(k-1)+1):stimIn(blockBoundaries(k));
        end
    else
        if k == length(blockBoundariesCell)
            blockBoundariesCell{k} = stimIn(blockBoundaries(k-1)+1):stimIn(end);
        elseif k-1 == 0
            if (adjBehavFlag && currBlock == 1)
                blockBoundariesCell{k} = 1:stimIn(blockBoundaries(k));
            else
                blockBoundariesCell{k} = stimIn(1):stimIn(blockBoundaries(k));
            end
        else
            blockBoundariesCell{k} = stimIn(blockBoundaries(k-1)+1):stimIn(blockBoundaries(k));
        end
    end

end

%% find timestamps of stim off regions and add 2s buffer
behavioralDataTimeFrames = behavioralData.timeframes;

for j = 1:length(blockBoundariesCell)
    currTimeFrames = behavioralDataTimeFrames(blockBoundariesCell{j});
    if addBuffer
        if currTimeFrames(1) == 0
            currTimeFrames = [currTimeFrames(1):1/30:(currTimeFrames(end)+2)];
        else
            currTimeFrames = [(currTimeFrames(1)-2):1/30:(currTimeFrames(end)+2)];
        end
    end
    blockBoundariesCellTimes{j} = currTimeFrames;
end

end