clc; clear;
addpath(genpath('/project/joelvoss/tmp-rostowsky/bosc/'));
addpath(genpath('/project/joelvoss/tmp-rostowsky/github/retinotopy'));

%%
indir = '/project/joelvoss/tmp-rostowsky/hpcData';
subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);
inFile = 'tfAnalysis.mat';

%%
for currSubject = 1:length(subjects)
    tic;
    fprintf([subjects{currSubject} '\n']);
    load([indir '/' subjects{currSubject} '/' inFile]);
    outdir = '/project/joelvoss/tmp-rostowsky/hpcData';

    %% find the refresh rate of the monitor and the ieeg recording. Find the number of frames (at the sampling
    % rate of the ieeg recording) per event
    eventPeriodinFrames = 50; % this was determined manually as the smallest number of frames for which a stimulus
    % passed over a given coordinate
    % here this is the bar stimulus which
    % continuously stimulates for 50 frames

    % eventPeriodinFrames = eventPeriodinFrames * 2; % Jim says an event should be defined as 2x the length for which a given coordinate is "on"
    eventPeriodinFrames = eventPeriodinFrames * 1; % 

    stimulusFrameLength = diff(hippocampalStimTimes); % this is the refresh rate of the monitor
    stimulusFrameLength = stimulusFrameLength(1);

    samplesPerFrame = diff(hippocampalTimes); % this is the sampling rate
    samplesPerFrame = samplesPerFrame(1);

    eventPeriodinSeconds = eventPeriodinFrames * stimulusFrameLength; % this is the total length of an event

    framesPerEvent = ceil(eventPeriodinSeconds / samplesPerFrame); % this is the number of frames (at the sampling rate) that fit into an event

    %%
    [onBlocks, offBlocks] = makeBlocks(sortedBlockInds, sortedBlockIndsData, sortedOnOffIdData);

    %% find indices from the on data for each event
    onBlockNames = fieldnames(onBlocks);
    onEventIndices = cell(length(onBlockNames), 1);
    onEventIndicesOff = cell(length(onBlockNames), 1);

    for j = 1:length(onBlockNames)
        %onEventIndices{j,1} = slidingWindowIndices(onBlocks.(onBlockNames{j}).allInds, framesPerEvent, 0);
        onEventIndices{j,1} = slidingWindowIndices(onBlocks.(onBlockNames{j}).onInds, framesPerEvent, 0, 1);
        onEventIndicesOff{j,1} = slidingWindowIndices(onBlocks.(onBlockNames{j}).offInds, framesPerEvent, 0, 1);
    end

    onEventTimes = cell(onEventIndices);
    for j = 1:length(onEventTimes)
        for jj = 1:size(onEventTimes{j},2)
            currIndices = onEventIndices{j}(:,jj);
            currIndices = currIndices(currIndices ~= 0);
            extractedTimes = hippocampalTimes(currIndices);
            numZeros = size(onEventTimes{j}(:,jj),1) - length(extractedTimes);
            if numZeros ~= 0
                onEventTimes{j}(:,jj) = [extractedTimes zeros(1, numZeros)];
            else
                onEventTimes{j}(:,jj) = extractedTimes;
            end
        end
    end

    onEventTimesOff = cell(onEventIndicesOff);
    for j = 1:length(onEventTimesOff)
        for jj = 1:size(onEventTimesOff{j},2)
            currIndices = onEventIndicesOff{j}(:,jj);
            currIndices = currIndices(currIndices ~= 0);
            extractedTimes = hippocampalTimes(currIndices);
            numZeros = size(onEventTimesOff{j}(:,jj),1) - length(extractedTimes);
            if numZeros ~= 0
                onEventTimesOff{j}(:,jj) = [extractedTimes zeros(1, numZeros)];
            else
                onEventTimesOff{j}(:,jj) = extractedTimes;
            end
        end
    end

    %% find indices from the off data for each event
    offBlockNames = fieldnames(offBlocks);
    offEventIndices = cell(length(offBlockNames), 1);
    for j = 1:length(offBlockNames)
        offEventIndices{j,1} = slidingWindowIndices(offBlocks.(offBlockNames{j}), framesPerEvent, 0, 0);
    end

    offEventTimes = cell(offEventIndices);
    for j = 1:length(offEventTimes)
        for jj = 1:size(offEventTimes{j},2)
            currIndices = offEventIndices{j}(:,jj);
            currIndices = currIndices(currIndices ~= 0);
            extractedTimes = hippocampalTimes(currIndices);
            numZeros = size(offEventTimes{j}(:,jj),1) - length(extractedTimes);
            if numZeros ~= 0
                offEventTimes{j}(:,jj) = [extractedTimes zeros(1, numZeros)];
            else
                offEventTimes{j}(:,jj) = extractedTimes;
            end
        end
    end

    %% plot a subsection of eventing for visualization
    plotSelection = [onEventIndices{:}];
    %numEvents = size([onEventIndices{:}],2);
    %plotSelection = onEventIndices{1}(:,1:numEvents);
    intervalLength = numel(unique(plotSelection(:)));

    currFig = figure('visible', 'off'); hold on;
    plot(1:length(hippocampalMaskIndData), hippocampalMaskIndData);
    yValMax = max(hippocampalMaskIndData);
    xcount = 0;
    for j = 1:size(plotSelection,2)
        currX = plotSelection(:,j);
        bufferCheck = currX(currX ~= 0);
        currX = [bufferCheck(1) bufferCheck(1) bufferCheck(end) bufferCheck(end) ];
        currY = [0 yValMax yValMax 0];
        fill(currX, currY, [0.3, 0.3, 0.3], 'EdgeColor', 'k', 'EdgeAlpha', 0.3, 'FaceAlpha', 0.3);
        x_label = currX(end);
        y_label = 1.1;
        text(x_label, y_label, ['event' num2str(j) '-end'], 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'Rotation', 90);
    end
    saveas(currFig, [figDir '/' subjects{currSubject} '/' 'eventFigure.png']);
    close();

    %% create design matrix

    % this holds descriptors of data
    fullDescriptionMat = zeros(length(hippocampalTimes), 5);
    fullDescriptionMatLabels = cell(1, size(fullDescriptionMat, 2));
    fullDescriptionMatLabels{1} = 'blockNumber';
    fullDescriptionMatLabels{2} = 'stimOnOff';
    fullDescriptionMatLabels{3} = 'maskIndex';
    fullDescriptionMatLabels{4} = 'imageIndex';
    fullDescriptionMatLabels{5} = 'time';

    fullDescriptionMat(:,1) = sortedBlockIndsData;
    fullDescriptionMat(:,2) = sortedOnOffIdData;
    fullDescriptionMat(:,3) = hippocampalMaskIndData;
    fullDescriptionMat(:,4) = hippocampalImageIndData;
    fullDescriptionMat(:,5) = hippocampalTimes;

    clearvars -except makeEvents onEventTimesOff onEventIndicesOff currSubject indir inFile outdir subjects hippocampalDataReref newContactList out onEventGammaData offEventGammaData gammaFreqs offMeans onEventData offEventData discharges fullDescriptionMat fullDescriptionMatLabels fullDataMat ...
        resultLabel res onEventTimes onEventIndices offEventTimes offEventIndices hippocampalImageIndData hippocampalMaskIndData

    %% make an array to hold spike data for each electrode
    timeCol = find(strcmpi(fullDescriptionMatLabels, 'time'));
    sampleFreq = abs(fullDescriptionMat(1, timeCol) - fullDescriptionMat(2, timeCol));

    % this gets the total number of events
    numEvents = 0;
    for j = 1:length(onEventTimes)
        currNumEvents = size(onEventTimes{j},2);
        numEvents = numEvents + currNumEvents;
    end
    for j = 1:length(onEventTimesOff)
        currNumEvents = size(onEventTimesOff{j},2);
        numEvents = numEvents + currNumEvents;
    end
    for j = 1:length(offEventTimes)
        currNumEvents = size(offEventTimes{j},2);
        numEvents = numEvents + currNumEvents;
    end

    spikeArray = zeros(numEvents, size(hippocampalDataReref,1));

    for currElect = 1:size(hippocampalDataReref)
        %hippocampalDischargeContactInds = find(ismember(newContactList, resultLabel));
        hippocampalDischargeContactInds = find(ismember(newContactList, resultLabel{currElect}));
        hippocampalDischargeContactInds = find(ismember(out.chan, hippocampalDischargeContactInds));
        spikeTimes = out.pos(hippocampalDischargeContactInds);

        % now we make an array that contains any spikes within the current hippocampal
        % channel
        spikeWindows = [];
        for j = 1:length(spikeTimes)
            currSpikeTime = spikeTimes(j);
            addTimes = spikeTimes(j) - 1:sampleFreq:spikeTimes(j) + 1;
            spikeWindows = [spikeWindows, addTimes];
        end
        spikeWindows = unique(round(spikeWindows,4));

        countEvent = 0;
        checkEpilepsy = zeros(numEvents,1);
        for j = 1:length(onEventTimes)
            for k = 1:size(onEventTimes{j},2)
                countEvent = countEvent + 1;
                currTimes = onEventTimes{j}(:,k);
                isSpikeAffected = sum(ismember(currTimes, spikeWindows));
                if isSpikeAffected
                    checkEpilepsy(countEvent) = 1;
                end
            end
        end
        for j = 1:length(onEventTimesOff)
            for k = 1:size(onEventTimesOff{j},2)
                countEvent = countEvent + 1;
                currTimes = onEventTimesOff{j}(:,k);
                isSpikeAffected = sum(ismember(currTimes, spikeWindows));
                if isSpikeAffected
                    checkEpilepsy(countEvent) = 1;
                end
            end
        end
        for j = 1:length(offEventTimes)
            for k = 1:size(offEventTimes{j},2)
                countEvent = countEvent + 1;
                currTimes = offEventTimes{j}(:,k);
                isSpikeAffected = sum(ismember(currTimes, spikeWindows));
                if isSpikeAffected
                    checkEpilepsy(countEvent) = 1;
                end
            end
        end
        spikeArray(:,currElect) = checkEpilepsy;
    end

    %
    descriptorMat = zeros(numEvents, 4); % event block, event on block/off block logic, first index(this is used for ordering the events), is full event
    descriptorMatLabels = cell(1, size(descriptorMat,2)); % this is what each column describes so you don't forget
    descriptorMatLabels{1} = 'blockNumber';
    descriptorMatLabels{2} = 'onOff';
    descriptorMatLabels{3} = 'firstIndex';
    descriptorMatLabels{4} = 'isFullEvent';

    addCount = 0;
    for j = 1:length(onEventTimes)
        eventBlock = ones(size(onEventTimes{j},2),1) .* j;
        eventOnOffBlock = ones(size(onEventTimes{j},2),1);
        eventFirstIndex = onEventIndices{j}(1,:)';
        descriptorMat(addCount + 1: addCount + length(eventBlock), 1) = eventBlock;
        descriptorMat(addCount + 1: addCount + length(eventOnOffBlock), 2) = eventOnOffBlock;
        descriptorMat(addCount + 1: addCount + length(eventOnOffBlock), 3) = eventFirstIndex;
        currCheckFull = zeros(size(onEventTimes{j},2), 1);
        for jj = 1:size(onEventTimes{j},2)
            currCheck = onEventIndices{j}(:,jj);
            if length(currCheck(currCheck ~= 0)) == length(currCheck)
                currCheckFull(jj) = 1;
            else
                currCheckFull(jj) = 0;
            end
        end
        descriptorMat(addCount + 1: addCount + length(eventOnOffBlock), 4) = currCheckFull;
        addCount = addCount + length(eventBlock);
    end
    for j = 1:length(onEventTimesOff)
        eventBlock = ones(size(onEventTimesOff{j},2),1) .* j;
        eventOnOffBlock = ones(size(onEventTimesOff{j},2),1) * 999; % this is just how I distinguish on/off within an on block
        eventFirstIndex = onEventIndicesOff{j}(1,:)';
        descriptorMat(addCount + 1: addCount + length(eventBlock), 1) = eventBlock;
        descriptorMat(addCount + 1: addCount + length(eventOnOffBlock), 2) = eventOnOffBlock;
        descriptorMat(addCount + 1: addCount + length(eventOnOffBlock), 3) = eventFirstIndex;
        currCheckFull = zeros(size(onEventTimesOff{j},2), 1);
        for jj = 1:size(onEventTimesOff{j},2)
            currCheck = onEventIndicesOff{j}(:,jj);
            if length(currCheck(currCheck ~= 0)) == length(currCheck)
                currCheckFull(jj) = 1;
            else
                currCheckFull(jj) = 0;
            end
        end
        descriptorMat(addCount + 1: addCount + length(eventOnOffBlock), 4) = currCheckFull;
        addCount = addCount + length(eventBlock);
    end
    for j = 1:length(offEventTimes)
        eventBlock = ones(size(offEventTimes{j},2),1) .* j;
        eventOnOffBlock = zeros(size(offEventTimes{j},2),1);
        eventFirstIndex = offEventIndices{j}(1,:)';
        descriptorMat(addCount + 1: addCount + length(eventBlock), 1) = eventBlock;
        descriptorMat(addCount + 1: addCount + length(eventOnOffBlock), 2) = eventOnOffBlock;
        descriptorMat(addCount + 1: addCount + length(eventOnOffBlock), 3) = eventFirstIndex;
        currCheckFull = zeros(size(offEventTimes{j},2), 1);
        for jj = 1:size(offEventTimes{j},2)
            currCheck = offEventIndices{j}(:,jj);
            if length(currCheck(currCheck ~= 0)) == length(currCheck)
                currCheckFull(jj) = 1;
            else
                currCheckFull(jj) = 0;
            end
        end
        descriptorMat(addCount + 1: addCount + length(eventOnOffBlock), 4) = currCheckFull;
        addCount = addCount + length(eventBlock);
    end

    %% make event stimuli
    makeEvents = 1;
    if makeEvents
    maskStimuli = load('/project/joelvoss/code/retinotopy/workspace_retinotopyCaltsmash.mat');
    maskStimuli = maskStimuli.maskimages;
    maskStimuli = maskStimuli(1:4:end, 1:4:end,:);
    blankMask = zeros(size(maskStimuli,1), size(maskStimuli,2));
    binaryMasks = zeros(size(maskStimuli,1), size(maskStimuli,2), size(maskStimuli,3));

    for j = 1:size(binaryMasks,3)
        currMask = maskStimuli(:,:,j);
        currNewMask = blankMask;
        nonZeroInds = find(currMask(:) ~= 0);
        currNewMask(nonZeroInds) = 1;
        binaryMasks(:, :, j) = currNewMask;
    end

    onEventImages = cell(size(onEventIndices));
    for j = 1:length(onEventImages)
        currIndices = onEventIndices{j};
        currBlockMasks = zeros(size(maskStimuli,1), size(maskStimuli,2), size(currIndices, 2));
        parfor jj = 1:size(currIndices, 2)
            currMaskInds = currIndices(:, jj);
            currMaskInds = currMaskInds(currMaskInds ~= 0);
            currMaskInds = hippocampalMaskIndData(currMaskInds);

            currEventMasks = zeros(size(maskStimuli,1), size(maskStimuli,2), length(currMaskInds));
            for k = 1:length(currMaskInds)
                if currMaskInds(k) == 0
                    currEventMasks(:,:,k) = blankMask;
                else
                    currEventMasks(:,:,k) = binaryMasks(:,:,currMaskInds(k));
                end
            end
            currEventMasks = mean(currEventMasks,3);
            currBlockMasks(:,:,jj) = currEventMasks;
            %figure; imagesc(currEventMasks);
        end
        onEventImages{j} = currBlockMasks;
    end
    end
    
    %%
    clearvars -except currSubject indir inFile outdir subjects hippocampalDataReref newContactList out onEventGammaData offEventGammaData gammaFreqs offMeans onEventData offEventData discharges fullDescriptionMat fullDescriptionMatLabels fullDataMat resultLabel ...
        onEventTimes onEventIndices onEventTimesOff onEventIndicesOff offEventTimes onEventImages offEventIndices descriptorMat descriptorMatLabels descriptorMatNonSpike powerMat powerMatNonSpike averageGammaPower makeEvents spikeArray
    
    makeEvents = 1;
    if makeEvents
        clearvars res
        save([outdir '/' subjects{currSubject} '/events.mat'], '-v7.3');
    end
    toc;
end

%%
function [outEvents] = slidingWindowIndices(inputArray, windowLength, overlap, onlyOn)
% this function will return the indices of the hippocampal data for each
% event

if onlyOn
    blockEnds = diff(inputArray); %you need to find stops within a block
    blockEnds = find(blockEnds ~= 1);
    % numEvents = ceil(length(inputArray) / windowLength);
    outEvents = zeros(windowLength, 1);
    eventCount = 0;
    for j = 1:length(blockEnds)+1
        if j == 1
            currOutEvents = buffer(inputArray(1:blockEnds(j)), windowLength, 0, 'nodelay');
            outEvents(:, (eventCount + 1):(eventCount + size(currOutEvents, 2))) = currOutEvents;
        elseif j == length(blockEnds)+1
            currOutEvents = buffer(inputArray(blockEnds(end)+1:end), windowLength, 0, 'nodelay');
            outEvents(:, (eventCount + 1):(eventCount + size(currOutEvents, 2))) = currOutEvents;
        else
            currOutEvents = buffer(inputArray(blockEnds(j-1)+1:blockEnds(j)), windowLength, 0, 'nodelay');
            outEvents(:, (eventCount + 1):(eventCount + size(currOutEvents, 2))) = currOutEvents;
        end
        eventCount = eventCount + size(currOutEvents, 2);
    end
else
    outEvents = buffer(inputArray, windowLength, ceil(windowLength*overlap), 'nodelay');
end

end