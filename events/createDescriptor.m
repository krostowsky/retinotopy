clc; clear;
addpath(genpath('/project/joelvoss/tmp-rostowsky/bosc/'));
addpath(genpath('/project/joelvoss/tmp-rostowsky/github/retinotopy'));

%%
indir = '/project/joelvoss/tmp-rostowsky/hpcData';
subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);
inFile = 'tfAnalysis.mat';

%%
for currSubject = 2:length(subjects)
    tic;
    fprintf([subjects{currSubject} '\n']);
    load([indir '/' subjects{currSubject} '/' inFile]);
    outdir = '/project/joelvoss/tmp-rostowsky/hpcData';

    %%
    [onBlocks, offBlocks] = makeBlocks(sortedBlockInds, sortedBlockIndsData, sortedOnOffIdData);

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

    spikeArray = zeros(size(hippocampalDataReref,2), size(hippocampalDataReref,1));

    for currElect = 1:size(hippocampalDataReref)
        hippocampalDischargeContactInds = find(ismember(newContactList, resultLabel{currElect}));
        hippocampalDischargeContactInds = find(ismember(out.chan, hippocampalDischargeContactInds));
        spikeTimes = out.pos(hippocampalDischargeContactInds);

        % now we make an array that contains any spikes within the current hippocampal
        % channel, spike times are anything within +- 1 second of the spike
        spikeWindows = [];
        for j = 1:length(spikeTimes)
            currSpikeTime = spikeTimes(j);
            addTimes = spikeTimes(j) - 1:sampleFreq:spikeTimes(j) + 1;
            spikeWindows = [spikeWindows, addTimes];
        end
        spikeWindows = unique(round(spikeWindows,4));
        [~, matchInd] = intersect(round(fullDescriptionMat(:,5), 4), round(spikeWindows, 4));
        if length(matchInd) ~= length(spikeWindows)
            error('weird');
        end
        spikeArray(matchInd, currElect) = 1;
    end

    %%
    save([outdir '/' subjects{currSubject} '/timeseries_descriptor.mat'], 'fullDescriptionMat', 'fullDescriptionMatLabels', 'res', 'resultLabel', 'spikeArray', '-v7.3');
    toc;
end
