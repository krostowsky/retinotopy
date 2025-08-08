function extractHpcData(indir, subjId)
indir2 = indir;
addpath(genpath(pwd));
load([indir2 '/rereferencedDataSpikeRemoved.mat' ], '-regexp', '^(?!outdir$|indir$|subjId$)...');

addpath(genpath('/scratch/midway3/rostowsky/'));
addpath(genpath('/project/joelvoss/tmp-rostowsky'));

%%
resultLabel = matchSubjectToHpcElectrodes(subjId);

%% extract data for hippocampal and nonhippocampal
resultLabelChannels = find(ismember(newContactList, resultLabel));
nonResultLabelChannels = find(~ismember(newContactList, resultLabel));
hippocampalDataReref = spikeRemovedRereferencedData;
nonHippocampalDataReref = spikeRemovedRereferencedData;

hippocampalDataReref = spikeRemovedRereferencedData(resultLabelChannels, :);
nonHippocampalDataReref = spikeRemovedRereferencedData(nonResultLabelChannels, :);

%% hippocampal on vs off blocks
[hippocampalDataReref, hippocampalMaskInd, hippocampalImageInds, hippocampalTimes, hippocampalStimTimes, hippocampalMaskIndData, hippocampalImageIndData, sortedOnOffIdData, sortedOnOffIndsData, rawOnOffIdData, rawOnOffIndsData, sortedBlockIndsData] = ...
    extractData(hippocampalDataReref, outStruct.time, builtTimes, outStruct.masks, outStruct.images, sortedOnOffId, sortedOnOffInds, rawOnOffId, rawOnOffInds, sortedBlockInds);

%%
save([ indir2 '/' 'hpcElectrodeData.mat']);
%delete([indir2 '/' 'rereferencedDataSpikeRemoved.mat']);
fprintf('\nfinished\n');
end

%%
function [hippocampalDataReref, hippocampalCurrCondMaskInd, hippocampalCurrCondImageInds, hippocampalCurrCondTimes, hippocampalCurrCondStimTimes, hippocampalCurrCondMaskIndData, hippocampalCurrCondImageIndData, sortedOnOffIdData, sortedOnOffIndsData, rawOnOffIdData, rawOnOffIndsData, sortedBlockIndsData] = ...
    extractData(hippocampalDataReref, hippocampalCurrCondTimes, hippocampalCurrCondStimTimes, hippocampalCurrCondMaskInd, hippocampalCurrCondImageInds, sortedOnOffId, sortedOnOffInds, rawOnOffId, rawOnOffInds, sortedBlockInds)

%
refreshRateStim = diff(hippocampalCurrCondStimTimes);
refreshRateStim = refreshRateStim(1);
hippocampalCurrCondStimTimes = 0:refreshRateStim:refreshRateStim*(length(hippocampalCurrCondStimTimes)-1);

refreshRateElectrode = diff(hippocampalCurrCondTimes);
refreshRateElectrode = refreshRateElectrode(1);
hippocampalCurrCondTimes = 0:refreshRateElectrode:refreshRateElectrode*(length(hippocampalCurrCondTimes)-1);

% initialize arrays to hold "expanded" mask and image data (creates timesseries of mask and image indices that are equal in length to sample data)
hippocampalCurrCondMaskIndData = zeros(1,length(hippocampalCurrCondTimes));
hippocampalCurrCondImageIndData = zeros(1,length(hippocampalCurrCondTimes));
sortedOnOffIdData = zeros(1,length(hippocampalCurrCondTimes));
sortedOnOffIndsData = zeros(1,length(hippocampalCurrCondTimes));
rawOnOffIdData = zeros(1,length(hippocampalCurrCondTimes));
rawOnOffIndsData = zeros(1,length(hippocampalCurrCondTimes));
sortedBlockIndsData = zeros(1,length(hippocampalCurrCondTimes));

% for each data point that the stimulus is on
for j = 1:length(hippocampalCurrCondTimes)

    % find the closest stimulus time to the current data time
    if find(hippocampalCurrCondTimes(j) == hippocampalCurrCondStimTimes)
        hippocampalCurrCondMaskIndData(j) = hippocampalCurrCondMaskInd(find(hippocampalCurrCondTimes(j) == hippocampalCurrCondStimTimes));
        hippocampalCurrCondImageIndData(j) = hippocampalCurrCondImageInds(find(hippocampalCurrCondTimes(j) == hippocampalCurrCondStimTimes));
        sortedOnOffIdData(j) = sortedOnOffId(find(hippocampalCurrCondTimes(j) == hippocampalCurrCondStimTimes));
        sortedOnOffIndsData(j) = sortedOnOffInds(find(hippocampalCurrCondTimes(j) == hippocampalCurrCondStimTimes));
        rawOnOffIdData(j) = rawOnOffId(find(hippocampalCurrCondTimes(j) == hippocampalCurrCondStimTimes));
        rawOnOffIndsData(j) = rawOnOffInds(find(hippocampalCurrCondTimes(j) == hippocampalCurrCondStimTimes));
        sortedBlockIndsData(j) = sortedBlockInds(find(hippocampalCurrCondTimes(j) == hippocampalCurrCondStimTimes));
    else
        greatestInd =  find(hippocampalCurrCondTimes(j) >= hippocampalCurrCondStimTimes);
        smallestInd = find(hippocampalCurrCondTimes(j) <= hippocampalCurrCondStimTimes);
        if hippocampalCurrCondTimes(j) > hippocampalCurrCondStimTimes(greatestInd(end)) && hippocampalCurrCondTimes(j) < hippocampalCurrCondStimTimes(smallestInd(1))
            hippocampalCurrCondMaskIndData(j) = hippocampalCurrCondMaskInd(smallestInd(1));
            hippocampalCurrCondImageIndData(j) = hippocampalCurrCondImageInds(smallestInd(1));
            sortedOnOffIdData(j) = sortedOnOffId(smallestInd(1));
            sortedOnOffIndsData(j) = sortedOnOffInds(smallestInd(1));
            rawOnOffIdData(j) = rawOnOffId(smallestInd(1));
            rawOnOffIndsData(j) = rawOnOffInds(smallestInd(1));
            sortedBlockIndsData(j) = sortedBlockInds(smallestInd(1));
        else
            error('ERROR');
        end
    end
end

end
