clc; clear;
addpath(genpath('/project/joelvoss/tmp-rostowsky/github/retinotopy/dependencies/tfAnalysis/bosc/'));
addpath(genpath('/project/joelvoss/tmp-rostowsky/github/retinotopy'));

%%
indir = '/project/joelvoss/tmp-rostowsky/hpcData';
subjects = struct2cell(dir(indir));
subjects = subjects(1,3:end);
inFile = 'timeseries_descriptor.mat';

subjectData = struct();
subjectData.('subject') = cell(length(subjects),1);

figDir = '/project/joelvoss/tmp-rostowsky/hpcDataFigs';

computeAPPerEvent = 1;

% I extend the bounds a bit to include peaks that would normally be
% excluded
lowTheta = [1 5];
highTheta = [5 11];

ltapfreqs = [1:19];
htapfreqs = [20:30];

gamma = [80, 200];
gammaFreqs = [67:80]; % these are indices from the TFA, assumes tfa from 1 to 200Hz, 80 freqs, log-spaced

fs = 500;
% apFreqsInds = 17:57;
apFreqsInds = 7:57;
spikeProportionThreshold = 0.6;

removeSpikes = 1; 

%%
subjectPVals = cell(length(subjects),1);

for j = 1:length(subjects)
    load([indir '/' subjects{j} '/' inFile]);
    load([indir '/' subjects{j} '/tfAnalysis.mat'], 'hippocampalDataReref');

    electrodeOnPeriodicData = cell(size(res.B,2), 2);
    electrodeOffPeriodicData = cell(size(res.B,2), 2);

    modOnStuff = cell(size(res.B, 3),1);
    modOffStuff = cell(size(res.B, 3),1);

    % for every electrode, compare on power to average off power for each event
    apFreqs = res.freqs(apFreqsInds);

    currResults = cell(size(res.B,2), 2);

    for k = 1:size(res.B, 3)
        if removeSpikes
            numSpikesOff = sum(spikeArray(find(fullDescriptionMat(:,2) == 0), k) == 1);
            numSpikesOn = sum(spikeArray(find(fullDescriptionMat(:,2) == 1), k) == 1);
            fprintf(['subject: ' num2str(j) ' electrode: ' num2str(k) ' has ' num2str(numSpikesOff) '/' num2str(size(res.B,2)) ' IED-affected events (off) \n'])
            fprintf(['subject: ' num2str(j) ' electrode: ' num2str(k) ' has ' num2str(numSpikesOn) '/' num2str(size(res.B,2)) ' IED-affected events (on) \n'])

            if (numSpikesOff + numSpikesOn) / (size(res.B,2) + size(res.B,2)) > spikeProportionThreshold
                fprintf(['too many spikes detected for electrode ' num2str(k) ', skipping...']);
                continue;
            end

            % this is a long conditional:
            % it finds the indices of "off" or "on",
            % uses those to index the spikeArray and
            % find where there are no spikes
            currData = res.B(:, find(spikeArray(:,k) == 0), k);
            modDescriptor = fullDescriptionMat(find(spikeArray(:,k) == 0), :);
        else 
            currData = res.B(:, :, k);
            modDescriptor = fullDescriptionMat;
        end
    
        lt = nanmean(log10(currData(1:24, :)),1) - nanmean(res.ap_ps(1:24, k),1);
        ht = nanmean(log10(currData(25:38, :)),1) - nanmean(res.ap_ps(25:38, k),1);
    
        autocorrColLt = nan(length(modDescriptor),1);
        autocorrColLt(1) = 0;
        autocorrColLt(2:end) = lt(1:end-1);
        [~,~,statslt] = glmfit([modDescriptor(:,2), autocorrColLt], lt);

        autocorrColHt = nan(length(modDescriptor),1);
        autocorrColHt(1) = 0;
        autocorrColHt(2:end) = ht(1:end-1);
        [~,~,statsHt] = glmfit([modDescriptor(:,2), autocorrColHt], ht);

        currResults{k,1} = statslt.p(2);
        currResults{k,2} = statsHt.p(2);

        %[acf,lags] = autocorr(statslt.resid);
        
    end
    subjectPVals{j} = currResults;
end
